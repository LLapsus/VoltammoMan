# VoltammoMan

<p align="center">
  <img src="docs/voltammoman.jpg" alt="VoltammoMan" width="600">
</p>

Batch evaluation of cyclic voltammograms. Point it at a folder of exported curves and get one
table with the peak potentials, peak currents, peak separation and reversibility of every curve,
together with a control plot of each evaluation.

- reads **NOVA (Metrohm Autolab) ASCII exports** and similar two-column text files,
- finds the anodic and cathodic peak and measures it against an automatically constructed baseline,
- estimates the **IpA/IpC ratio by Nicholson's method**, which is more reliable than a straight
  baseline for the second peak,
- flags curves that need a manual check (broad or weak peaks, peaks at a vertex, no peak),
- command line (`vom`) and a simple graphical interface in the web browser (`vom gui`),
- includes a small **CV simulator** that writes curves in the NOVA format (used for the tests and
  the example data).

## Installation

Python ≥ 3.10 is required.

```bash
pip install "voltammoman[gui] @ git+https://github.com/LLapsus/VoltammoMan.git"
```

or from a local copy:

```bash
git clone https://github.com/LLapsus/VoltammoMan.git
cd VoltammoMan
python -m venv .venv
source .venv/bin/activate                # Windows: .venv\Scripts\activate
pip install -e ".[gui]"
```

The `vom` command is installed into the virtual environment, so activate it (the `activate` line
above) in every new terminal before using `vom`. Without `[gui]` everything except the graphical
interface works.

## Quick start

```bash
vom examples demo              # writes a few synthetic curves into ./demo
vom analyze demo --plots       # -> voltammoman_results.xlsx + voltammoman_plots/*.png
```

In a folder full of curves just run

```bash
vom analyze
```

which evaluates all `*.txt` files in the current folder and writes `voltammoman_results.xlsx`.

Graphical interface (opens in the web browser, runs only on your computer):

```bash
vom gui
```

From Python:

```python
from voltammoman import analyze_files, save_table
df = analyze_files("path/to/curves", baseline="flattest")
save_table(df, "results.xlsx")
```

A step-by-step tutorial of the Python interface (single curves, batches, flags, baseline
methods, the simulator) is in the notebook [`examples/tutorial.ipynb`](examples/tutorial.ipynb).

## Command line

| Command | What it does |
|---|---|
| `vom analyze [PATHS ...]` | evaluate files and/or folders (default: `*.txt` in the current folder) |
| `vom gui` | start the graphical interface |
| `vom examples [DIR]` | write synthetic example curves |
| `vom simulate OUT.txt` | simulate one curve (see `vom simulate -h` for parameters) |

Main options of `vom analyze`:

| Option | Meaning |
|---|---|
| `-o FILE` | output table, `.xlsx` (default `voltammoman_results.xlsx`) or `.csv` |
| `--plots [DIR]` | save one control plot per curve (default folder `voltammoman_plots`) |
| `-b {flattest,linear,tangent}` | baseline method (default `flattest`) |
| `-p PATTERN`, `-r` | file pattern in folders (default `*.txt`), search subfolders |
| `--cycle N` | cycle to evaluate in multi-cycle files (default: last complete) |
| `--edge`, `--gap`, `--width`, `--base-from`, `--base-to`, `--smooth` | peak search settings, see below |

## Input

A text file with a header line and columns of potential (V) and current (A), e.g. the NOVA export

```
Potential applied (V);Smoothed WE(1).Current (A)
0.00244140625;-3.44506730248918E-06
...
```

Separators `;`, `,` or tab and a decimal comma are recognised; columns are found by names
containing "potential" and "current" (otherwise the first two columns are used). A column named
"Scan" is used to select the cycle. The curve should contain at least one complete cycle (two
changes of the scan direction); it may start at a vertex or inside the window, in either direction.

## Output

One row per file (the XLSX file also has a sheet describing the columns). Potentials in V,
currents in A.

| Column | Meaning |
|---|---|
| `EpA`, `EpC` | anodic / cathodic peak potential |
| `dEp` | peak separation EpA − EpC |
| `E_half` | (EpA + EpC) / 2 |
| `IpA`, `IpC` | peak current above / below the baseline (IpC < 0) |
| `ratio` | IpA / \|IpC\| from the baseline-corrected currents |
| `ratio_nicholson` | IpA / \|IpC\| by Nicholson's formula (recommended for reversibility) |
| `wA`, `wC` | half-peak widths \|Ep − Ep/2\| |
| `IpA_raw`, `IpC_raw` | current at the peak without baseline |
| `first_peak` | which peak was recorded first (the forward peak) |
| `E_min`, `E_max`, `n_points` | potential window and number of points of the evaluated cycle |
| `baseline` | baseline method used |
| `flags` | problems, `A:`/`C:` for the anodic/cathodic peak: `broad_peak`, `weak_peak`, `no_peak`, `peak_at_vertex`, `nonpositive_height`, `no_baseline` |
| `error` | why a file could not be evaluated (the batch continues) |

**Always look at the control plots of flagged curves.**

## How it works

1. **Cycle and branches.** The (last complete) cycle is split at its vertices into the anodic
   (increasing potential) and cathodic branch. If the cycle starts inside the window, the
   interrupted branch is joined from the end and the beginning of the file (a steady-state cycle).
2. **Peak.** On each branch, the most prominent local extremum at least `edge` V (default 0.03)
   from the vertices is the peak candidate.
3. **Baseline.** A straight line is fitted to the background before the peak and extrapolated
   under it; the peak is then refined to the extremum of the baseline-corrected current within
   ±20 mV and the baseline is rebuilt. Methods:
   - `flattest` (default): the straightest window of width `width` (0.1 V) in the region ending
     `gap` (0.1 V) before the peak, judged by the curvature of a parabola fitted to it; among
     windows that are equally straight within the noise, the one closest to the peak,
   - `linear`: a fixed window from `base_from` (0.35 V) to `base_to` (0.15 V) before the peak,
   - `tangent`: the tangent at the point where the curve rises least steeply towards the peak.
4. **Ratio.** `ratio` divides the two baseline-corrected peak currents. For the peak recorded
   second, the true baseline is the decaying current of the first peak, not a straight line, so
   this ratio drifts away from 1 when the peaks are broad (slow electron transfer).
   `ratio_nicholson` uses Nicholson's formula (R. S. Nicholson, *Anal. Chem.* 38 (1966) 1406),

   i_rev / i_fwd = (i_rev)₀ / i_fwd + 0.485 (i_sp)₀ / i_fwd + 0.086,

   where the reverse peak current (i_rev)₀ and the current at the switching potential (i_sp)₀ are
   measured from the current level of the respective branch at the start of the forward sweep,
   which removes the charging current and offsets. It is reported as IpA/|IpC| regardless of the
   scan direction.

## Accuracy

On simulated curves (fast to slow electron transfer, with charging current, sloped background,
offset and noise) the default method finds the same peak potentials as without the background
(within a few mV) and peak currents that typically change by only a few per cent. For a chemically
stable couple `ratio_nicholson` stays within 0.9–1.1, while `ratio` drops to ≈ 0.8 for slow
electron transfer (ΔEp ≈ 200 mV). On ≈ 500 measured
curves of 3D-printed carbon electrodes it agrees with manual evaluation in NOVA within a median of
2.4 mV (peak potentials); the median peak currents differ by 2–4 %.

## Limitations

- One redox couple per curve; multiple peaks or shoulders are not resolved (the most prominent
  peak is taken).
- Linear baselines only; very broad peaks (e.g. poorly activated electrodes) have an ill-defined
  peak current and are flagged `broad_peak`. For slow electron transfer (ΔEp ≳ 300 mV) the wave
  starts so early that every linear baseline is uncertain; peak currents can then differ by
  10–40 % between methods (`flattest` and `linear` tend to give lower values than `tangent`).
- Nicholson's formula was derived for reversible waves; for slow electron transfer it is an
  approximation (±10 %).
- No iR correction, no scan-rate studies (Randles–Ševčík, Laviron) and no k⁰ estimation yet.
- The scan rate is not needed (and not read).

## Simulator

`voltammoman.simulate.simulate_cv` simulates a one-electron couple with Butler–Volmer kinetics
and planar diffusion (implicit finite differences), double-layer charging, uncompensated
resistance, a linear background, offset and noise. Its default potential program reproduces NOVA
exports (2.44 mV steps, 0 V → first vertex → second vertex → 0 V). It reproduces the
Randles–Ševčík peak current and Nicholson's ΔEp(ψ) table within a few mV.

```bash
vom simulate slow.txt --k0 1e-3 --Cdl 300 --noise 0.03 --seed 1
```

## Development

```bash
pip install -e ".[gui,test]"
pytest
```

To run the tutorial notebook: `pip install jupyter`, then `jupyter lab examples/tutorial.ipynb`.

## License

MIT, see [LICENSE](LICENSE).
