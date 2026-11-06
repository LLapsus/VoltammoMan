"""Graphical interface (Gradio) running locally in the web browser: ``vom gui``."""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from .analysis import RESULT_COLUMNS
from .batch import analyze_file, find_files, save_table
from .peaks import BASELINES, PeakSettings

SHOWN = ["file", "EpA", "EpC", "dEp", "IpA", "IpC", "ratio", "ratio_nicholson", "flags", "error"]


def _evaluate(uploads, folder, pattern, recursive, baseline, edge, gap, width, base_from, base_to, use_examples):
    files = []
    if use_examples:
        from .examples import write_examples
        files += write_examples(Path(tempfile.mkdtemp(prefix="vom_examples_")))
    if uploads:
        files += [Path(u) for u in uploads]
    if folder and folder.strip():
        files += find_files(folder.strip(), pattern or "*.txt", recursive)
    if not files:
        raise ValueError("Choose files, a folder or the example curves first.")
    settings = PeakSettings(edge=edge, gap=gap, width=width, base_from=base_from, base_to=base_to)
    rows, results = [], {}
    for f in files:
        try:
            res = analyze_file(f, baseline, settings)
            rows.append(res.values)
            results[f.name] = res
        except Exception as exc:  # noqa: BLE001
            row = dict.fromkeys(RESULT_COLUMNS, np.nan)
            row.update(file=f.name, baseline=baseline, flags="", error=f"{type(exc).__name__}: {exc}")
            rows.append(row)
    df = pd.DataFrame(rows, columns=RESULT_COLUMNS)
    out = Path(tempfile.mkdtemp(prefix="vom_results_")) / "voltammoman_results.xlsx"
    save_table(df, out)
    shown = df[SHOWN].copy()
    for c in ("EpA", "EpC", "dEp"):
        shown[c] = (shown[c] * 1e3).round(1)          # mV
    for c in ("IpA", "IpC"):
        shown[c] = (shown[c] * 1e6).round(3)          # µA
    shown[["ratio", "ratio_nicholson"]] = shown[["ratio", "ratio_nicholson"]].round(3)
    shown = shown.rename(columns={"EpA": "EpA [mV]", "EpC": "EpC [mV]", "dEp": "ΔEp [mV]",
                                  "IpA": "IpA [µA]", "IpC": "IpC [µA]", "ratio": "IpA/IpC",
                                  "ratio_nicholson": "IpA/IpC Nicholson"})
    import gradio as gr
    names = list(results)
    return results, shown, str(out), gr.update(choices=names, value=names[0] if names else None)


def _plot(results, name):
    import matplotlib
    matplotlib.use("Agg", force=False)
    import matplotlib.pyplot as plt

    from .plot import plot_result
    plt.close("all")
    if not results or name not in results:
        return None
    return plot_result(results[name])


def build_app():
    import gradio as gr

    d = PeakSettings()
    with gr.Blocks(title="VoltammoMan") as app:
        gr.Markdown("# VoltammoMan\nBatch evaluation of cyclic voltammograms (NOVA ASCII exports). "
                    "Everything runs on this computer; no data are sent anywhere.")
        state = gr.State({})
        with gr.Row():
            with gr.Column(scale=1):
                uploads = gr.File(label="Curve files", file_count="multiple", type="filepath")
                folder = gr.Textbox(label="…or a folder on this computer", placeholder="/path/to/curves")
                with gr.Row():
                    pattern = gr.Textbox(label="File pattern", value="*.txt")
                    recursive = gr.Checkbox(label="Include subfolders", value=False)
                use_examples = gr.Checkbox(label="Add synthetic example curves", value=False)
                baseline = gr.Dropdown(list(BASELINES), value="flattest", label="Baseline method")
                with gr.Accordion("Peak search settings (V)", open=False):
                    edge = gr.Number(d.edge, label="edge: reject peaks this close to a vertex")
                    gap = gr.Number(d.gap, label="gap: baseline region ends this far before the peak")
                    width = gr.Number(d.width, label="width: baseline window (flattest)")
                    base_from = gr.Number(d.base_from, label="base_from (linear)")
                    base_to = gr.Number(d.base_to, label="base_to (linear)")
                run = gr.Button("Evaluate", variant="primary")
            with gr.Column(scale=2):
                table = gr.Dataframe(label="Results (potentials in mV, currents in µA)", interactive=False, wrap=True)
                download = gr.File(label="Download results (XLSX, values in V and A)")
                pick = gr.Dropdown(label="Show curve", choices=[], interactive=True)
                plot = gr.Plot(label="Control plot")

        run.click(_evaluate,
                  [uploads, folder, pattern, recursive, baseline, edge, gap, width, base_from, base_to, use_examples],
                  [state, table, download, pick])
        pick.change(_plot, [state, pick], plot)
    return app


def launch(port=None, share=False, inbrowser=True):
    """Start the interface; blocks until interrupted."""
    try:
        import gradio  # noqa: F401
    except ImportError as exc:
        raise SystemExit("The GUI needs gradio: pip install 'voltammoman[gui]'") from exc
    build_app().launch(server_port=port, share=share, inbrowser=inbrowser)
