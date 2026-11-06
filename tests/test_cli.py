import pandas as pd
import pytest

from voltammoman.cli import main
from voltammoman.examples import EXAMPLES


def test_examples_analyze_and_plots(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert main(["examples", "ex"]) == 0
    (tmp_path / "ex" / "broken.txt").write_text("not a curve\n")
    assert main(["analyze", "ex", "-o", "res.csv", "--plots", "png", "-q"]) == 0
    df = pd.read_csv("res.csv")
    assert len(df) == len(EXAMPLES) + 1
    assert df.set_index("file").loc["broken.txt", "error"].startswith("ValueError")
    ok = df[df["file"].str.startswith("sim-ox-k0_1e-1")].iloc[0]
    assert 0.05 < ok["dEp"] < 0.09
    assert len(list((tmp_path / "png").glob("*.png"))) == len(EXAMPLES)


def test_analyze_xlsx_and_simulate(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert main(["simulate", "one.txt", "--k0", "0.01", "--Cdl", "300", "--noise", "0.03", "--seed", "1"]) == 0
    assert main(["analyze", "one.txt", "-o", "r.xlsx", "-q"]) == 0
    x = pd.read_excel("r.xlsx", sheet_name=None)
    assert set(x) == {"results", "columns"}
    assert x["results"].loc[0, "dEp"] == pytest.approx(0.081, abs=0.01)


def test_gui_builds():
    pytest.importorskip("gradio")
    from voltammoman.gui import build_app
    assert build_app() is not None
