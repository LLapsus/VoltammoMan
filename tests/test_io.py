import numpy as np
import pytest

from voltammoman.io import NOVA_HEADER, format_current, format_potential, read_curve, write_curve
from voltammoman.simulate import potential_program


def test_number_formats_match_nova():
    assert format_potential(0.0) == "0"
    assert format_potential(10 / 4096) == "0.00244140625"
    assert format_potential(0.8984375) == "0.8984375"
    assert format_potential(-0.30029296875) == "-0.30029296875"
    assert format_current(-3.44506730248918e-06) == "-3.44506730248918E-06"
    assert format_current(1.3838242039536e-05) == "1.3838242039536E-05"
    assert format_current(1.161088159074e-05) == "1.161088159074E-05"


def test_write_produces_nova_layout(tmp_path):
    E = potential_program()
    I = np.sin(E) * 1e-5
    p = write_curve(tmp_path / "c.txt", E, I)
    raw = p.read_bytes()
    lines = raw.split(b"\r\n")
    assert lines[0].decode() == ";".join(NOVA_HEADER)
    assert lines[1].startswith(b"0.00244140625;")
    assert lines[-2].startswith(b"0;")          # last point is the start potential
    assert raw.endswith(b"\r\n") and b"\n" not in raw.replace(b"\r\n", b"")
    assert len(lines) - 2 == 984                # same length as a NOVA export of this window


def test_roundtrip(tmp_path):
    E = potential_program(0, 0.3, -0.5)
    I = np.random.default_rng(0).normal(0, 1e-5, len(E))
    cv = read_curve(write_curve(tmp_path / "c.txt", E, I))
    np.testing.assert_array_equal(cv["E"], E)
    np.testing.assert_allclose(cv["I"], I, rtol=1e-13)


@pytest.mark.parametrize("sep,dec", [(";", ","), ("\t", "."), (",", ".")])
def test_read_variants(tmp_path, sep, dec):
    E = np.linspace(-0.2, 0.6, 50)
    I = E * 1e-5
    rows = [f"Time (s){sep}WE(1).Potential (V){sep}WE(1).Current (A){sep}Scan"]
    for k, (e, i) in enumerate(zip(E, I)):
        rows.append(sep.join([f"{k * 0.1:.1f}", f"{e:.6f}", f"{i:.6e}", "1"]).replace(".", dec)
                    if dec == "," else sep.join([f"{k * 0.1:.1f}", f"{e:.6f}", f"{i:.6e}", "1"]))
    p = tmp_path / "v.txt"
    p.write_text("\n".join(rows))
    cv = read_curve(p)
    np.testing.assert_allclose(cv["E"], E, atol=1e-6)
    np.testing.assert_allclose(cv["I"], I, rtol=1e-5)
    assert "scan" in cv


def test_read_rejects_garbage(tmp_path):
    p = tmp_path / "x.txt"
    p.write_text("hello\nworld\n")
    with pytest.raises(ValueError):
        read_curve(p)
