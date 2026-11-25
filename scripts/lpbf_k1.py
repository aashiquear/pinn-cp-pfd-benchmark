#!/usr/bin/env python3
"""Part V, case K1: AM Bench 2022 single tracks on alloy 718 (melt-pool cross section).

usage: scripts/lpbf_k1.py <study deck> targets | decks | run [case ...] | collect | figure | fit

  targets  read the NIST cross-section workbook (weaver 2022 trackdata) and write the
           per-case table (laser settings; mean and sample standard deviation of depth
           and width over the six sections of the three single tracks of a case)
  decks    write one standalone PINNACLE deck per case from the template
  run      run the decks (all, or the named cases) one after another, log the wall time
  collect  read every case's `process_state.npz` and write the comparison table
  figure   depth and width, computed against measured
  fit      one Gauss-Newton step of (absorptivity, penetration depth) on the fit cases

The cross section of a track is taken from the peak temperature each element reached:
the isotherm of the liquidus (and of the solidus, as the bracket) along the depth below
the track axis and across the top face, linearly interpolated between element centers,
averaged over the middle part of the track.

Runs need PINNACLE importable (PYTHONPATH or an installed package); `PINNACLE_PYTHON`
selects the interpreter (default: python3.11 when present, else the current one).
"""

from __future__ import annotations

import csv
import html
import json
import os
import re
import shutil
import subprocess
import sys
import time
import zipfile
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
PINNACLE_SRC = os.environ.get("PINNACLE_SRC", "/mnt/extended/projects/pinnacle/src")
SOLIDUS, LIQUIDUS = 1528.15, 1610.15
DOMAIN_UM = (700.0, 300.0, 350.0)
TRACK_UM = 500.0
DIFFUSIVITY = 33.13 / (8190.0 * 775.0)  # liquid at 1800 C: the largest in the tables
K_LIQUID = [(1610.15, 26.80), (1673.15, 27.70), (1773.15, 29.10), (1873.15, 30.45),
            (1973.15, 31.80), (2073.15, 33.13)]


def read_deck(path: Path) -> dict:
    sys.path.insert(0, PINNACLE_SRC)
    from pinnacle.io.i_parser import InputParser
    from pinnacle.io.translator import Translator

    cfg = Translator().translate(InputParser().parse(path.read_text()))
    for key in ("study", "cases"):
        if key not in cfg:
            raise SystemExit(f"{path}: missing [{key.capitalize()}] block")
    return cfg


def fit_cases(study: dict) -> list[str]:
    """Names of the fitted cases (the deck value arrives as a string or as a token list)."""
    raw = study.get("fit_cases", "")
    return [str(v) for v in raw] if isinstance(raw, (list, tuple)) else str(raw).split()


def resolve(p: str, deck: Path) -> Path:
    q = Path(str(p))
    return q if q.is_absolute() else (deck.parent / q).resolve()


def read_xlsx(path: Path) -> list[list[str]]:
    """Rows of the first sheet of an .xlsx file (standard library only)."""
    with zipfile.ZipFile(path) as z:
        shared = []
        if "xl/sharedStrings.xml" in z.namelist():
            text = z.read("xl/sharedStrings.xml").decode("utf8")
            for item in re.findall(r"<si>(.*?)</si>", text, re.S):
                shared.append(html.unescape("".join(re.findall(r"<t[^>]*>(.*?)</t>", item, re.S))))
        sheet = z.read("xl/worksheets/sheet1.xml").decode("utf8")
    rows = []
    for row in re.findall(r"<row [^>]*>(.*?)</row>", sheet, re.S):
        cells: dict[str, str] = {}
        for attrs, body in re.findall(r"<c ([^>]*?)(?:/>|>(.*?)</c>)", row, re.S):
            value = re.search(r"<v>(.*?)</v>", body or "")
            if not value:
                continue
            column = re.search(r'r="([A-Z]+)\d+"', attrs).group(1)
            cells[column] = shared[int(value.group(1))] if 't="s"' in attrs else value.group(1)
        rows.append([cells.get(c, "") for c in "ABCDEFGHIJK"])
    return rows


def cmd_targets(deck: Path, cfg: dict) -> None:
    study = cfg["study"]
    rows = read_xlsx(Path(str(study["workbook"])))
    header, body = rows[0], rows[1:]
    print("workbook columns:", header)
    groups: dict[str, list[tuple[float, float, float, float, float]]] = {}
    for r in body:
        match = re.fullmatch(r"Line ([0-9.]+)_(\d)", r[3].strip())
        if not match or r[9].strip() not in ("N/A", ""):
            continue  # pads and overlap rows carry overlap measures
        groups.setdefault(match.group(1), []).append(
            (float(r[5]), float(r[4]), float(r[6]), float(r[7]), float(r[8]))
        )
    out = resolve(study["targets"], deck)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["# AM Bench 2022 single tracks on alloy 718 (weaver 2022 trackdata, NIST mds2-2718):"
                    " mean and sample SD over the sections of the three tracks of a case; derived"])
        w.writerow(["line", "power_W", "speed_mm_s", "d4sigma_um", "sections",
                    "depth_um", "depth_sd_um", "width_um", "width_sd_um"])
        for line, vals in sorted(groups.items(), key=lambda kv: float(kv[0])):
            a = np.asarray(vals)
            assert np.all(a[:, :3] == a[0, :3]), line
            w.writerow([line, a[0, 0], a[0, 1], a[0, 2], len(vals),
                        f"{a[:, 3].mean():.2f}", f"{a[:, 3].std(ddof=1):.2f}",
                        f"{a[:, 4].mean():.2f}", f"{a[:, 4].std(ddof=1):.2f}"])
            print(f"line {line}: P {a[0,0]:.0f} W, v {a[0,1]:.0f} mm/s, D {a[0,2]:.0f} um, "
                  f"n {len(vals)}, depth {a[:,3].mean():.1f} +/- {a[:,3].std(ddof=1):.1f}, "
                  f"width {a[:,4].mean():.1f} +/- {a[:,4].std(ddof=1):.1f}")
    print("wrote", out)


def targets(deck: Path, cfg: dict) -> dict[str, dict[str, float]]:
    path = resolve(cfg["study"]["targets"], deck)
    out = {}
    with path.open() as fh:
        rows = [r for r in csv.reader(fh) if r and not r[0].startswith("#")]
    head = rows[0]
    for r in rows[1:]:
        out[r[0]] = {k: float(v) for k, v in zip(head[1:], r[1:])}
    return out


def case_table(deck: Path, cfg: dict) -> list[dict]:
    meas = targets(deck, cfg)
    cases = []
    for name, entry in cfg["cases"].items():
        line = str(entry["line"])
        line = line[:-2] if line.endswith(".0") else line
        if line not in meas:
            raise SystemExit(f"case {name}: line {line} is not in the targets table {sorted(meas)}")
        cases.append({"name": name, "line": line, **meas[line]})
    return cases


def paths(deck: Path, cfg: dict, name: str) -> tuple[Path, Path]:
    study = cfg["study"]
    tag = str(study["tag"])
    case_deck = resolve(study["case_dir"], deck) / tag / f"{name}.i"
    out_dir = resolve(study["out_dir"], deck) / tag / name
    return case_deck, out_dir


def cmd_decks(deck: Path, cfg: dict, overrides: dict | None = None, tag: str | None = None) -> None:
    study = {**cfg["study"], **(overrides or {})}
    template = resolve(study["template"], deck).read_text()
    h = float(study["mesh_um"]) * 1e-6
    n = [int(round(d * 1e-6 / h)) for d in DOMAIN_UM]
    dt = 0.8 * h * h / (6.0 * DIFFUSIVITY * float(study.get("liquid_conductivity_factor", 1.0)))
    factor = float(study.get("liquid_conductivity_factor", 1.0))
    k_liquid = " ".join(f"{t:.2f} {k * factor:.3f}" for t, k in K_LIQUID)
    for case in case_table(deck, cfg):
        case_deck, out_dir = paths(deck, {**cfg, "study": {**study, "tag": tag or study["tag"]}}, case["name"])
        case_deck.parent.mkdir(parents=True, exist_ok=True)
        fields = {
            "BACKEND": str(study["backend"]),
            "USE_GPU": "true" if study["use_gpu"] else "false",
            "OUT_DIR": os.path.relpath(out_dir, case_deck.parent),
            "NX": n[0], "NY": n[1], "NZ": n[2],
            "POWER": f"{case['power_W']:.1f}",
            "SPEED": f"{case['speed_mm_s'] * 1e-3:.4f}",
            "RADIUS": f"{0.5 * case['d4sigma_um'] * 1e-6:.4e}",
            "ABSORPTIVITY": f"{float(study['absorptivity']):.5f}",
            "PENETRATION": f"{float(study['penetration_um']) * 1e-6:.5e}",
            "BAND": f"{h:.3e}",
            "DT": f"{dt:.4e}",
            "END_TIME": f"{float(study['end_time']):.4e}",
            "K_LIQUID": k_liquid,
        }
        text = template
        for key, value in fields.items():
            text = text.replace("{" + key + "}", str(value))
        left = re.findall(r"\{[A-Z_]+\}", text.split("[Global]", 1)[1])
        if left:
            raise SystemExit(f"template fields not filled: {sorted(set(left))}")
        note = (f"# GENERATED by scripts/lpbf_k1.py from {deck.name}: case {case['name']} "
                f"(NIST line {case['line']}), mesh {float(study['mesh_um']):g} um.\n")
        case_deck.write_text(note + text)
        print("wrote", case_deck.relative_to(REPO))


def python_bin() -> str:
    return os.environ.get("PINNACLE_PYTHON") or shutil.which("python3.11") or sys.executable


def run_case(case_deck: Path, out_dir: Path) -> float:
    out_dir.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env.setdefault("PYTHONPATH", PINNACLE_SRC)
    env.setdefault("TORCH_BLAS_PREFER_HIPBLASLT", "0")
    env.setdefault("TORCH_ROCM_AOTRITON_ENABLE_EXPERIMENTAL", "1")
    env.setdefault("OMP_NUM_THREADS", "8")
    tic = time.time()
    with (out_dir / "run.log").open("w") as log:
        rc = subprocess.call([python_bin(), "-W", "ignore", "-m", "pinnacle", "run", case_deck.name],
                             cwd=case_deck.parent, env=env, stdout=log, stderr=subprocess.STDOUT)
    wall = time.time() - tic
    (out_dir / "wall.json").write_text(json.dumps({"exit": rc, "wall_s": round(wall, 1)}))
    if rc:
        raise SystemExit(f"{case_deck.name}: exit {rc}, see {out_dir / 'run.log'}")
    return wall


def cmd_run(deck: Path, cfg: dict, names: list[str]) -> None:
    for case in case_table(deck, cfg):
        if names and case["name"] not in names:
            continue
        case_deck, out_dir = paths(deck, cfg, case["name"])
        wall = run_case(case_deck, out_dir)
        print(f"{case['name']}: {wall:.0f} s", flush=True)


def cross_section(record: Path, isotherm: float) -> tuple[float, float, float]:
    """Depth, width (m) of the isotherm of the peak temperature, and the peak temperature."""
    d = np.load(record)
    c, peak = d["centroid"], d["peak_temperature"]
    xs, ys, zs = (np.unique(np.round(c[:, i], 12)) for i in range(3))
    order = np.lexsort((c[:, 0], c[:, 1], c[:, 2]))
    grid = peak[order].reshape(zs.size, ys.size, xs.size)
    top = zs[-1] + 0.5 * (zs[1] - zs[0])
    middle = (xs > xs.mean() - 0.2 * TRACK_UM * 1e-6) & (xs < xs.mean() + 0.2 * TRACK_UM * 1e-6)
    depths, widths = [], []
    jc = ys.size // 2
    for ix in np.nonzero(middle)[0]:
        section = grid[:, :, ix]
        column = np.maximum(section[::-1, jc - 1], section[::-1, jc])  # from the top down
        zr = zs[::-1]
        if column[0] < isotherm:
            depths.append(0.0); widths.append(0.0); continue
        k = int(np.argmax(column < isotherm))
        zc = zr[k - 1] + (isotherm - column[k - 1]) * (zr[k] - zr[k - 1]) / (column[k] - column[k - 1])
        depths.append(top - zc)
        row = section[-1, :]
        above = np.nonzero(row >= isotherm)[0]
        j0, j1 = above[0], above[-1]
        yl = ys[j0 - 1] + (isotherm - row[j0 - 1]) * (ys[j0] - ys[j0 - 1]) / (row[j0] - row[j0 - 1])
        yr = ys[j1] + (isotherm - row[j1]) * (ys[j1 + 1] - ys[j1]) / (row[j1 + 1] - row[j1])
        widths.append(yr - yl)
    return float(np.mean(depths)), float(np.mean(widths)), float(peak.max())


def collect(deck: Path, cfg: dict) -> list[dict]:
    rows = []
    for case in case_table(deck, cfg):
        _, out_dir = paths(deck, cfg, case["name"])
        record = out_dir / "process_state.npz"
        if not record.is_file():
            continue
        depth_l, width_l, peak = cross_section(record, LIQUIDUS)
        depth_s, width_s, _ = cross_section(record, SOLIDUS)
        wall = json.loads((out_dir / "wall.json").read_text())["wall_s"] if (out_dir / "wall.json").is_file() else float("nan")
        rows.append({
            "case": case["name"], "power_W": case["power_W"], "speed_mm_s": case["speed_mm_s"],
            "d4sigma_um": case["d4sigma_um"],
            "depth_meas_um": case["depth_um"], "width_meas_um": case["width_um"],
            "depth_liquidus_um": 1e6 * depth_l, "depth_solidus_um": 1e6 * depth_s,
            "width_liquidus_um": 1e6 * width_l, "width_solidus_um": 1e6 * width_s,
            "depth_um": 0.5e6 * (depth_l + depth_s), "width_um": 0.5e6 * (width_l + width_s),
            "peak_temperature_K": peak, "wall_s": wall,
        })
    for r in rows:
        r["depth_error_pct"] = 100.0 * (r["depth_um"] / r["depth_meas_um"] - 1.0)
        r["width_error_pct"] = 100.0 * (r["width_um"] / r["width_meas_um"] - 1.0)
    return rows


def cmd_collect(deck: Path, cfg: dict) -> list[dict]:
    rows = collect(deck, cfg)
    if not rows:
        raise SystemExit("no finished case")
    study = cfg["study"]
    out = resolve(study["out_dir"], deck) / str(study["tag"]) / "k1_cross_sections.csv"
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.3f}" if isinstance(v, float) else v) for k, v in r.items()})
    fit = fit_cases(study)
    for r in rows:
        role = "fit " if r["case"] in fit else "held"
        print(f"{r['case']:8s} {role} depth {r['depth_um']:6.1f} / {r['depth_meas_um']:6.1f} um "
              f"({r['depth_error_pct']:+5.1f} %)  width {r['width_um']:6.1f} / {r['width_meas_um']:6.1f} um "
              f"({r['width_error_pct']:+5.1f} %)  peak {r['peak_temperature_K']:.0f} K  {r['wall_s']:.0f} s")
    err = np.array([[r["depth_error_pct"], r["width_error_pct"]] for r in rows])
    print(f"rms error over {len(rows)} cases: depth {np.sqrt((err[:,0]**2).mean()):.1f} %, "
          f"width {np.sqrt((err[:,1]**2).mean()):.1f} %")
    print("wrote", out)
    return rows


def cmd_figure(deck: Path, cfg: dict) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = collect(deck, cfg)
    study = cfg["study"]
    fit = fit_cases(study)
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.0))
    for ax, key, label in zip(axes, ("depth", "width"), ("Melt-pool depth", "Melt-pool width")):
        meas = np.array([r[f"{key}_meas_um"] for r in rows])
        sim = np.array([r[f"{key}_um"] for r in rows])
        lo = np.array([min(r[f"{key}_liquidus_um"], r[f"{key}_solidus_um"]) for r in rows])
        hi = np.array([max(r[f"{key}_liquidus_um"], r[f"{key}_solidus_um"]) for r in rows])
        lim = [0.0, 1.15 * max(meas.max(), hi.max())]
        ax.plot(lim, lim, color="0.3", lw=0.8)
        ax.plot(lim, [0.9 * v for v in lim], color="0.6", lw=0.6, ls="--")
        ax.plot(lim, [1.1 * v for v in lim], color="0.6", lw=0.6, ls="--")
        for r, m, s, a, b in zip(rows, meas, sim, lo, hi):
            fitted = r["case"] in fit
            ax.errorbar(m, s, yerr=[[s - a], [b - s]], fmt="s" if fitted else "o", ms=6,
                        color="C3" if fitted else "C0", capsize=2)
            ax.annotate(r["case"].replace("case", ""), (m, s), textcoords="offset points",
                        xytext=(5, -9), fontsize=7)
        ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect("equal")
        ax.set_xlabel(f"{label}, measured ($\\mu$m)"); ax.set_ylabel(f"{label}, computed ($\\mu$m)")
    axes[0].plot([], [], "s", color="C3", label="fitted cases")
    axes[0].plot([], [], "o", color="C0", label="held-out cases")
    axes[0].legend(frameon=False, fontsize=8, loc="upper left")
    fig.tight_layout()
    out = REPO / "figures" / "lpbf" / f"k1_cross_sections_{study['tag']}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200)
    print("wrote", out.relative_to(REPO))


def cmd_fit(deck: Path, cfg: dict) -> None:
    """One Gauss-Newton step of (absorptivity, penetration depth) on the fit cases.

    Residuals: relative errors of depth and width of the fit cases. The Jacobian is by
    forward differences (two more runs per fit case); the runs go to tags <tag>_fit_*.
    """
    study = cfg["study"]
    fit = fit_cases(study)
    base = {"absorptivity": float(study["absorptivity"]), "penetration_um": float(study["penetration_um"])}
    steps = {"absorptivity": 0.05, "penetration_um": 8.0}

    def residual(tag: str, overrides: dict) -> np.ndarray:
        local = {**cfg, "study": {**study, **overrides, "tag": tag}}
        cmd_decks(deck, local)
        cmd_run(deck, local, fit)
        rows = {r["case"]: r for r in collect(deck, local)}
        return np.array([v for c in fit for v in (rows[c]["depth_error_pct"] / 100.0,
                                                  rows[c]["width_error_pct"] / 100.0)])

    tag = str(study["tag"])
    r0 = residual(f"{tag}_fit_base", base)
    jac = np.zeros((r0.size, 2))
    for j, key in enumerate(base):
        r1 = residual(f"{tag}_fit_{key}", {**base, key: base[key] + steps[key]})
        jac[:, j] = (r1 - r0) / steps[key]
    delta, *_ = np.linalg.lstsq(jac, -r0, rcond=None)
    new = {k: base[k] + d for k, d in zip(base, delta)}
    predicted = r0 + jac @ delta
    report = {"base": base, "residual_base": r0.tolist(), "jacobian": jac.tolist(),
              "step": dict(zip(base, delta.tolist())), "next": new,
              "rms_base_pct": float(100 * np.sqrt((r0**2).mean())),
              "rms_predicted_pct": float(100 * np.sqrt((predicted**2).mean()))}
    out = resolve(study["out_dir"], deck) / tag / "k1_fit_step.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=1))
    print(json.dumps(report, indent=1))


def main() -> None:
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    deck = Path(sys.argv[1]).resolve()
    command, args = sys.argv[2], sys.argv[3:]
    cfg = read_deck(deck)
    if command == "targets":
        cmd_targets(deck, cfg)
    elif command == "decks":
        cmd_decks(deck, cfg)
    elif command == "run":
        cmd_run(deck, cfg, args)
    elif command == "collect":
        cmd_collect(deck, cfg)
    elif command == "figure":
        cmd_figure(deck, cfg)
    elif command == "fit":
        cmd_fit(deck, cfg)
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
