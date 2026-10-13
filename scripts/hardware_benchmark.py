#!/usr/bin/env python3
"""Hardware and time benchmark of the routes on one plate deck.

Deck-driven (PINNACLE ``.i`` syntax):

  [Benchmark]
    out_dir  = ../../results/cp_study/hardware
    figures  = ../../figures/report
    prefix   = hardware_h40
    n_steps  = 16                       # steps of every short run (same load path start)
    repeats  = 1
  [Cases]
    [./fe_numpy_t8]   deck = ../cp_study/hybrid/h40_linear.i   backend = numpy  threads = 8  predictor = none  []
    [./fe_numpy_t1]   deck = ../cp_study/hybrid/h40_linear.i   backend = numpy  threads = 1  predictor = none  []
    [./fe_torch_cpu]  deck = ../cp_study/hybrid/h40_linear.i   backend = torch  use_gpu = false threads = 8 predictor = none []
    [./fe_torch_gpu]  deck = ../cp_study/hybrid/h40_linear.i   backend = torch  use_gpu = true  threads = 8 predictor = none []
    [./linear_t8]     deck = ../cp_study/hybrid/h40_linear.i   backend = numpy  threads = 8  []
    [./grain_gpu]     deck = ../cp_study/hybrid/h40_fem_solver_grain.i backend = numpy threads = 8 predictor_device = cuda []
    [./grain_cpu]     deck = ../cp_study/hybrid/h40_fem_solver_grain.i backend = numpy threads = 8 predictor_device = cpu  []
  []

Every case takes a PIPELINE deck of the study (already staged: its case
directory must hold ``cpfe.i``), rewrites the generated PINNACLE deck into
``<out_dir>/<case>/cpfe.i`` with the requested backend / device / thread
count / predictor device and ``n_steps`` (the load program is rescaled so
the strain per step is unchanged), runs it with ``OMP_NUM_THREADS`` set,
and records from ``solver_history.csv``: wall time per step, Newton time
per iteration, predictor time per step, assemblies per step, plus the peak
resident memory of the process and the peak GPU memory when a GPU is used
(``nvidia-smi`` / ``rocm-smi`` polled). ``predictor = none`` strips the
``[Utilities]`` block (the cold FE run).

Writes ``<out_dir>/<prefix>.csv`` and a bar figure of time per step by case.

Usage: hardware_benchmark.py decks/benchmark/hardware_h40.i [case ...]
"""
from __future__ import annotations

import csv
import os
import re
import resource
import subprocess
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, "/mnt/extended/projects/pinnacle/src")
STUDY = Path(__file__).resolve().parent.parent
PINNACLE = os.environ.get("PINNACLE_BIN", "pinnacle")


def read_deck(path: Path) -> dict:
    from pinnacle.io.i_parser import InputParser
    from pinnacle.io.translator import Translator
    cfg = Translator().translate(InputParser().parse(path.read_text()))
    for key in ("benchmark", "cases"):
        if key not in cfg:
            raise SystemExit(f"{path}: missing [{key.capitalize()}] block")
    return cfg


def _resolve(p: str, deck: Path) -> Path:
    q = Path(p)
    return q if q.is_absolute() else (deck.parent / q).resolve()


def _pipeline_out_dir(pipe: Path) -> Path:
    m = re.search(r"^\s*out_dir\s*=\s*(\S+)", pipe.read_text(), re.M)
    if not m:
        raise SystemExit(f"{pipe}: no out_dir in [Mechanics]")
    return _resolve(m.group(1), pipe)


def rewrite(src: Path, dst: Path, case: dict, n_steps: int) -> None:
    text = src.read_text()
    def setkey(block: str, key: str, value: str) -> None:
        nonlocal text
        m = re.search(rf"^\[{block}\]", text, re.M)  # a block header at line start, not a mention in a comment
        if not m:
            raise SystemExit(f"{src}: no [{block}] block")
        i = m.start(); j = text.index("\n[]", i)
        seg = text[i:j]
        pat = re.compile(rf"^(\s*{re.escape(key)}\s*=\s*)\S+", re.M)
        if pat.search(seg):
            seg = pat.sub(rf"\g<1>{value}", seg, count=1)
        else:
            seg = seg + f"\n  {key} = {value}"
        text = text[:i] + seg + text[j:]
    setkey("Global", "backend", str(case.get("backend", "numpy")))
    setkey("Global", "use_gpu", "true" if bool(case.get("use_gpu", False)) else "false")
    setkey("Global", "output_dir", str(dst.parent / "solver"))
    # Shorten the load path keeping the strain per step.
    m_n = re.search(r"^\s*n_steps\s*=\s*(\S+)", text, re.M); m_t = re.search(r"^\s*t_end\s*=\s*(\S+)", text, re.M)
    n_full, t_full = int(float(m_n.group(1))), float(m_t.group(1))
    setkey("Executioner", "n_steps", str(n_steps))
    setkey("Executioner", "t_end", f"{t_full * n_steps / n_full:.12g}")
    setkey("Executioner", "snapshots_every", "0")
    if str(case.get("predictor", "")).lower() == "none":
        text = re.sub(r"^\[Utilities\].*?\n\[\]\n", "", text, flags=re.S | re.M)
        setkey("Executioner", "corrector", "newton")
    elif case.get("predictor_device"):
        text = re.sub(r"(^\[Utilities\].*?^\s*device\s*=\s*)\S+", rf"\g<1>{case['predictor_device']}", text, count=1, flags=re.S | re.M)
    dst.write_text(text)


def gpu_peak_mb_poll(proc: subprocess.Popen, interval: float = 2.0) -> float:
    peak = 0.0
    while proc.poll() is None:
        for cmd, pat in (("nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits", r"(\d+)"),
                         ("rocm-smi --showmemuse --csv", r"(\d+)")):
            try:
                out = subprocess.run(cmd.split(), capture_output=True, text=True, timeout=5).stdout
                vals = [float(v) for v in re.findall(pat, out)]
                if vals:
                    peak = max(peak, max(vals) if "nvidia" in cmd else max(vals) / 1e6)
                    break
            except Exception:
                continue
        time.sleep(interval)
    return peak


def run_case(name: str, case: dict, bench: dict, deck: Path, out_dir: Path) -> dict:
    pipe = _resolve(str(case["deck"]), deck)
    src = _pipeline_out_dir(pipe) / "cpfe.i"
    if not src.exists():
        raise SystemExit(f"{name}: {src} not found (stage the pipeline deck first)")
    cdir = out_dir / name; cdir.mkdir(parents=True, exist_ok=True)
    dst = cdir / "cpfe.i"
    n_steps = int(bench.get("n_steps", 16))
    rewrite(src, dst, case, n_steps)
    # Point files and orientation inputs are referenced relative to the source directory.
    for extra in src.parent.glob("*.csv"):
        if not (cdir / extra.name).exists():
            os.symlink(extra, cdir / extra.name)
    env = dict(os.environ, OMP_NUM_THREADS=str(int(case.get("threads", 8))), MKL_NUM_THREADS=str(int(case.get("threads", 8))))
    results = []
    for rep in range(int(bench.get("repeats", 1))):
        tic = time.perf_counter()
        before = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        with (cdir / "run.log").open("w") as log:
            proc = subprocess.Popen([PINNACLE, "run", "cpfe.i"], cwd=cdir, stdout=log, stderr=subprocess.STDOUT, env=env)
            gpu_peak = gpu_peak_mb_poll(proc) if bool(case.get("use_gpu", False)) or case.get("predictor_device") not in (None, "", "cpu") else 0.0
            proc.wait()
        wall = time.perf_counter() - tic
        rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024.0
        hist = [r for r in csv.DictReader((cdir / "solver" / "solver_history.csv").open()) if str(r.get("converged", "1")).strip() in ("1", "1.0", "true", "True")] if (cdir / "solver" / "solver_history.csv").exists() else []
        it = np.array([float(r.get("iterations", 0) or 0) for r in hist]); nt = np.array([float(r.get("newton_time_s", 0) or 0) for r in hist])
        pt = np.array([float(r.get("predictor_time_s", 0) or 0) if r.get("predictor_time_s") not in ("", None, "nan") else 0.0 for r in hist])
        asm = np.array([float(r.get("predictor_assemblies", 0) or 0) if r.get("predictor_assemblies") not in ("", None, "nan") else 0.0 for r in hist])
        results.append({"case": name, "repeat": rep, "exit": proc.returncode, "backend": case.get("backend", "numpy"), "use_gpu": bool(case.get("use_gpu", False)),
                        "threads": int(case.get("threads", 8)), "predictor_device": case.get("predictor_device", ""), "steps": len(hist), "wall_s": wall,
                        "wall_per_step_s": wall / max(len(hist), 1), "newton_iterations": int(it.sum()), "newton_s_per_iteration": float(nt.sum() / max(it.sum(), 1)),
                        "newton_s": float(nt.sum()), "predictor_s": float(pt.sum()), "predictor_s_per_step": float(pt.sum() / max(len(hist), 1)),
                        "assemblies": int(asm.sum()), "steps_without_newton": int(np.sum(it == 0)), "peak_rss_mb": rss, "peak_gpu_mb": gpu_peak})
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in results[-1].items()})
    return results


def main(argv):
    if len(argv) < 2:
        print(__doc__); return 2
    deck = Path(argv[1]).resolve(); only = set(argv[2:])
    cfg = read_deck(deck); bench = cfg["benchmark"]
    out_dir = _resolve(str(bench["out_dir"]), deck); out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = _resolve(str(bench.get("figures", bench["out_dir"])), deck); fig_dir.mkdir(parents=True, exist_ok=True)
    prefix = str(bench.get("prefix", "hardware"))
    rows = []
    for name, case in cfg["cases"].items():
        if only and name not in only:
            continue
        rows += run_case(name, case, bench, deck, out_dir)
    path = out_dir / f"{prefix}.csv"
    existing = list(csv.DictReader(path.open())) if path.exists() and only else []
    keep = [r for r in existing if r["case"] not in {x["case"] for x in rows}]
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(keep); w.writerows(rows)
    allrows = keep + rows
    plt.rcParams.update({"font.size": 9, "axes.labelsize": 10, "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True})
    fig, ax = plt.subplots(figsize=(7.5, 3.6), dpi=200)
    names = [r["case"] for r in allrows]; x = np.arange(len(names))
    ax.bar(x, [float(r["newton_s"]) / max(int(r["steps"]), 1) for r in allrows], 0.6, color="#1f4e9c", label="Newton per step")
    ax.bar(x, [float(r["predictor_s_per_step"]) for r in allrows], 0.6, bottom=[float(r["newton_s"]) / max(int(r["steps"]), 1) for r in allrows], color="#d1495b", label="predictor per step")
    ax.set_xticks(x); ax.set_xticklabels(names, rotation=45, ha="right", fontsize=7); ax.set_ylabel("time per load step (s)"); ax.grid(False); ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(fig_dir / f"{prefix}_time_per_step.png"); plt.close(fig)
    print("wrote", path)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
