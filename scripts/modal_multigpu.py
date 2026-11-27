#!/usr/bin/env python3
"""Part V: scaling of the PINNACLE device pool over several NVIDIA GPUs, run on Modal.

usage: modal run scripts/modal_multigpu.py [--gpu A10G] [--counts 1,2,4] [--steps 200] [--mesh 5] [--no-mechanics] [--profile]

Run with the client of `.venv-modal/` (the author's token). The container asks for 8 CPU
cores: with Modal's default share of a core the host-side work of a step (services,
Python) hides the GPUs (first run of 2026-11-25). The script sends the PINNACLE source and the
decks of this repository to Modal and runs them there; nothing else leaves the machine.

What it measures, on one container with the largest count of `--counts` GPUs:

  explicit   the single track with melt-pool flow (decks/lpbf/in718/k1/cases/flow_h10/
             case0.i, mesh refined to `--mesh` um by an override), `--steps` steps, with
             `[Global] devices` set to the first 1, 2, 4, ... GPUs: seconds per step
             and the largest difference of the final temperature against one GPU
  mechanics  the tensile cube of case K3 (k3_cube_n6_seed1.i refined to 16 x 16 x 16
             elements, 32,768 quadrature points, 3 load steps) with the crystal-plasticity
             points on the same device lists: wall time

It prints one JSON document; copy it into results/lpbf/hardware/ and add a ledger row.
Cost: the container bills all its GPUs for the whole run; with the defaults (4 x A10G,
about 15 minutes) a few US dollars. Check the current prices of Modal before a larger run.
"""

from __future__ import annotations

import json
from pathlib import Path

import modal

REPO = Path(__file__).resolve().parent.parent
PINNACLE = REPO.parent / "pinnacle"

image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install("torch>=2.4", "numpy>=1.26", "scipy>=1.11", "jsonschema>=4.0", "pyyaml", "h5py")
    .add_local_dir(str(PINNACLE / "src"), "/work/pinnacle/src", copy=True)
    .add_local_dir(str(PINNACLE / "schemas"), "/work/pinnacle/schemas", copy=True)
    .add_local_dir(str(REPO / "decks" / "lpbf"), "/work/bench/decks/lpbf", copy=True)
    .add_local_dir(str(REPO / "scripts"), "/work/bench/scripts", copy=True)
    .add_local_dir(str(REPO / "data" / "lpbf" / "in718"), "/work/bench/data/lpbf/in718", copy=True)
)
app = modal.App("pinnacle-device-pool-scaling", image=image)


def _run(script: str, deck: str, overrides: list[str], out: str, profile: bool = False,
         visible: str | None = None) -> dict:
    import os
    import subprocess
    import time

    env = dict(os.environ, PYTHONPATH="/work/pinnacle/src")
    if visible is not None:
        env["CUDA_VISIBLE_DEVICES"] = visible
    tic = time.time()
    proc = subprocess.run(
        ["python", *(["-m", "cProfile", "-o", f"{out}.prof"] if profile else []),
         f"/work/bench/scripts/{script}", *deck.split(), *overrides, f"global.output_dir={out}"],
        cwd="/work/bench", env=env, capture_output=True, text=True,
    )
    wall = time.time() - tic
    report = {"wall_with_setup_s": round(wall, 1), "returncode": proc.returncode}
    wall_file = Path(out) / "wall.json"
    if wall_file.exists():
        report.update(json.loads(wall_file.read_text()))
    else:
        report["stderr_tail"] = proc.stderr[-2000:]
    if profile and Path(f"{out}.prof").exists():
        import io
        import pstats

        text = io.StringIO()
        pstats.Stats(f"{out}.prof", stream=text).sort_stats("cumulative").print_stats(70)
        report["profile_main_thread"] = text.getvalue()[-12000:]
    return report


def _measure(counts: list[int], steps: int, mesh_um: float, mechanics: bool = True,
             profile: bool = False) -> dict:
    import numpy as np
    import torch

    names = [torch.cuda.get_device_name(i) for i in range(torch.cuda.device_count())]
    report: dict = {"gpus": names, "torch": str(torch.__version__), "explicit": {}, "mechanics": {}}
    n = [int(round(v / mesh_um)) for v in (700.0, 300.0, 350.0)]
    dt = 5.2869e-07 * (mesh_um / 10.0) ** 2
    track = "decks/lpbf/in718/k1/cases/flow_h10/case0.i"
    cube = "run decks/lpbf/in718/k3/k3_cube_n6_seed1.i"
    reference = None
    for count in counts:
        devices = json.dumps(" ".join(f"cuda:{i}" for i in range(count)))
        out = f"/tmp/explicit_{count}"
        row = _run(
            "lpbf_run.py", track,
            [f"mesh.divisions={json.dumps(n)}", f"executioner.dt={dt}",
             f"executioner.end_time={dt * steps}", "executioner.snapshots_every=0",
             f"kernels.carried.mesh_size={mesh_um * 1e-6}", f"kernels.momentum.mesh_size={mesh_um * 1e-6}",
             f"kernels.surface.band_thickness={mesh_um * 1e-6}",
             f"kernels.marangoni.band_thickness={mesh_um * 1e-6}",
             f"global.devices={devices}", "global.pool_min_batch=1", "global.num_threads=8"],
            out, profile,
        )
        fields = Path(out) / "fields.npz"
        if fields.exists():
            temperature = np.load(fields)["temperature"].ravel()
            if reference is None:
                reference = temperature
            row["max_abs_difference_to_one_gpu_K"] = float(np.abs(temperature - reference).max())
        report["explicit"][str(count)] = row
        if not mechanics:
            continue
        report["mechanics"][str(count)] = _run(
            "lpbf_k3.py", cube,
            ["global.use_gpu=true", f"global.devices={devices}", "global.pool_min_batch=1",
             "mesh.divisions=[16,16,16]", "executioner.n_steps=3", "executioner.t_end=6.0",
             "bcs.pull.values=[7.2e-7]", "global.num_threads=8"],
            f"/tmp/mechanics_{count}",
        )
    return report


@app.function(gpu="A10G:4", cpu=8.0, memory=32768, timeout=3600)
def measure_a10g(counts: list[int], steps: int, mesh_um: float, mechanics: bool = True,
                 profile: bool = False) -> dict:
    return _measure(counts, steps, mesh_um, mechanics, profile)


@app.function(gpu="L4:4", cpu=8.0, memory=32768, timeout=3600)
def measure_l4(counts: list[int], steps: int, mesh_um: float, mechanics: bool = True,
                 profile: bool = False) -> dict:
    return _measure(counts, steps, mesh_um, mechanics, profile)


@app.function(gpu="H100:8", cpu=8.0, memory=32768, timeout=1800)
def measure_h100(counts: list[int], steps: int, mesh_um: float, mechanics: bool = True,
                 profile: bool = False) -> dict:
    return _measure(counts, steps, mesh_um, mechanics, profile)


@app.function(gpu="A10G:4", cpu=8.0, memory=32768, timeout=2400)
def processes(steps: list[int], mesh_um: float) -> dict:
    """What one process per card could reach: the track with flow as one process on one
    card (whole mesh), as one process on one card with a quarter of the elements, and as
    four such processes at the same time, one per card. No exchange between the four:
    the time of the last is the floor of a split into four processes that each own a
    quarter of the mesh. Two step counts separate the set-up from the time per step."""
    from concurrent.futures import ThreadPoolExecutor

    n = [int(round(v / mesh_um)) for v in (700.0, 300.0, 350.0)]
    dt = 5.2869e-07 * (mesh_um / 10.0) ** 2
    track = "decks/lpbf/in718/k1/cases/flow_h10/case0.i"

    def one(tag: str, divisions: list[int], count: int, card: str, threads: int) -> dict:
        row = _run(
            "lpbf_run.py", track,
            [f"mesh.divisions={json.dumps(divisions)}", f"executioner.dt={dt}",
             f"executioner.end_time={dt * count}", "executioner.snapshots_every=0",
             f"kernels.carried.mesh_size={mesh_um * 1e-6}", f"kernels.momentum.mesh_size={mesh_um * 1e-6}",
             f"kernels.surface.band_thickness={mesh_um * 1e-6}",
             f"kernels.marangoni.band_thickness={mesh_um * 1e-6}",
             f"global.num_threads={threads}"],
            f"/tmp/{tag}_{count}_{card}", visible=card,
        )
        return {k: row.get(k) for k in ("wall_s", "steps", "elements", "returncode", "stderr_tail")}

    quarter = [n[0] // 4, n[1], n[2]]
    report: dict = {"mesh_um": mesh_um, "divisions": n, "quarter_divisions": quarter}
    for count in steps:
        report[f"whole_one_card_{count}"] = one("whole", n, count, "0", 8)
        report[f"quarter_one_card_alone_{count}"] = one("quarter", quarter, count, "0", 2)
        with ThreadPoolExecutor(4) as ex:
            rows = list(ex.map(lambda c: one("together", quarter, count, str(c), 2), range(4)))
        report[f"four_quarters_together_{count}"] = rows
    return report


@app.function(gpu="L4:4", cpu=8.0, memory=49152, timeout=1500, retries=0)
def decomposed_l4(steps: list[int], mesh_um: float, rank_counts: list[int], partition: str) -> dict:
    return decomposed.local(steps, mesh_um, rank_counts, partition)


@app.function(gpu="A100-40GB:4", cpu=8.0, memory=49152, timeout=1200, retries=0)
def decomposed_a100(steps: list[int], mesh_um: float, rank_counts: list[int], partition: str) -> dict:
    return decomposed.local(steps, mesh_um, rank_counts, partition)


@app.function(gpu="T4:4", cpu=8.0, memory=49152, timeout=1500, retries=0)
def decomposed_t4(steps: list[int], mesh_um: float, rank_counts: list[int], partition: str) -> dict:
    return decomposed.local(steps, mesh_um, rank_counts, partition)


@app.function(gpu="A10G:4", cpu=8.0, memory=49152, timeout=1500, retries=0)
def decomposed(steps: list[int], mesh_um: float, rank_counts: list[int], partition: str) -> dict:
    """Decomposed runs (PINNACLE ADR 0057, library gpumpi): the track with flow as one
    process on one card, and with the mesh cut into parts, one process and one card per
    part, exchanging the nodal sums at the cuts. Two step counts separate the start of
    the processes and the set-up from the time per step; the fields of the longer runs
    are compared with the single run."""
    import numpy as np

    n = [int(round(v / mesh_um)) for v in (700.0, 300.0, 350.0)]
    dt = 5.2869e-07 * (mesh_um / 10.0) ** 2
    track = "decks/lpbf/in718/k1/cases/flow_h10/case0.i"

    def one(tag: str, count: int, extra: list[str]) -> tuple[dict, str]:
        import os
        import subprocess
        import time

        out = f"/tmp/{tag}_{count}"
        # under the memory guard, with its lines in the log of the container: a run that
        # eats the memory is stopped and reported in place of a killed container
        command = [
            "python", "/work/bench/scripts/memory_guard.py", "6", "600",
            "python", "/work/bench/scripts/lpbf_run.py", track,
            f"mesh.divisions={json.dumps(n)}", f"executioner.dt={dt}",
            f"executioner.end_time={dt * count}", "executioner.snapshots_every=0",
            f"kernels.carried.mesh_size={mesh_um * 1e-6}", f"kernels.momentum.mesh_size={mesh_um * 1e-6}",
            f"kernels.surface.band_thickness={mesh_um * 1e-6}",
            f"kernels.marangoni.band_thickness={mesh_um * 1e-6}", *extra, f"global.output_dir={out}",
        ]
        print(f"=== {tag} {count} steps", flush=True)
        tic = time.time()
        proc = subprocess.run(command, cwd="/work/bench", env=dict(os.environ, PYTHONPATH="/work/pinnacle/src"))
        row = {"wall_with_setup_s": round(time.time() - tic, 1), "returncode": proc.returncode}
        wall_file = Path(out) / "wall.json"
        if wall_file.exists():
            row.update(json.loads(wall_file.read_text()))
        elif Path("/tmp/guard_out.txt").exists():
            row["stderr_tail"] = Path("/tmp/guard_out.txt").read_text()[-1500:]
            print(row["stderr_tail"], flush=True)
        keep = ("wall_s", "wall_with_setup_s", "wall_with_start_and_merge_s", "steps", "elements",
                "returncode", "stderr_tail", "ranks")
        return {k: row.get(k) for k in keep if row.get(k) is not None}, out

    report: dict = {"mesh_um": mesh_um, "divisions": n, "partition": partition}
    reference = None
    for count in steps:
        row, out = one("single", count, ["global.num_threads=8"])
        report[f"single_{count}"] = row
        if count == max(steps) and (Path(out) / "fields.npz").exists():
            reference = {k: v for k, v in np.load(Path(out) / "fields.npz").items()}
        for ranks in rank_counts:
            devices = json.dumps(" ".join(f"cuda:{i}" for i in range(ranks)))
            row, out = one(f"ranks{ranks}", count,
                           [f"global.rank_devices={devices}", f"global.partition={json.dumps(partition)}",
                            f"global.num_threads={max(1, 8 // ranks)}"])
            fields = Path(out) / "fields.npz"
            if reference is not None and count == max(steps) and fields.exists():
                got = np.load(fields)
                row["nodes_equal"] = bool(np.array_equal(got["nodes"], reference["nodes"]))
                row["shapes_equal"] = all(got[k].shape == reference[k].shape for k in reference)
                n_nodes = reference["nodes"].shape[0]  # same shape before subtracting: no broadcast
                row["max_abs_difference_to_single"] = {
                    k: float(np.abs(got[k].astype(float).reshape(n_nodes, -1)
                                    - reference[k].astype(float).reshape(n_nodes, -1)).max())
                    for k in ("temperature", "velocity", "pressure")
                }
                row["field_maxima"] = {k: float(np.abs(reference[k]).max()) for k in ("temperature", "velocity", "pressure")}
            report[f"ranks{ranks}_{count}"] = row
    return report


@app.function(gpu="A100-40GB:4", cpu=8.0, memory=65536, timeout=1800, retries=0)
def chain_parts(steps: list[int], cube: list[int]) -> dict:
    """The grain field and the stress stage on parts (PINNACLE ADR 0057, addendum).

    grains  the wide build column with flow and 30 grain fields at 10 um (518,400
            elements), as one process on one card and cut into four parts on four cards,
            for each step count of `steps`: time per step, device memory, differences
    stress  the polycrystal tensile cube (two load steps) at 16^3 elements as one process
            (direct solve) and on four cards (GMRES over the parts), and at the `cube`
            sizes on four cards only: wall time, linear iterations, differences
    Every run goes through scripts/memory_guard.py, with its lines in the log."""
    import os
    import subprocess
    import time

    cards = "cuda:0 cuda:1 cuda:2 cuda:3"
    env = dict(os.environ, PYTHONPATH="/work/pinnacle/src")

    def check(tag: str, limit: int, arguments: list[str]) -> dict:
        print(f"=== {tag}", flush=True)
        tic = time.time()
        proc = subprocess.run(
            ["python", "/work/bench/scripts/memory_guard.py", "6", str(limit),
             "python", "/work/bench/scripts/lpbf_parts_check.py", *arguments],
            cwd="/work/bench", env=env,
        )
        text = Path("/tmp/guard_out.txt").read_text() if Path("/tmp/guard_out.txt").exists() else ""
        rows = [line for line in text.splitlines() if line.startswith("PARTS-CHECK ")]
        row = json.loads(rows[-1][len("PARTS-CHECK "):]) if rows else {"output_tail": text[-1500:]}
        row["returncode"] = proc.returncode
        row["wall_of_the_check_s"] = round(time.time() - tic, 1)
        if not rows:
            print(text[-1500:], flush=True)
        return row

    report: dict = {"cards": cards}
    dt = 5.0e-7
    column = "decks/lpbf/in718/k2/k2_column_flow_wide_h10.i"
    for count in steps:
        report[f"grains_{count}"] = check(f"grains {count} steps", 700, [
            column, cards, f"executioner.end_time={dt * count}", "executioner.snapshots_every=0",
            "utilities.state.series_interval=1.0", "global.num_threads=2",
            "utilities.state.record_file=/tmp/state.npz", "utilities.state.series_file=/tmp/series.npz",
            "utilities.state.state_series_file=/tmp/consolidated.npz", "global.output_dir=/tmp/column",
        ])
    pull = ["global.use_gpu=true", "executioner.n_steps=2", "executioner.t_end=4.0",
            "bcs.pull.values=[4.8e-7]", "global.num_threads=2", "global.output_dir=/tmp/cube"]
    tensile = "decks/lpbf/in718/k3/k3_cube_n6_seed1.i"
    report["stress_16"] = check("stress 16^3", 700, [tensile, cards, "mesh.divisions=[16,16,16]", *pull])
    for n in cube:
        report[f"stress_{n}_parts_only"] = check(f"stress {n}^3, parts only", 600, [
            "--no-single", tensile, cards, f"mesh.divisions=[{n},{n},{n}]", *pull])
    return report


@app.local_entrypoint()
def chain(steps: str = "20,60", cube: str = "24") -> None:
    out = chain_parts.remote([int(v) for v in steps.split(",")], [int(v) for v in cube.split(",") if v])
    print("CHAIN-RESULT " + json.dumps(out))


@app.local_entrypoint()
def cut(steps: str = "40,120", mesh: float = 5.0, ranks: str = "2,4", partition: str = "slabs",
        gpu: str = "A10G") -> None:
    function = {"A10G": decomposed, "L4": decomposed_l4, "A100": decomposed_a100,
                "T4": decomposed_t4}[gpu]
    print(json.dumps(function.remote([int(v) for v in steps.split(",")], mesh,
                                       [int(v) for v in ranks.split(",")], partition), indent=1))


@app.local_entrypoint()
def split(steps: str = "40,120", mesh: float = 5.0) -> None:
    print(json.dumps(processes.remote([int(v) for v in steps.split(",")], mesh), indent=1))


@app.function(gpu="A10G:4", cpu=8.0, memory=32768, timeout=900)
def diagnose() -> dict:
    """Where do the arrays of a pool live, and do four GPUs work at the same time at all?"""
    import sys
    import time
    from concurrent.futures import ThreadPoolExecutor

    import torch

    sys.path.insert(0, "/work/pinnacle/src")
    from pinnacle.backend import get_backend
    from pinnacle.backend.device_pool import DevicePool

    backend = get_backend("torch", precision="float64", device="cuda:0")
    devices = [f"cuda:{i}" for i in range(torch.cuda.device_count())]
    pool = DevicePool(backend, devices)
    report: dict = {
        "pool_devices": devices,
        "backend_tensor_devices": [str(b.tensor([0.0]).device) for b in pool.backends],
        "namespace_asarray_devices": [str(ns.asarray([0.0]).device) for ns in pool.namespaces],
        "namespace_index_devices": [str(ns.asarray([0], integer=True).device) for ns in pool.namespaces],
    }

    def work(device: str, n: int) -> float:
        x = torch.rand(n, 8, 3, dtype=torch.float64, device=device)
        total = torch.zeros((), dtype=torch.float64, device=device)
        for _ in range(300):
            total = total + (x * x).sum(dim=1).exp().sum()
        torch.cuda.synchronize(device)
        return float(total)

    n = 2_000_000
    for device in devices:
        work(device, 1000)
    tic = time.time(); work("cuda:0", n); one = time.time() - tic
    tic = time.time()
    for device in devices:
        work(device, n // len(devices))
    serial = time.time() - tic
    tic = time.time()
    with ThreadPoolExecutor(len(devices)) as ex:
        list(ex.map(lambda d: work(d, n // len(devices)), devices))
    threaded = time.time() - tic
    report["micro_benchmark_s"] = {"all_on_one_gpu": one, "quarters_one_after_another": serial,
                                   "quarters_in_threads": threaded}
    return report


@app.local_entrypoint()
def check() -> None:
    print(json.dumps(diagnose.remote(), indent=1))


@app.local_entrypoint()
def main(gpu: str = "A10G", counts: str = "1,2,4", steps: int = 200, mesh: float = 5.0,
         mechanics: bool = True, profile: bool = False) -> None:
    wanted = [int(v) for v in counts.split(",")]
    function = {"A10G": measure_a10g, "L4": measure_l4, "H100": measure_h100}[gpu]
    print(json.dumps(function.remote(wanted, steps, mesh, mechanics, profile), indent=1))
