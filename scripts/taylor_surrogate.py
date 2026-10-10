#!/usr/bin/env python3
"""Taylor material-point surrogate: OAT sensitivities and a bounded fit to the experiment.

Reads a surrogate deck (see decks/cp_study/taylor_screen.i): fixed CP values,
free parameters as 'lower upper start', the experiment markers, and the
aggregate settings. Stages:

  oat   one-at-a-time: each free parameter at its lower bound, start and
        upper bound (others at start); writes oat.csv with the engineering
        stress at every marker strain and the RMSE against the markers.
  fit   damped Gauss-Newton in log-parameter space with box bounds, finite-
        difference Jacobian evaluated in parallel; writes fit.json (path,
        final set, per-marker errors) and fit_curve.csv.

The aggregate is the archived first-campaign routine (isostrain Taylor
average over N random grains with the lateral stretch solved so the
transverse stress vanishes; engineering stress = sigma_zz * lat^2), ported
here so the current study carries no live dependency on archive/.

Usage: taylor_surrogate.py DECK.i oat|fit|all
"""
from __future__ import annotations

import csv
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

STUDY = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(STUDY / "scripts"))
from polycrystal_plate import _resolve  # noqa: E402


def read_deck(path: Path) -> dict:
    from pinnacle.io.i_parser import InputParser
    from pinnacle.io.translator import Translator
    cfg = Translator().translate(InputParser().parse(path.read_text()))
    for k in ("surrogate", "fixed", "free"):
        if k not in cfg:
            raise SystemExit(f"{path}: missing [{k.capitalize()}]")
    return cfg


def random_orientations(n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    phi1 = rng.uniform(0.0, 2.0 * np.pi, n)
    Phi = np.arccos(rng.uniform(-1.0, 1.0, n))
    phi2 = rng.uniform(0.0, 2.0 * np.pi, n)
    return np.degrees(np.column_stack([phi1, Phi, phi2]))


def load_markers(path: Path):
    xs, ys = [], []
    with path.open() as fh:
        for row in csv.DictReader(r for r in fh if not r.startswith("#")):
            xs.append(float(row["strain_pct"]) / 100.0)
            ys.append(float(row["stress_MPa"]))
    return np.asarray(xs), np.asarray(ys)


def taylor_curve(params: dict, eulers: np.ndarray, sample_at: np.ndarray, edot: float,
                 n_increments: int, lat_tol: float = 1.0e-3, max_iter: int = 6):
    """Engineering stress (MPa) at ``sample_at``; NaN if the path diverged."""
    import warnings
    warnings.filterwarnings("ignore")
    from pinnacle import registry
    from pinnacle.backend import get_backend
    from pinnacle.plugins import discover_components
    discover_components(force=True)
    backend = get_backend("numpy", precision="float64", device="cpu", seed=0)
    mats = [registry.create("material", "dislocation_density_cp", backend=backend,
                            config={"material": dict(params),
                                    "crystal": {"structure": "fcc", "euler_angles": [float(x) for x in e]}})
            for e in eulers]
    states = [m.initial_state() for m in mats]
    n = len(mats)

    def aggregate(F, states_in, dt):
        Ft = backend.tensor(F, dtype="float64")
        acc = np.zeros((3, 3)); out = []
        for m, st in zip(mats, states_in):
            cauchy, new_st, _ = m.compute(Ft, st, dt=dt)
            acc += np.asarray(backend.to_numpy(cauchy), dtype=float); out.append(new_st)
        return acc / n, out

    eps_max = float(np.max(sample_at))
    path = np.linspace(eps_max / n_increments, eps_max, n_increments)
    dt = (eps_max / n_increments) / edot
    he, heng = [0.0], [0.0]
    lat = 1.0
    for eps in path:
        lam_z = 1.0 + eps
        trial = states; sig = None
        for _ in range(max_iter):
            sig, trial = aggregate(np.diag([lat, lat, lam_z]), states, dt)
            r = 0.5 * (sig[0, 0] + sig[1, 1])
            if abs(r) < lat_tol * max(1.0, abs(sig[2, 2])):
                break
            dl = 1.0e-6 * max(lat, 1.0)
            sig2, _ = aggregate(np.diag([lat + dl, lat + dl, lam_z]), states, dt)
            drdl = (0.5 * (sig2[0, 0] + sig2[1, 1]) - r) / dl
            if not np.isfinite(drdl) or abs(drdl) < 1.0e-12:
                break
            step = float(np.clip(r / drdl, -0.05 * lat, 0.05 * lat))
            lat = max(lat - step, 0.2)
        if sig is None or not np.isfinite(sig).all():
            break
        states = trial
        he.append(float(eps)); heng.append(float(sig[2, 2]) * lat * lat)
    he_a = np.asarray(he)
    if he_a.size < 3 or he_a[-1] < 0.98 * eps_max:
        return np.full(sample_at.shape, np.nan)
    return np.interp(sample_at, he_a, np.asarray(heng))


def _eval(args):
    params, eulers, sample, edot, ninc = args
    return taylor_curve(params, eulers, sample, edot, ninc)


class Surrogate:
    def __init__(self, cfg: dict, deck: Path):
        s = cfg["surrogate"]
        self.cfg = cfg
        self.out = _resolve(str(s["out_dir"]), deck); self.out.mkdir(parents=True, exist_ok=True)
        self.ex_e, self.ex_s = load_markers(_resolve(str(s["markers"]), deck))
        self.eulers = random_orientations(int(s.get("n_grains", 30)), int(s.get("seed", 0)))
        self.edot = float(s.get("strain_rate", 1e-3)); self.ninc = int(s.get("n_increments", 200))
        self.workers = int(s.get("workers", 8))
        self.fixed = {k: float(v) for k, v in cfg["fixed"].items()}
        self.free = {}
        for k, v in cfg["free"].items():
            lo, hi, x0 = (float(t) for t in (v if isinstance(v, list) else str(v).split()))
            self.free[k] = (lo, hi, x0)
        self.names = list(self.free)

    def params(self, x: dict) -> dict:
        p = dict(self.fixed); p.update(x); return p

    def evaluate_many(self, sets: list[dict]) -> list[np.ndarray]:
        jobs = [(self.params(x), self.eulers, self.ex_e, self.edot, self.ninc) for x in sets]
        with ProcessPoolExecutor(max_workers=self.workers) as ex:
            return list(ex.map(_eval, jobs))

    def rmse(self, curve: np.ndarray) -> float:
        return float(np.sqrt(np.nanmean((curve - self.ex_s) ** 2)))

    def oat(self):
        start = {k: v[2] for k, v in self.free.items()}
        sets, labels = [start], ["start"]
        for k, (lo, hi, x0) in self.free.items():
            for tag, val in (("lo", lo), ("hi", hi)):
                x = dict(start); x[k] = val; sets.append(x); labels.append(f"{k}={val:g}")
        t0 = time.time(); curves = self.evaluate_many(sets); wall = time.time() - t0
        rows = []
        for lab, x, c in zip(labels, sets, curves):
            rows.append({"case": lab, **x, "rmse_MPa": self.rmse(c),
                         **{f"S_{e*100:.2g}pct": float(v) for e, v in zip(self.ex_e, c)}})
        with (self.out / "oat.csv").open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
        print(f"OAT: {len(sets)} curves in {wall:.0f} s -> {self.out/'oat.csv'}")
        for r in rows:
            print(f"  {r['case']:>28}  rmse {r['rmse_MPa']:7.1f}  " + " ".join(f"{r[k]:6.0f}" for k in rows[0] if k.startswith("S_")))
        print("  experiment:                     " + " ".join(f"{v:6.0f}" for v in self.ex_s))

    def fit(self, max_iter: int = 6, rel_step: float = 0.15):
        lo = np.log([self.free[k][0] + 1e-12 for k in self.names]); hi = np.log([self.free[k][1] for k in self.names])
        z = np.log([max(self.free[k][2], 1e-12) for k in self.names])
        path = []
        for it in range(max_iter):
            sets = [dict(zip(self.names, np.exp(z)))]
            for j in range(len(z)):
                zj = z.copy(); zj[j] = min(zj[j] + rel_step, hi[j]) if zj[j] + rel_step <= hi[j] else zj[j] - rel_step
                sets.append(dict(zip(self.names, np.exp(zj))))
            curves = self.evaluate_many(sets)
            f0 = curves[0]; r0 = f0 - self.ex_s
            if not np.all(np.isfinite(f0)):
                raise SystemExit("surrogate diverged at the current point")
            J = np.zeros((len(self.ex_s), len(z)))
            for j in range(len(z)):
                dz = np.log(sets[j + 1][self.names[j]]) - z[j]
                J[:, j] = (curves[j + 1] - f0) / dz
            rmse = float(np.sqrt(np.mean(r0 ** 2)))
            path.append({"iter": it, "params": dict(zip(self.names, np.exp(z).tolist())), "rmse_MPa": rmse,
                         "curve": f0.tolist()})
            print(f"iter {it}: rmse {rmse:.1f}  " + "  ".join(f"{k}={v:.4g}" for k, v in zip(self.names, np.exp(z))))
            # Levenberg-Marquardt step with bounds by clipping, accepted if it improves.
            lam = 1e-3 * np.trace(J.T @ J) / len(z)
            step = -np.linalg.solve(J.T @ J + lam * np.eye(len(z)), J.T @ r0)
            step = np.clip(step, -0.7, 0.7)
            z_new = np.clip(z + step, lo, hi)
            if np.max(np.abs(z_new - z)) < 1e-3:
                break
            z = z_new
        final = dict(zip(self.names, np.exp(z).tolist()))
        c = self.evaluate_many([final])[0]
        rmse = self.rmse(c)
        rec = {"objective": "RMSE of engineering stress at the experiment markers (MPa)",
               "final": final, "rmse_MPa": rmse,
               "per_marker": [{"strain": float(e), "experiment": float(s), "surrogate": float(v),
                               "rel_err": float((v - s) / s)} for e, s, v in zip(self.ex_e, self.ex_s, c)],
               "path": path, "fixed": self.fixed, "bounds": {k: v[:2] for k, v in self.free.items()},
               "n_grains": len(self.eulers), "strain_rate": self.edot, "n_increments": self.ninc}
        (self.out / "fit.json").write_text(json.dumps(rec, indent=1))
        print(f"final rmse {rmse:.1f} MPa; {final}; wrote {self.out/'fit.json'}")


def main(argv):
    if len(argv) != 2:
        print(__doc__); return 2
    deck = Path(argv[0]).resolve(); sur = Surrogate(read_deck(deck), deck)
    for st in (["oat", "fit"] if argv[1] == "all" else [argv[1]]):
        getattr(sur, st)()
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
