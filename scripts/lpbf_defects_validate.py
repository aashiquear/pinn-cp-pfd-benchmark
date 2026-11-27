#!/usr/bin/env python3
"""Part V, defect indicators against the downloaded data sets (what each table allows).

usage: scripts/lpbf_defects_validate.py <command> [<raw data dir>] [<out dir>] [<figure dir>]

  gan      keyhole number of the published keyhole-depth rows (gan 2021 keyhole,
           Supplementary Data 1, measured absorptivity): how the rows fall into the
           regimes of the paper (conduction below 1.4, keyhole above 6), and the depth
           over spot radius in each regime
  modes    keyhole number and normalized enthalpy against the melt-pool mode labels of a
           literature compilation (akbari 2022 meltpoolnet, classification table): per
           material the area under the ROC curve (independent of the absorptivity), the
           best threshold, and the confusion counts at the paper's threshold for several
           assumed absorptivities
  overlap  the overlap index (H/W)^2 + (L/D)^2 from MEASURED pool width and depth of the
           same compilation (regression table) against its lack-of-fusion label and
           against relative density
  density  relative density of built parts (dovgyy 2021 lpbfdata compilation) against
           the keyhole number, alloy 718 and 316L (the table has no pool dimensions, so
           the overlap rule cannot be formed)
  pores    per-pore lists from X-ray CT of three build conditions of alloy 718 (NIST
           mds2-3083): counts and sizes beside the indicators of the three conditions
  all      every command; writes data/lpbf/defects/validation_summary.json

Defaults: raw data in data/lpbf/defects/raw, results in data/lpbf/defects, figures in
figures/lpbf. The raw tables are not tracked (sources: raw/*/SOURCE.md). The indicators
are PINNACLE's functions (`pinnacle.components.additive.defects`), which must be
importable; the workbooks need openpyxl, the figures matplotlib.

ASSUMPTIONS, stated because the tables do not give these values:
  - substrate temperature 298.15 K unless a build-plate temperature is given;
  - "beam D" / "Beam size" is the beam diameter 2 r0 (the definition is not recorded);
  - absorptivity: not measured in the compilations. The ROC area and the ranking inside
    one material do not depend on it (the keyhole number is proportional to it); every
    count at a fixed threshold is given for 0.35, 0.50 and 0.70;
  - the compilation writes the melting temperature of a material in kelvin in some rows
    and in degrees Celsius in others (two values 273 apart): the kelvin value is used,
    and a material with one value only is left out because its unit cannot be told;
  - alloy 718 in the density table: properties of the decks of this study at the
    liquidus; 316L: the constants of the compilation of (akbari 2022 meltpoolnet).
The labels of the compilation are melt-pool modes assigned by its compilers from pool
shapes (their notes: depth over width for keyhole, depth over layer for lack of fusion,
length over width for balling), not measured porosity.
"""

from __future__ import annotations

import csv
import json
import sys
from collections import Counter
from itertools import pairwise
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
ROOM = 298.15
ETAS = (0.35, 0.50, 0.70)
KE_KEYHOLE, KE_CONDUCTION = 6.0, 1.4
# alloy 718 at the liquidus (decks of this study, docs/16_lpbf_notebook.md, section 2)
IN718 = {"capacity": 6.3473e6, "conductivity": 26.80, "liquidus": 1610.15, "boiling": 3186.0}
# 316L: the constants of the compilation (density 8000, c_p 500, k 16.3, melting 1688 K)
SS316L = {"capacity": 8000.0 * 500.0, "conductivity": 16.3, "liquidus": 1688.0}


def num(text: object) -> float | None:
    try:
        value = float(str(text).strip())
    except (TypeError, ValueError):
        return None
    return value if np.isfinite(value) else None


def indicators(eta: float, power: float, speed: float, radius: float, liquidus: float,
               capacity: float, conductivity: float, substrate: float = ROOM) -> tuple[float, float]:
    """Keyhole number and normalized enthalpy from PINNACLE's functions (SI units)."""
    from pinnacle.components.additive.defects import keyhole_number, normalized_enthalpy

    diffusivity = conductivity / capacity
    return (
        keyhole_number(eta, power, speed, radius, liquidus, substrate, capacity, diffusivity),
        normalized_enthalpy(eta, power, speed, radius, liquidus, capacity, diffusivity),
    )


def roc(score: np.ndarray, positive: np.ndarray) -> dict:
    """ROC of `score` (larger = positive): area, best threshold (Youden), its rates."""
    order = np.argsort(-score, kind="stable")
    s, p = score[order], positive[order].astype(bool)
    n_pos, n_neg = int(p.sum()), int((~p).sum())
    if n_pos == 0 or n_neg == 0:
        return {"n_positive": n_pos, "n_negative": n_neg}
    # area by the rank statistic (ties count one half)
    ranks = np.empty(score.size)
    sorter = np.argsort(score, kind="stable")
    sorted_scores = score[sorter]
    start = 0
    while start < score.size:
        stop = start
        while stop + 1 < score.size and sorted_scores[stop + 1] == sorted_scores[start]:
            stop += 1
        ranks[sorter[start:stop + 1]] = 0.5 * (start + stop) + 1.0
        start = stop + 1
    area = (ranks[positive.astype(bool)].sum() - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)
    tpr = np.cumsum(p) / n_pos
    fpr = np.cumsum(~p) / n_neg
    last = np.r_[s[1:] != s[:-1], True]  # thresholds between distinct scores only
    youden = np.where(last, tpr - fpr, -np.inf)
    best = int(np.argmax(youden))
    return {
        "n_positive": n_pos, "n_negative": n_neg, "area": round(float(area), 3),
        "best_threshold": float(s[best]), "tpr_at_best": round(float(tpr[best]), 3),
        "fpr_at_best": round(float(fpr[best]), 3),
        "accuracy_at_best": round(float((tpr[best] * n_pos + (1 - fpr[best]) * n_neg) / (n_pos + n_neg)), 3),
        "curve": [fpr[last].round(4).tolist(), tpr[last].round(4).tolist()],
    }


def confusion(flag: np.ndarray, truth: np.ndarray) -> dict:
    flag, truth = flag.astype(bool), truth.astype(bool)
    tp, fp = int((flag & truth).sum()), int((flag & ~truth).sum())
    fn, tn = int((~flag & truth).sum()), int((~flag & ~truth).sum())
    return {"true_positive": tp, "false_positive": fp, "false_negative": fn, "true_negative": tn,
            "accuracy": round((tp + tn) / max(tp + fp + fn + tn, 1), 3)}


def plotter():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


# ----------------------------------------------------------------------------- gan

def cmd_gan(raw: Path, out: Path, figs: Path) -> dict:
    sys.path.insert(0, str(REPO / "scripts"))
    from lpbf_defects import gan_rows, with_ke

    rows = gan_rows(raw / "gan2021_keyhole_scaling" / "41467_2021_22704_MOESM4_ESM.xlsx")
    ke = with_ke(rows)
    depth = np.array([r["depth"] / r["radius"] for r in rows])
    regime = np.where(ke < KE_CONDUCTION, "conduction", np.where(ke <= KE_KEYHOLE, "transition", "keyhole"))
    report: dict = {"rows": len(rows), "note": "every row is a measured vapor depression; the workbook has no porosity or stability label",
                    "regimes": {}, "by_set": {}}
    for name in ("conduction", "transition", "keyhole"):
        pick = regime == name
        report["regimes"][name] = {
            "rows": int(pick.sum()),
            "depth_over_radius_min_median_max": [round(float(v), 3) for v in
                                                 (depth[pick].min(), np.median(depth[pick]), depth[pick].max())] if pick.any() else None,
        }
    for label in sorted({r["set"] for r in rows}):
        pick = np.array([r["set"] == label for r in rows])
        report["by_set"][label] = {"rows": int(pick.sum()), "Ke_min_max": [round(float(ke[pick].min()), 2), round(float(ke[pick].max()), 2)],
                                   "regimes": dict(Counter(regime[pick].tolist()))}
    # the two statements of the regime picture that the depths can test
    report["rows_deeper_than_one_radius"] = {
        "below_Ke_6": int(((depth > 1.0) & (ke <= KE_KEYHOLE)).sum()), "of": int((ke <= KE_KEYHOLE).sum()),
        "above_Ke_6": int(((depth > 1.0) & (ke > KE_KEYHOLE)).sum()), "of_above": int((ke > KE_KEYHOLE).sum()),
    }
    report["depth_over_diameter_at_least_one"] = roc(ke, depth >= 2.0)
    report["depth_over_diameter_at_least_one"].pop("curve", None)
    report["Ke_above_16_rows"] = int((ke > 16.0).sum())
    report["Ke_above_30_rows"] = int((ke > 30.0).sum())
    (out / "validation_gan_regimes.json").write_text(json.dumps(report, indent=1))
    return report


# --------------------------------------------------------------------------- modes

def meltpool_rows(raw: Path) -> tuple[list[dict], dict]:
    path = raw / "meltpoolnet_akbari2022" / "github_BaratiLab_MeltpoolNet" / "meltpoolnet_classification.csv"
    with path.open(encoding="utf-8", errors="replace") as fh:
        table = list(csv.reader(fh))
    head = [h.strip() for h in table[0]]
    col = {name: head.index(name) for name in ("Material", "Sub-process", "Power", "Velocity", "Hatch spacing",
                                                "layer thickness", "beam D", "density", "Cp", "k", "melting T",
                                                "minimal absorptivity", "meltpool shape")}
    melting: dict[str, set] = {}
    for r in table[1:]:
        value = num(r[col["melting T"]])
        if value:
            melting.setdefault(r[col["Material"]].strip(), set()).add(value)
    kelvin = {}
    for material, values in melting.items():
        pairs = [hi for hi in values for lo in values if abs(hi - lo - 273.0) < 1.5]
        if pairs:
            kelvin[material] = max(pairs)
    rows, skipped = [], Counter()
    for r in table[1:]:
        material = r[col["Material"]].strip()
        if r[col["Sub-process"]].strip() != "SLM":
            skipped["not laser (EBM)"] += 1
            continue
        if material not in kelvin:
            skipped["unit of the melting temperature not established"] += 1
            continue
        values = {k: num(r[col[k]]) for k in ("Power", "Velocity", "beam D", "density", "Cp", "k")}
        if any(v is None or v <= 0.0 for v in values.values()):
            skipped["power, speed, beam diameter or a property missing"] += 1
            continue
        rows.append({
            "material": material, "label": r[col["meltpool shape"]].strip(),
            "power": values["Power"], "speed": values["Velocity"] * 1e-3, "radius": 0.5e-6 * values["beam D"],
            "capacity": values["density"] * values["Cp"], "conductivity": values["k"], "liquidus": kelvin[material],
            "listed_absorptivity": num(r[col["minimal absorptivity"]]),
        })
    return rows, {"melting_temperature_K": kelvin, "skipped": dict(skipped), "rows_in_table": len(table) - 1}


def cmd_modes(raw: Path, out: Path, figs: Path) -> dict:
    rows, info = meltpool_rows(raw)
    unit = np.array([indicators(1.0, r["power"], r["speed"], r["radius"], r["liquidus"], r["capacity"], r["conductivity"]) for r in rows])
    ke1, enthalpy1 = unit[:, 0], unit[:, 1]
    label = np.array([r["label"] for r in rows])
    material = np.array([r["material"] for r in rows])
    report: dict = {"table": "meltpoolnet_classification.csv", **info, "rows_used": len(rows),
                    "labels_used": dict(Counter(label.tolist())), "materials": {}, "pooled": {}}
    curves = {}
    for name in sorted(set(material.tolist()), key=lambda m: -int((material == m).sum())):
        pick = material == name
        entry: dict = {"rows": int(pick.sum()), "labels": dict(Counter(label[pick].tolist())),
                       "Ke_per_unit_absorptivity_min_median_max": [round(float(v), 2) for v in
                                                                   (ke1[pick].min(), np.median(ke1[pick]), ke1[pick].max())]}
        for title, positive, negative in (
            ("keyhole_against_desirable", "keyhole", ("desirable",)),
            ("keyhole_against_all_other_modes", "keyhole", ("desirable", "LOF", "balling")),
            ("desirable_against_lack_of_fusion", "desirable", ("LOF",)),
        ):
            sub = pick & np.isin(label, (positive, *negative))
            result = roc(ke1[sub], label[sub] == positive)
            if "area" not in result:
                entry[title] = result
                continue
            if title == "keyhole_against_desirable":
                curves[name] = result["curve"]
            result.pop("curve")
            result["best_threshold_Ke_at_absorptivity"] = {f"{eta:.2f}": round(result["best_threshold"] * eta, 2) for eta in ETAS}
            result["best_threshold_Ke_per_unit_absorptivity"] = round(result.pop("best_threshold"), 2)
            if positive == "keyhole":
                result["at_Ke_6"] = {f"{eta:.2f}": confusion(ke1[sub] * eta > KE_KEYHOLE, label[sub] == positive) for eta in ETAS}
            entry[title] = result
        listed = {r["listed_absorptivity"] for r, p in zip(rows, pick) if p and r["listed_absorptivity"]}
        entry["absorptivity_values_listed_in_table"] = sorted(round(v, 3) for v in listed)[:4] + (["..."] if len(listed) > 4 else [])
        report["materials"][name] = entry
    # Alloy 718 again with the properties of this study at the liquidus: the compilation
    # lists room-temperature values (rho c_p 3.56e6 J/(m3 K), k 11.2 W/(m K), melting
    # 1673 K), and the keyhole number of a row is larger with them by a constant factor.
    pick = material == "IN718"
    if pick.any():
        first = next(r for r in rows if r["material"] == "IN718")
        factor = (indicators(1.0, first["power"], first["speed"], first["radius"], IN718["liquidus"], IN718["capacity"], IN718["conductivity"])[0]
                  / indicators(1.0, first["power"], first["speed"], first["radius"], first["liquidus"], first["capacity"], first["conductivity"])[0])
        sub = pick & np.isin(label, ("keyhole", "desirable"))
        truth = label[sub] == "keyhole"
        entry = report["materials"]["IN718"]
        entry["with_liquidus_properties_of_this_study"] = {
            "Ke_factor_against_the_table_properties": round(float(factor), 3),
            "best_threshold_Ke_at_absorptivity": {f"{eta:.2f}": round(entry["keyhole_against_desirable"]["best_threshold_Ke_per_unit_absorptivity"] * factor * eta, 2) for eta in ETAS},
            "keyhole_against_desirable_at_Ke_6": {f"{eta:.2f}": confusion(ke1[sub] * factor * eta > KE_KEYHOLE, truth) for eta in ETAS},
        }
    for eta in ETAS:
        sub = np.isin(label, ("keyhole", "desirable"))
        truth = label[sub] == "keyhole"
        by_ke = roc(ke1[sub] * eta, truth)
        by_h = roc(enthalpy1[sub] * eta, truth)
        report["pooled"][f"{eta:.2f}"] = {
            "keyhole_against_desirable_rows": int(sub.sum()),
            "keyhole_number": {k: by_ke[k] for k in ("area", "best_threshold", "accuracy_at_best")} | {"at_Ke_6": confusion(ke1[sub] * eta > KE_KEYHOLE, truth)},
            "normalized_enthalpy": {k: by_h[k] for k in ("area", "best_threshold", "accuracy_at_best")},
        }
    (out / "validation_meltpool_modes.json").write_text(json.dumps(report, indent=1))
    with (out / "validation_meltpool_modes.csv").open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["material", "label", "power_W", "speed_m_s", "beam_radius_m", "Ke_per_unit_absorptivity", "normalized_enthalpy_per_unit_absorptivity"])
        for r, k, h in zip(rows, ke1, enthalpy1):
            writer.writerow([r["material"], r["label"], r["power"], r["speed"], r["radius"], f"{k:.5g}", f"{h:.5g}"])

    plt = plotter()
    names = list(report["materials"])
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(11.0, 4.4), gridspec_kw={"width_ratios": [1.5, 1.0]})
    colors = {"LOF": "tab:blue", "balling": "tab:purple", "desirable": "tab:green", "keyhole": "tab:red"}
    rng = np.random.default_rng(0)
    for i, name in enumerate(names):
        for j, mode in enumerate(("LOF", "balling", "desirable", "keyhole")):
            pick = (material == name) & (label == mode)
            y = i + (j - 1.5) * 0.19 + rng.uniform(-0.06, 0.06, int(pick.sum()))
            ax.scatter(ke1[pick] * 0.35, y, s=7, color=colors[mode], alpha=0.6, linewidths=0, label=mode if i == 0 else None)
    ax.axvline(KE_KEYHOLE, color="k", lw=0.8, ls="--")
    ax.axvline(KE_CONDUCTION, color="0.5", lw=0.8, ls=":")
    ax.set_xscale("log")
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names)
    ax.set_xlabel("keyhole number Ke at an absorptivity of 0.35 (dashed: 6, dotted: 1.4)")
    ax.legend(title="mode label of the compilation", fontsize=8, title_fontsize=8, loc="lower right")
    for name, (fpr, tpr) in curves.items():
        area = report["materials"][name]["keyhole_against_desirable"]["area"]
        bx.plot([0.0, *fpr], [0.0, *tpr], lw=1.2, label=f"{name} (area {area:.2f})")
    bx.plot([0, 1], [0, 1], color="0.7", lw=0.7)
    bx.set_xlabel("share of 'desirable' rows flagged")
    bx.set_ylabel("share of 'keyhole' rows flagged")
    bx.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(figs / "defects_validation_modes.png", dpi=180)
    return report


# ------------------------------------------------------------------------- overlap

def cmd_overlap(raw: Path, out: Path, figs: Path) -> dict:
    from pinnacle.components.additive.defects import overlap_index

    path = raw / "meltpoolnet_akbari2022" / "github_BaratiLab_MeltpoolNet" / "meltpoolnet_regression.csv"
    with path.open(encoding="utf-8", errors="replace") as fh:
        table = list(csv.reader(fh))
    head = [h.strip() for h in table[0]]
    col = {k: head.index(k) for k in ("Material", "Process", "layer thickness", "Hatch spacing", "depth of meltpool",
                                      "width of melt pool", "meltpool shape", "relative density")}
    rows = []
    for r in table[1:]:
        if r[col["Process"]].strip() != "PBF":
            continue
        v = {k: num(r[col[k]]) for k in ("layer thickness", "Hatch spacing", "depth of meltpool", "width of melt pool")}
        if any(x is None or x <= 0.0 for x in v.values()):
            continue
        rows.append({
            "material": r[col["Material"]].strip(), "label": r[col["meltpool shape"]].strip(),
            "density": num(r[col["relative density"]]),
            "layer": v["layer thickness"], "hatch": v["Hatch spacing"], "depth": v["depth of meltpool"], "width": v["width of melt pool"],
            "index": overlap_index(v["Hatch spacing"], v["width of melt pool"], v["layer thickness"], v["depth of meltpool"]),
        })
    index = np.array([r["index"] for r in rows])
    label = np.array([r["label"] for r in rows])
    density = np.array([np.nan if r["density"] is None else r["density"] for r in rows])
    depth_over_layer = np.array([r["depth"] / r["layer"] for r in rows])
    hatch_over_width = np.array([r["hatch"] / r["width"] for r in rows])
    report: dict = {"table": "meltpoolnet_regression.csv", "rows_with_hatch_layer_width_depth": len(rows),
                    "materials": dict(Counter(r["material"] for r in rows).most_common()),
                    "labels": dict(Counter(label.tolist()))}
    labelled = label != ""
    truth = label[labelled] == "LOF"
    report["label_test"] = {
        "rows": int(labelled.sum()),
        "index_above_1_against_LOF_label": confusion(index[labelled] > 1.0, truth),
        "roc_of_the_index": {k: v for k, v in roc(index[labelled], truth).items() if k != "curve"},
        "LOF_rows_index_min_median_max": [round(float(v), 2) for v in (index[label == "LOF"].min(), np.median(index[label == "LOF"]), index[label == "LOF"].max())],
        "other_rows_index_min_median_max": [round(float(v), 2) for v in (index[labelled & (label != "LOF")].min(), np.median(index[labelled & (label != "LOF")]), index[labelled & (label != "LOF")].max())],
        "circularity": {
            "note": "the compilers' note gives depth over layer below 1.5 as their lack-of-fusion rule; the (L/D) term of the index is the same ratio",
            "LOF_rows_with_depth_over_layer_below_1.5": int(((depth_over_layer < 1.5) & (label == "LOF")).sum()),
            "LOF_rows": int((label == "LOF").sum()),
            "other_labelled_rows_with_depth_over_layer_below_1.5": int(((depth_over_layer < 1.5) & labelled & (label != "LOF")).sum()),
            "LOF_rows_flagged_by_the_hatch_term_alone_(H/W)^2_above_1": int(((hatch_over_width > 1.0) & (label == "LOF")).sum()),
        },
    }
    has = np.isfinite(density)
    dense = density[has] >= 99.0
    inside = index[has] <= 1.0

    def rank(a: np.ndarray) -> np.ndarray:
        return np.argsort(np.argsort(a)).astype(float)

    report["density_test"] = {
        "rows": int(has.sum()), "materials": dict(Counter(r["material"] for r, h in zip(rows, has) if h)),
        "index_at_most_1": {"rows": int(inside.sum()), "density_at_least_99": int((dense & inside).sum()),
                            "median_density": round(float(np.median(density[has][inside])), 2) if inside.any() else None},
        "index_above_1": {"rows": int((~inside).sum()), "density_at_least_99": int((dense & ~inside).sum()),
                          "median_density": round(float(np.median(density[has][~inside])), 2) if (~inside).any() else None},
        "rank_correlation_of_index_and_density": round(float(np.corrcoef(rank(index[has]), rank(density[has]))[0, 1]), 3),
        "index_min_max": [round(float(index[has].min()), 3), round(float(index[has].max()), 3)],
        "density_min_max": [float(density[has].min()), float(density[has].max())],
    }
    with (out / "validation_overlap_index.csv").open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["material", "label", "layer_um", "hatch_um", "depth_um", "width_um", "overlap_index", "relative_density_pct"])
        for r in rows:
            writer.writerow([r["material"], r["label"], r["layer"], r["hatch"], r["depth"], r["width"], f"{r['index']:.4g}", "" if r["density"] is None else r["density"]])
    (out / "validation_overlap_index.json").write_text(json.dumps(report, indent=1))

    plt = plotter()
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(10.5, 4.2))
    colors = {"LOF": "tab:blue", "balling": "tab:purple", "desirable": "tab:green", "keyhole": "tab:red"}
    for mode, color in colors.items():
        pick = label == mode
        ax.scatter(hatch_over_width[pick] ** 2, (1.0 / depth_over_layer[pick]) ** 2, s=14, color=color, alpha=0.7, linewidths=0, label=f"{mode} ({int(pick.sum())})")
    line = np.linspace(0.0, 1.0, 100)
    ax.plot(line, 1.0 - line, color="k", lw=0.9, label="index = 1")
    ax.set_xscale("symlog", linthresh=0.05)
    ax.set_yscale("symlog", linthresh=0.05)
    ax.set_xlabel("(hatch / measured width)$^2$")
    ax.set_ylabel("(layer / measured depth)$^2$")
    ax.legend(title="mode label of the compilation", fontsize=8, title_fontsize=8)
    for name in sorted({r["material"] for r, h in zip(rows, has) if h}):
        pick = has & np.array([r["material"] == name for r in rows])
        bx.scatter(index[pick], density[pick], s=16, alpha=0.8, linewidths=0, label=f"{name} ({int(pick.sum())})")
    bx.axvline(1.0, color="k", lw=0.9)
    bx.axhline(99.0, color="0.6", lw=0.7, ls="--")
    bx.set_xlabel("overlap index from measured width and depth")
    bx.set_ylabel("relative density (%)")
    bx.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(figs / "defects_validation_overlap.png", dpi=180)
    return report


# ------------------------------------------------------------------------- density

def cmd_density(raw: Path, out: Path, figs: Path) -> dict:
    import openpyxl

    book = openpyxl.load_workbook(raw / "zenodo_6514897_lpbf_parameters_properties" / "Process Parameters and Properties LPBF v1.1.xlsx",
                                  read_only=True, data_only=True)
    table = list(book.worksheets[0].iter_rows(values_only=True))
    head = [str(h).strip() if h is not None else "" for h in table[1]]
    col = {k: next(i for i, h in enumerate(head) if h.startswith(k)) for k in
           ("Material", "Laser power", "Laser speed", "Layer thickness", "Hatch spacing", "Beam size", "Build plate temperature", "Density/Consolidation", "Link to paper")}
    report: dict = {"table": "Process Parameters and Properties LPBF v1.1.xlsx", "rows_in_table": len(table) - 2,
                    "pool_width_or_depth_columns": "none: the overlap rule cannot be formed on this table", "materials": {}}
    plt = plotter()
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), sharey=True)
    lines = []
    for ax, (name, key, prop) in zip(axes, (("alloy 718", "IN718", IN718), ("316L", "316L", SS316L))):
        ke, dens, papers = [], [], set()
        counts = Counter()
        for r in table[2:]:
            if str(r[col["Material"]]).strip() != key:
                continue
            counts["rows"] += 1
            d = num(r[col["Density/Consolidation"]])
            p, v, beam = num(r[col["Laser power"]]), num(r[col["Laser speed"]]), num(r[col["Beam size"]])
            if d is None or d <= 0.0 or d > 100.0:
                counts["no usable density"] += 1
                continue
            if None in (p, v, beam) or min(p, v, beam) <= 0.0:
                counts["power, speed or beam size missing"] += 1
                continue
            plate = num(r[col["Build plate temperature"]])
            substrate = ROOM if plate is None else plate + 273.15
            if substrate >= prop["liquidus"]:
                continue
            ke.append(indicators(0.35, p, v * 1e-3, 0.5e-6 * beam, prop["liquidus"], prop["capacity"], prop["conductivity"], substrate)[0])
            dens.append(d)
            papers.add(str(r[col["Link to paper"]]))
            lines.append([key, p, v, beam, "" if plate is None else plate, f"{ke[-1]:.4g}", d])
        ke_a, dens_a = np.array(ke), np.array(dens)
        edges = [0.0, 3.0, 6.0, 10.0, 16.0, 30.0, np.inf]
        bins = []
        for lo, hi in pairwise(edges):
            pick = (ke_a >= lo) & (ke_a < hi)
            if pick.any():
                bins.append({"Ke_from": lo, "Ke_to": None if not np.isfinite(hi) else hi, "rows": int(pick.sum()),
                             "median_density": round(float(np.median(dens_a[pick])), 2),
                             "share_below_99": round(float(np.mean(dens_a[pick] < 99.0)), 3)})
        report["materials"][name] = {"counts": dict(counts), "rows_used": int(ke_a.size), "papers": len(papers),
                                     "absorptivity_assumed": 0.35, "bins_of_Ke": bins,
                                     "density_at_least_99": {"rows": int((dens_a >= 99.0).sum()),
                                                             "Ke_p05_median_p95": [round(float(v), 1) for v in np.percentile(ke_a[dens_a >= 99.0], [5, 50, 95])] if (dens_a >= 99.0).any() else None},
                                     "density_below_95": {"rows": int((dens_a < 95.0).sum()),
                                                          "Ke_p05_median_p95": [round(float(v), 1) for v in np.percentile(ke_a[dens_a < 95.0], [5, 50, 95])] if (dens_a < 95.0).any() else None}}
        ax.scatter(ke_a, dens_a, s=9, alpha=0.5, linewidths=0, color="tab:blue")
        centers = [np.sqrt(max(b["Ke_from"], 1.0) * (b["Ke_to"] or 60.0)) for b in bins]
        ax.plot(centers, [b["median_density"] for b in bins], color="k", marker="s", ms=4, lw=1.0, label="median per bin of Ke")
        ax.axvline(KE_KEYHOLE, color="k", lw=0.7, ls="--")
        ax.axhline(99.0, color="0.6", lw=0.7, ls="--")
        ax.set_xscale("log")
        ax.set_ylim(80.0, 100.5)
        ax.set_xlabel(f"keyhole number Ke, {name}, absorptivity 0.35 assumed")
        ax.legend(fontsize=8, loc="lower right")
    axes[0].set_ylabel("relative density (%), values below 80 not shown")
    fig.tight_layout()
    fig.savefig(figs / "defects_validation_density.png", dpi=180)
    with (out / "validation_density_keyhole_number.csv").open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["material", "power_W", "speed_mm_s", "beam_size_um", "plate_C", "Ke_at_absorptivity_0.35", "relative_density_pct"])
        writer.writerows(lines)
    (out / "validation_density_keyhole_number.json").write_text(json.dumps(report, indent=1))
    return report


# --------------------------------------------------------------------------- pores

def cmd_pores(raw: Path, out: Path, figs: Path) -> dict:
    from pinnacle.components.additive.defects import keyhole_depth

    base = raw / "nist_mds2-3083_in718_asbuilt_psp"
    with (base / "Paper1" / "Surface Roughness" / "sr.csv").open(encoding="latin-1") as fh:
        table = list(csv.reader(fh))
    head = table[0]
    col = {k: head.index(k) for k in ("Label", "Power (Fill)", "Scan Speed (Fill)", "Layer Thickness", "Build Orientation", "Spot Size")}
    stated = {"P1": 62.0, "P2": 62.0, "P4": 38.0}  # J/mm3, Paper2/paper2_readme.txt
    report: dict = {
        "record": "NIST mds2-3083, Paper2/Micro-CT data/Data",
        "columns": "no header in the pore files and no description in the record: column 3 equals the diameter of a sphere of the volume in column 2 "
                   "(units of the reconstruction, pixel 0.998 um per the scan read-me, an inference); column 1 is not interpreted",
        "parameters": "fill power, speed, layer and spot of the labels P1, P2, P4 in Paper1/Surface Roughness/sr.csv; with a hatch of 0.10 mm they give the "
                      "energy densities the Paper2 read-me states for P1, P2, P4 (62, 62, 38 J/mm3), which is the evidence that the labels are the same "
                      "conditions and that the hatch was 0.10 mm (the hatch itself is not in the downloaded files)",
        "conditions": {},
    }
    plt = plotter()
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(10.5, 4.2))
    for name, color in (("P1", "tab:green"), ("P2", "tab:blue"), ("P4", "tab:red")):
        row = next(r for r in table[2:] if r[col["Label"]] == name)
        power, speed = float(row[col["Power (Fill)"]]), float(row[col["Scan Speed (Fill)"]])
        layer, spot = float(row[col["Layer Thickness"]]), float(row[col["Spot Size"]])
        entry: dict = {"power_W": power, "speed_mm_s": speed, "layer_um": layer, "spot_um": spot,
                       "orientation_in_file_deg": row[col["Build Orientation"]],
                       "energy_density_stated_J_mm3": stated[name],
                       "energy_density_with_hatch_100um_J_mm3": round(power / (speed * 0.10 * layer * 1e-3), 1)}
        for eta in ETAS:
            ke, h = indicators(eta, power, speed * 1e-3, 0.5e-6 * spot, IN718["liquidus"], IN718["capacity"], IN718["conductivity"])
            entry[f"absorptivity_{eta:.2f}"] = {"Ke": round(ke, 1), "normalized_enthalpy": round(h, 1),
                                                "depression_depth_by_the_law_um": round(1e6 * keyhole_depth(ke, 0.5e-6 * spot), 0),
                                                "depth_by_the_law_over_layer": round(1e6 * keyhole_depth(ke, 0.5e-6 * spot) / layer, 1)}
        for kind, style in (("LOF", "-"), ("gas", ":")):
            with (base / "Paper2" / "Micro-CT data" / "Data" / f"{name}_{kind}.csv").open() as fh:
                data = np.array([[float(x) for x in line] for line in csv.reader(fh) if line])
            diameter, volume = data[:, 2], data[:, 1]
            entry[f"{kind}_pores"] = {"count": int(data.shape[0]), "column2_sum": float(volume.sum()),
                                      "column3_median": round(float(np.median(diameter)), 1), "column3_p95": round(float(np.percentile(diameter, 95)), 1),
                                      "column3_max": round(float(diameter.max()), 1),
                                      "column3_is_sphere_diameter_of_column2_within": round(float(np.max(np.abs(diameter / np.cbrt(6.0 * volume / np.pi) - 1.0))), 4)}
            target = ax if kind == "LOF" else bx
            target.plot(np.sort(diameter), np.arange(diameter.size, 0, -1), style if kind == "gas" else "-", color=color,
                        label=f"{name}: {power:.0f} W, {speed:.0f} mm/s, {stated[name]:.0f} J/mm$^3$ ({diameter.size} pores)")
        report["conditions"][name] = entry
    for target, kind in ((ax, "pores labelled lack of fusion"), (bx, "pores labelled gas")):
        target.set_xscale("log")
        target.set_yscale("log")
        target.set_xlabel(f"{kind}: column 3 of the pore list (equivalent diameter, units of the scan)")
        target.set_ylabel("number of pores at least this large")
        target.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(figs / "defects_validation_pores.png", dpi=180)
    lof = {n: report["conditions"][n]["LOF_pores"]["column2_sum"] for n in report["conditions"]}
    report["lack_of_fusion_volume_ratio_P4_over_P2"] = round(lof["P4"] / lof["P2"], 1)
    report["lack_of_fusion_volume_ratio_P4_over_P1"] = round(lof["P4"] / lof["P1"], 1)
    report["what_the_indicators_can_say"] = (
        "the keyhole number is lowest for P4, the condition with the lack-of-fusion porosity (6.91 % by volume per the read-me), but it is above the keyhole "
        "threshold for all three conditions at every assumed absorptivity: it does not flag lack of fusion. The overlap index needs the pool width and depth, "
        "which the record does not give for these conditions; it was not computed from a model pool."
    )
    (out / "validation_nist_pores.json").write_text(json.dumps(report, indent=1))
    return report


def cracking_note(raw: Path) -> dict:
    hits = []
    for path in raw.rglob("*"):
        if path.suffix.lower() in (".csv", ".txt") and path.is_file() and path.stat().st_size < 5_000_000:
            try:
                if "crack" in path.read_text(encoding="latin-1").lower():
                    hits.append(str(path.relative_to(raw)))
            except OSError:
                pass
    return {"tables_with_the_word_crack": hits,
            "statement": "no downloaded table carries a cracking label or measure (the word occurs in descriptive read-me text only); the cracking "
                         "indicators stay verified against the published Al-Cu curve of (rappaz 1999 hottearing) and are not validated"}


def main() -> None:
    args = sys.argv[1:]
    if not args or args[0] not in ("gan", "modes", "overlap", "density", "pores", "all"):
        raise SystemExit(__doc__)
    raw = Path(args[1]) if len(args) > 1 else REPO / "data" / "lpbf" / "defects" / "raw"
    out = Path(args[2]) if len(args) > 2 else REPO / "data" / "lpbf" / "defects"
    figs = Path(args[3]) if len(args) > 3 else REPO / "figures" / "lpbf"
    out.mkdir(parents=True, exist_ok=True)
    figs.mkdir(parents=True, exist_ok=True)
    commands = {"gan": cmd_gan, "modes": cmd_modes, "overlap": cmd_overlap, "density": cmd_density, "pores": cmd_pores}
    if args[0] == "all":
        summary = {name: function(raw, out, figs) for name, function in commands.items()}
        summary["cracking"] = cracking_note(raw)
        (out / "validation_summary.json").write_text(json.dumps(summary, indent=1))
        print(json.dumps(summary, indent=1))
    else:
        print(json.dumps(commands[args[0]](raw, out, figs), indent=1))


if __name__ == "__main__":
    main()
