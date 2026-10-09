#!/usr/bin/env python3
"""Summarize a MOOSE Exodus output: timesteps, variables, extent, and the
nominal stress-strain history (mean element stress vs. top-face displacement).

Usage: exodus_summary.py FILE.e [--pull-axis y] [--length L]
"""
import argparse
import netCDF4 as nc
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("file")
ap.add_argument("--pull-axis", default="y")
ap.add_argument("--length", type=float, default=None, help="gauge length along the pull axis (default: extent)")
ap.add_argument("--csv", default=None, help="write step,t,nominal_strain,mean_stress for every written step to this CSV")
a = ap.parse_args()
d = nc.Dataset(a.file)
t = d.variables["time_whole"][:]
n_written = 1 + int(np.count_nonzero(np.asarray(t[1:]) > 0)) if len(t) > 1 else len(t)
t = t[:n_written]
def names(k):
    if k not in d.variables:
        return []
    return ["".join(x.decode() if isinstance(x, bytes) else str(x) for x in row).strip("\x00 -") for row in d.variables[k][:]]
nv, ev, gv = names("name_nod_var"), names("name_elem_var"), names("name_glo_var")
print("timesteps written", len(t), "of", d.dimensions["time_step"].size, "t:", t[:4], "...", t[-4:])
print("nodes", d.dimensions["num_nodes"].size, "elems", d.dimensions["num_elem"].size)
print("nodal", nv); print("global", gv); print("elem", len(ev), ev)
c = {ax: d.variables["coord" + ax][:] for ax in "xyz"}
print("extent", {ax: (float(c[ax].min()), float(c[ax].max())) for ax in "xyz"})
ax = a.pull_axis
L = a.length or float(c[ax].max() - c[ax].min())
sname = f"stress_{ax}{ax}"
if sname in ev and f"disp_{ax}" in nv:
    s = d.variables[f"vals_elem_var{ev.index(sname)+1}eb1"][:]
    u = d.variables[f"vals_nod_var{nv.index('disp_'+ax)+1}"][:]
    top = c[ax] == c[ax].max()
    print(f"step  t  eps_{ax}{ax}_nominal  mean_{sname}")
    idx = list(range(0, len(t), max(1, len(t) // 20))) + [len(t) - 1]
    for k in idx:
        print(k, round(float(t[k]), 4), round(float(u[k][top].mean()) / L, 5), round(float(s[k].mean()), 2))
    if a.csv:
        import csv
        with open(a.csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["step", "t", f"nominal_strain_{ax}{ax}", f"mean_{sname}"])
            for k in range(len(t)):
                stress_k = float(s[k].mean())
                if k > 0 and stress_k == 0.0:
                    break  # preallocated, unwritten step
                w.writerow([k + 1, float(t[k]), float(u[k][top].mean()) / L, stress_k])
        print("wrote", a.csv)
