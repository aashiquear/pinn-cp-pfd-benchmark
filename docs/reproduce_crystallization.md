# Reproducing Part IV (phase-field crystallization), stages A, B and C (first part)

PINNACLE `cd2288f` or later on `feat/pinn_benchmark` (ADR 0037 to 0039). Python 3.11 with numpy,
scipy and torch for PINNACLE; matplotlib and Pillow for the figures. All commands from the
repository root. Nothing here needs a GPU. Stage A takes about 10 minutes on 4 threads with
the guess scan, stage B about 35 minutes (30 of them the noise run).

## 1. A single case

Every deck under `decks/crystallization/single_field/` is a complete input:

```
python3.11 -m pinnacle run decks/crystallization/single_field/a3_ni_nucleus_x0.200_fem.i
```

It writes `results/crystallization/single_field/<deck name>/`: `solver_history.csv` (FE: one
row per step with `t`, `energy`, `mass`, Newton iterations), `postprocessors.csv` (the free
energy of the final field), `fields.vtk`, and the series `fields_step*.vtk` + `fields.pvd`.
The header of each deck states the physics, the units, the parameters with their provenance
and the expected result.

| case | decks | expected |
|---|---|---|
| A1 planar interface | `a1_planar_interface_{fem,pinn,hybrid}.i` | free energy 0.30200 aJ/nm^2, logistic profile of width 1.000 nm |
| A2 front velocity | `a2_front_velocity_{tanh,standard}_fem.i` | d(mass)/dt = 0.6723 and 0.6445 for t > 5 |
| A3 Ni nucleus | `a3_ni_nucleus_x<0.100 ... 0.300>_fem.i`, `a3_ni_nucleus_x0.200_{pinn,hybrid}.i` | W* in the deck header (0.91370 aJ = 47.87 kT at 0.200) |
| A4 other functions | `a4_ni_nucleus_<well>_<interpolation>_x<...>_fem.i` | W* in the deck header |
| A5 seeds | `a5_seed_{sub,super}critical_fem.i` | the smaller seed melts, the larger grows |

## 2. The whole stage

```
# reference solutions (no PINNACLE) and the digitized published curve
python3 scripts/cryst_reference.py decks/crystallization/reference/ni_nucleus_3d.i
python3 scripts/cryst_reference.py decks/crystallization/reference/ni_front.i
pdfimages -png -f 9 -l 10 data/literature/crystallization_review.pdf results/crystallization/digitize/rev
python3 scripts/cryst_digitize_fig3c.py          # needs the review PDF (not shipped)

# the A3 / A4 decks other than the hand-written x0.200 base (already in the repository)
python3 scripts/cryst_nucleus_decks.py

# runs
for d in decks/crystallization/single_field/*.i; do
  OMP_NUM_THREADS=4 python3.11 -m pinnacle run $d
done
python3 scripts/cryst_seed_runs.py 5             # network routes over five seeds

# tables, figures, animations
python3 scripts/cryst_nucleus_collect.py fem
python3 scripts/cryst_single_field_figures.py decks/crystallization/figures/single_field.i
```

Outputs: `results/crystallization/single_field/nucleus_fem_table.csv`, `front_fem_table.csv`,
`seeds/seed_table.csv`; figures and GIFs in `figures/crystallization/`.

## 3. Notes

- `[Functions] kind = parsed` expressions are written unquoted and without blanks.
- The stationary FE solve is one `bdf1` step with `dt = 1e12`; it converges to the nucleus
  from a droplet near it (the A3 and A4 decks carry a droplet fitted to the reference
  profile; at 0.200 the classical droplet also works).
- The digitized curve `data/digitized/cryst_review_fig3c_standard_pf.csv` is shipped; the
  raster it was read from is not.

## 4. Stage B

| case | decks | expected |
|---|---|---|
| B1 wall nucleus, Model I | `decks/crystallization/walls/b1_wall_model1_dT<010 ... 090>_fem.i` | W = 2 (free_energy + wall_energy); W / W_cap in the report's table |
| B1 wall nucleus, Model II | `b1_wall_model2_dT<010 ... 080>_{nucleus,layer}_fem.i` | W = 2 (F_nucleus - F_layer) |
| B2 noise and kinetics | `decks/crystallization/noise/b2_jmak_2d_fem.i` | Avrami exponent close to 3 (2.89 +/- 0.05 for the deck's seed) |

```
# stage A additions
python3 scripts/cryst_nucleus_collect.py pinn
python3 scripts/cryst_nucleus_collect.py hybrid
python3 scripts/cryst_guess_scan.py              # 27 runs, about 4 minutes

# B1: decks other than the hand-written Model I 40 K base (already in the repository)
python3 scripts/cryst_wall_decks.py
for d in decks/crystallization/walls/*.i; do
  OMP_NUM_THREADS=2 python3.11 -m pinnacle run $d
done
pdfimages -png -f 11 -l 11 data/literature/crystallization_review.pdf results/crystallization/digitize/p11
python3 scripts/cryst_digitize_fig6.py           # needs the review PDF (not shipped)
python3 scripts/cryst_wall_collect.py            # results/crystallization/walls/wall_table.csv
python3 scripts/cryst_wall_figures.py

# B2
OMP_NUM_THREADS=6 python3.11 -m pinnacle run decks/crystallization/noise/b2_jmak_2d_fem.i
python3 scripts/cryst_noise_figures.py           # kinetics csv, fit json, figures, GIF
```

Notes for stage B: the wall decks use `pf_fem_saddle` and start from the classical circular
cap; the Model II nucleus decks name a window function as `constraint_weight` (the wall layer
would otherwise be counted). The noise run writes 71 field files of 103,041 nodes (about
400 MB under `results/`).

## 5. Stage C, first part (binary model)

| case | decks | expected |
|---|---|---|
| C1 orientational free energy | none (`scripts/cryst_orientation_energy.py`) | maximum at 21.4 deg, minimum at 30.0 deg |
| C2 check | `decks/crystallization/binary/c2_check_iso_s1.0_fem.i` (200 x 200, 500 steps, about 1 min) | equivalent radius 100 cells (reference 93.0; 97.1 at half the spacing) |
| C2 Fig. 4 a, b, d, e | `decks/crystallization/binary/c2_fig4{a,b,d,e}_*_fem.i` (500 x 500; 26 to 42 min each) | growth forms of the report's section 5.3 |

```
python3 scripts/cryst_orientation_energy.py

# independent finite-difference reference (no PINNACLE), 2 to 8 min each
for d in decks/crystallization/reference/c2_*.i; do python3 scripts/cryst_binary_reference.py $d; done

# PINNACLE
for d in decks/crystallization/binary/c2_*_fem.i; do
  OMP_NUM_THREADS=3 python3.11 -m pinnacle run $d &
done; wait
for r in c2_fig4a_iso_s0.8_fem c2_fig4b_aniso_s0.8_fem c2_fig4d_iso_s1.0_fem c2_fig4e_aniso_s1.0_fem; do
  python3 scripts/cryst_binary_collect.py $r          # fields_step*.vtk -> frames.npz, history.csv
done
R=results/crystallization/binary
python3 scripts/cryst_binary_figures.py c2_pinnacle $R/c2_fig4a_iso_s0.8_fem $R/c2_fig4b_aniso_s0.8_fem \
    $R/c2_fig4d_iso_s1.0_fem $R/c2_fig4e_aniso_s1.0_fem
R=results/crystallization/binary_reference
python3 scripts/cryst_binary_figures.py c2_reference $R/fig4a_iso_s0.8 $R/fig4b_aniso_s0.8 \
    $R/fig4d_iso_s1.0 $R/fig4e_aniso_s1.0
```

Each 500 x 500 run writes 21 field files of 251,001 nodes (about 500 MB under `results/`).
The four C2 decks differ from the check deck in mesh size, supersaturation, anisotropy, noise
and step count only.

## 6. Stage C, orientation field (C3)

On the GPU host (Torch with ROCm; `TORCH_BLAS_PREFER_HIPBLASLT=0`):

```bash
for d in c3_fig4b_s0.8_chi360 c3_fig4c_s0.8_chi54 c3_fig4f_s1.0_chi54; do
  python scripts/cryst_orientation_reference.py decks/crystallization/reference/$d.i   # reference
  pinnacle run decks/crystallization/orientation/${d}_fem.i                              # PINNACLE
done
```

Then, where matplotlib is installed:

```bash
R=results/crystallization
for d in c3_fig4b_s0.8_chi360_fem c3_fig4c_s0.8_chi54_fem c3_fig4f_s1.0_chi54_fem; do
  python3 scripts/cryst_binary_collect.py orientation/$d
done
python3 scripts/cryst_orientation_figures.py c3_pinnacle $R/orientation/c3_fig4b_s0.8_chi360_fem \
  $R/orientation/c3_fig4c_s0.8_chi54_fem $R/orientation/c3_fig4f_s1.0_chi54_fem
python3 scripts/cryst_orientation_figures.py c3_reference $R/orientation_reference/fig4b_s0.8_chi360 \
  $R/orientation_reference/fig4c_s0.8_chi54 $R/orientation_reference/fig4f_s1.0_chi54
```

`use_gpu = false` in the deck runs the same case on the CPU. Expected numbers: report section 5.4.

## 7. Later cases (report sections 6.5.1, 7.4.4, 7.5, 8 to 11)

Every deck is complete: `pinnacle run <deck>` on the GPU host, or with `use_gpu = false`
(and `backend = numpy`) on a CPU. A case is a set of decks, then one figure script. The
PINNACLE revision, host and wall time of every set are in `docs/run_ledger.md` (rows
`cryst-*`). The parameter series of sections 7.5.1, 7.5.2 and 7.5.4 are run the same way;
their decks are in `decks/crystallization/orientation/{series,fibril,pinning}/` and the
`cryst_orientation_series.py` call of each figure is in its ledger row.

| case (report section) | decks under `decks/crystallization/` | deck generator | after the runs |
|---|---|---|---|
| E1 planar front (6.5.1) | `trapping/e1_planar_s*_fem.i` | `cryst_trapping_decks.py` | `cryst_trapping_figures.py` |
| C5 kinetics (7.4.4) | `orientation/kinetics/c5_n500_*_fem.i`; reference `reference/c5_kinetics_*.i` | none (the decks are shipped) | `cryst_binary_collect.py orientation/kinetics/<run>` per run, then `cryst_kinetics_figures.py <tag> [t_start=0.01425] <run dirs>` |
| C6 gallery (7.5.3) | `orientation/gallery/c6_fig10{a..j}_fem.i` | `cryst_orientation_decks.py` | `cryst_binary_collect.py orientation/gallery/<run>` per run, then `cryst_gallery_figures.py` |
| thermal histories (7.5.6) | `orientation/heat_treatment/c6_fig9{d_direct,e_short,e_mild,e_deep}_fem.i` | none (the three quench decks differ from the direct one in the two driving functions and the step count) | `cryst_binary_collect.py orientation/heat_treatment/<run>` per run, then `cryst_heat_treatment_figures.py` |
| E5 rough wall (7.5.5) | `walls_rough/e5_rough_wall_std*_fem.i` | mesh: `cryst_rough_wall_mesh.py` | `cryst_rough_wall_figures.py` |
| F1 hard spheres (8) | `hard_sphere/f1_nucleus_*_fem.i` | `cryst_hard_sphere_decks.py`; reference `cryst_hard_sphere_reference.py` | `cryst_hard_sphere_figures.py` |
| F4 quaternion field (9) | `quaternion/f4_*_fem.i` | `cryst_quaternion_decks.py spherulite n=64 S=1 chi=<chi> gamma=0.2 steps=<n>` | `cryst_quaternion_figures.py f4_spherulite max_fraction=0.3 <run dirs>` |
| F2 eutectic (10.4.1 to 10.4.3) | `eutectic/f2_lamellar_*_fem.i`, `f2_long100_*_fem.i` | `cryst_eutectic_decks.py <T> 0.345 <lambda ...> steps=300000 length=60` | `cryst_eutectic_figures.py` |
| F2b eutectic grain (10.4.4) | `eutectic/grain/f2b_grain_*_fem.i` | `cryst_eutectic_decks.py grain ...` | `cryst_eutectic_grain_figures.py` |
| F3 ternary (11) | `ternary/f3_ternary_*_fem.i`, `f3_long_*_fem.i` | `cryst_ternary_decks.py fs=<...> lam=<...> width=64 length=192 steps=<n>` | `cryst_ternary_figures.py`, `cryst_ternary_figures.py long` |

Notes:

- The quaternion decks need `energy_shape = linear` (PINNACLE `4bad49b` or later). With
  the default saturated shape a crystal does not order a liquid of random orientations.
- `cryst_binary_collect.py` takes one run directory per call, relative to
  `results/crystallization/`. It needs no matplotlib and can run on the GPU host; only
  `frames.npz` and the CSV files are then copied to the machine that makes the figures.
- The VTK series of the eutectic, ternary and rough-wall cases are read directly by their
  figure scripts (150 to 450 MB per run under `results/`).
- Gallery panel (e) stops with a non-finite field (report section 7.5.3).

## 8. Stage D: phase-field crystal (report section 12)

Needs PINNACLE with ADR 0046 (kernels `sh_potential`, `weak_laplacian`; `linear_solver = gmres_lu`).

| step | command | output |
|---|---|---|
| decks | `python3 scripts/cryst_pfc_decks.py` | `decks/crystallization/pfc/*.i` (26 decks) |
| reference: phase diagram | `python3 scripts/cryst_pfc_reference.py diagram` (about 10 min, one core) | `data/reference/cryst_pfc_phase_diagram.csv` |
| reference: one state | `python3 scripts/cryst_pfc_reference.py cell 0.3 -0.25` | free-energy density, peak height and cell size of the three phases (printed) |
| reference: growth | the command printed in the header of every `d2_growth_*.i` deck, e.g. `python3 scripts/cryst_pfc_reference.py grow 0.3 -0.34 26 15 12 1500 0.125 -0.505164` | `results/crystallization/pfc/reference_grow_*.npz` |
| runs | `pinnacle run decks/crystallization/pfc/<deck>.i` for every deck (two threads each are enough: `OMP_NUM_THREADS=2`) | `results/crystallization/pfc/<deck>/` |
| collect | `python3 scripts/cryst_pfc_collect.py` (NumPy only; can run on the compute host) | `<run>/frames.npz` |
| figures | `python3 scripts/cryst_pfc_figures.py all` | `figures/crystallization/d1_phase_diagram.png`, `d2_growth.png`, `d2_density_series.png`, `d2_eps_series.png`, `d3_polycrystal.png`, the GIFs, `results/crystallization/pfc/*_summary.csv` |

Notes.
- The growth decks solve a quarter of the crystal (mirror sides); the reference solves the full periodic box, which is the same problem.
- The seed amplitude is the one-mode amplitude at the solidus density of `data/reference/cryst_pfc_phase_diagram.csv`; the generator holds those densities as constants (`LIQUIDUS`, `SOLIDUS`), so regenerate the decks if the diagram is recomputed with other settings.
- The reference time step matters: the first-order semi-implicit scheme at dt = 0.5 overestimates the transformed fraction by about 3 % at t = 1500; the commands in the decks use 0.125.

## 9. Network routes of cases A2, A5 and B1 (report sections 3.4.2, 3.4.5 and 4.4)

Needs PINNACLE with ADR 0045 (`pf_pinn_stationary residual_form = fem`, `deflation_function`,
the per-step predictor of `pf_pinn_fem_hybrid`, `pf_pinn` on a spherical mesh) and Torch
(CPU is enough). The FE decks of stages A and B must have been run first: the wall decks
are derived from them and the figures compare with their results.

| step | command | output |
|---|---|---|
| wall decks | `python3 scripts/cryst_wall_network_decks.py` | `decks/crystallization/walls/b1_wall_*_{pinn,hybrid,hybrid150}_seed*.i`, `*_newton.i` (21 decks) |
| wall runs | `scripts/cryst_run_list.sh 2 results/crystallization/walls/runs_network.txt decks/crystallization/walls/b1_wall_*_newton.i decks/crystallization/walls/b1_wall_*_hybrid_seed*.i decks/crystallization/walls/b1_wall_*_pinn_seed*.i` and the same with `runs_network_20k.txt` for the `hybrid150` deck | one result directory per deck; 6.2 h of run time in total on 4 threads per run |
| wall table and figures | `python3 scripts/cryst_wall_network_figures.py` | `results/crystallization/walls/network_table.csv`, `figures/crystallization/b1_wall_network.png`, `b1_wall_newton_control.png` |
| Newton from saved fields | `python3 scripts/cryst_wall_network_newton_check.py` (about 2 min) | `results/crystallization/walls/newton_check.csv` |
| transient decks | `python3 scripts/cryst_transient_network_decks.py` | `decks/crystallization/single_field/{a2_front_velocity_standard,a5_seed_subcritical,a5_seed_supercritical}_{window_fem,pinn_seed{0,1,2},hybrid_seed0}.i` (15 decks) |
| transient runs | `scripts/cryst_run_list.sh 3 results/crystallization/single_field/runs_network_a2.txt decks/crystallization/single_field/a2_front_velocity_standard_{window_fem,hybrid_seed0,pinn_seed*}.i`; the same with `runs_network_a5.txt` and the `a5_seed_*` decks | 13 to 25 min per network run |
| transient table and figures | `python3 scripts/cryst_transient_network_figures.py` | `results/crystallization/single_field/network_transient_table.csv`, `figures/crystallization/a2_network.png`, `a5_network.png` |

Notes:

- The wall time of a run comes from the run list's log, so it depends on what else the
  machine does. Run at most 16 threads in total.
- A network run is reproducible for its seed on the same Torch version and thread count;
  the hybrid's corrected state does not depend on the seed when it reaches the nucleus.
- Model II, hybrid: check the end state. A field that is uniform along the wall is the wall
  layer, a thicker layer or the grown crystal, not the nucleus (report Table 17).
