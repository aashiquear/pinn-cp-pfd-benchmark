# Reproducing Part IV (phase-field crystallization), stage A

PINNACLE `ddaeb91` or later on `feat/pinn_benchmark` (ADR 0037). Python 3.11 with numpy,
scipy and torch for PINNACLE; matplotlib and Pillow for the figures. All commands from the
repository root. Nothing here needs a GPU; the whole stage takes about 5 minutes on 4 threads.

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
