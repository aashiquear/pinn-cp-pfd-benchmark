# Training targets of the Ni specimen surrogate (docs/12 §8, docs/13 P7; ADR 0030 addendum):
# the 30 synthetic specimens of decks/ni/final/tl/surrogate_design_h020.i, run at h020 by
# ni_study.py on results/ni_study/tl/final/surrogate/data/study.i, converted into the
# target.npz files PINNACLE's parametric_surrogate reads. No high-fidelity subset: the h040
# shortcut of P7.4 is not used (h040 runs 6-20% above h020 on the curve, docs/12 §5.2).
# Read by scripts/ni_surrogate_targets.py.
[Targets]
  design       = ../../../../../results/ni_study/tl/final/surrogate/data/design.csv
  runs         = ../../../../../results/ni_study/tl/final/surrogate/data/runs
  out          = ../../../../../results/ni_study/tl/final/surrogate/targets_h020
  abscissa     = gauge_strain
  columns      = eng_stress
  scalars      = 'failure_gauge_strain location_fraction peak_eng_stress'
  up_to        = failure
  orientations = ../../../../../data/extracted/ni_table1_orientations.csv
  thickness    = ../../../../../data/extracted/ni_table2_thickness.csv
[]
