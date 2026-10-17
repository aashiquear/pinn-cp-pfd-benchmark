# Ni Part III stage 3: mesh sensitivity on specimen L.1 with the calibrated CP set
# (docs/11 §5; decks/ni/cp/calibrated_L1_h040.i). The h040 reference is the calibrated
# run itself (results/ni_study/cal/gn4/base, linked as <out_root>/base); the cases change
# only the mesh: gage cell 0.2 and 0.1 mm with one element through the thickness, and
# 0.2 mm with two. Grip strain 30% in 300 steps (the calibration's 1e-3 increment and
# 5e-2/s bar rate; the curves end at the first cutback, 54-60% gauge strain at h040).
# Threads per case: 2 / 8 / 4 (14 of the host's 16 cores), all side by side.
[Study]
  base       = calibrated_L1_h040.i
  out_root   = ../../../results/ni_study/mesh_study/L1
  markers    = ../../../data/digitized/ni_fig3_experiment_markers.csv
  concurrent = 3
[]
[Cases]
  [./h020]
    mesh    = ../../../results/ni_study/mesh/h020/L.1/mesh.e
    strain  = 0.30
    t_end   = 6.0
    n_steps = 300
    threads = 2
  []
  [./h010]
    mesh    = ../../../results/ni_study/mesh/h010/L.1/mesh.e
    strain  = 0.30
    t_end   = 6.0
    n_steps = 300
    threads = 8
  []
  [./h020z2]
    mesh    = ../../../results/ni_study/mesh/h020z2/L.1/mesh.e
    strain  = 0.30
    t_end   = 6.0
    n_steps = 300
    threads = 4
  []
[]
