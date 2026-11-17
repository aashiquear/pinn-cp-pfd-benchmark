# Recheck of 2026-11-16 (recalibrated set, Bunge angles, total-Lagrangian balance, implicit update).
# Surrogate design for the CP+PFD plate at h40 (docs/09 §5.1 step 6, PINNACLE
# ADR 0030 data source): 24 LHS samples of the two fracture parameters and the
# two CP parameters that set yield and hardening, each a full coupled run of
# the calibrated base deck (scripts/pfd_study.py), targets the engineering
# stress and the maximum damage along the path (the fracture strain is the
# first crossing of 0.8) and the final damage field.
[Design]
  base       = ../base_h40.i
  driver     = pfd
  out_root   = ../../../../../results/pfd_study/cp/bunge/surrogate/h40_data
  samples    = 24
  sampler    = lhs
  seed       = 0
  concurrent = 8
  columns    = 'eng_stress max_damage'
  abscissa   = nominal_strain
  fields     = 'damage'
[]
[Parameters]
  [./Gc]
    path  = plate/Gc
    range = '561.33 2245.32'
    scale = log
  []
  [./l0]
    path  = plate/l0
    range = '4 10'
    scale = log
  []
  [./rho0]
    path  = plate/rho0
    range = '26.9653 107.861'
    scale = log
  []
  [./tau_c0]
    path  = plate/tau_c0
    range = '2.69478 10.7791'
    scale = log
  []
[]
