# Equivalent plastic strain on the deformed plate for the four meshes of the
# mesh study, animated over the loading (author's request 2026-10-11). Read by
# scripts/plastic_strain_animation.py. The series come from the field reruns
# decks/cp_study/fields_h{20,40,80,160}.i (snapshots_every = 4: 64 frames to
# 20% strain).
[Animation]
  out_gif    = ../../figures/report/plastic_strain_mesh.gif
  out_stills = ../../figures/report/plastic_strain_mesh_{strain}pct.png
  stills     = '2 8 15 20'
  strain_max = 0.20
  vmax       = 0.6
  scale      = 1.0
  fps        = 8
[]
[Series]
  [./h20]
    dir   = ../../results/cp_study/fields/h20/solver
    label = 'h20, h = 5 um'
  []
  [./h40]
    dir   = ../../results/cp_study/fields/h40/solver
    label = 'h40, h = 2.5 um'
  []
  [./h80]
    dir   = ../../results/cp_study/fields/h80/solver
    label = 'h80, h = 1.25 um'
  []
  [./h160]
    dir   = ../../results/cp_study/fields/h160/solver
    label = 'h160, h = 0.625 um'
  []
[]
