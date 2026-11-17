# Figure deck of Part IV, stage A (single-field model). Consumed by
# scripts/cryst_single_field_figures.py; every path is relative to the repository root.
[Paths]
  results   = results/crystallization/single_field
  reference = results/crystallization/reference
  digitized = data/digitized/cryst_review_fig3c_standard_pf.csv
  figures   = figures/crystallization
[]
[Material]
  # Ni at the reference undercooling of the transient cases (Delta T / T_f = 0.2)
  delta_nm   = 0.2275598
  gamma      = 0.302
  delta_f    = -0.4923638
  w_star_aJ  = 0.9137052      # shooting solution
[]
[Animations]
  # run directory : frame duration in ms
  a1_planar_interface_fem   = 120
  a2_front_velocity_standard_fem = 80
  a5_seed_subcritical_fem   = 80
  a5_seed_supercritical_fem = 80
[]
