# J2 isotropic-hardening fit to the Vaughan (2024) Al2219 curve (digitized in
# Rezwan 2025 Fig. 3), for the PFD-on-J2 stage of the benchmark (docs/09 §5.1
# step 1). Elastic constants are the Voigt average of the CP study's cubic
# constants so that the J2 and CP plates share one elastic response; the fit
# is done in true stress and log strain because the J2 plate runs at finite
# strain (Hencky).
[Fit]
  curve       = ../../data/digitized/al2219_fig3_calibration_stress_strain.csv
  series      = experiment_Vaughan2024
  markers     = ../../data/digitized/al2219_fig3_experiment_markers.csv
  out_dir     = ../../results/pfd_study/j2_fit
  figure      = ../../figures/pfd/j2_fit.png
  C11         = 108.0e3
  C12         = 61.3e3
  C44         = 28.5e3
  measure     = true
  plastic_min = 5.0e-4
  hardening   = 'power piecewise'
  table_points = 24
[]
