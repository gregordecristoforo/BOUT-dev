bout_add_integrated_test(
  test_pardiv_fci
  SOURCES test_pardiv_fci.cxx
  USE_DATA_BOUT_INP
  REQUIRES BOUT_HAS_PETSC
  REQUIRES BOUT_USE_METRIC_3D
  EXTRA_FILES rotating-ellipse.py plot_grid.py
)
