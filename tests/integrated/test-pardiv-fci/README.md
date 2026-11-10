# test-pardiv-fci

Standalone test of the `fci` InvertParDiv backend on a small FCI grid
(rotating ellipse, 20 x 8 x 20).

## Python scripts

Needs a Python environment with `zoidberg`, `netCDF4` and `matplotlib`,
e.g. `source ~/git/BOUT-dev/bout-venv/bin/activate`.

1. Create the grid (writes the PetscOperators metadata via `zoidberg.weights.modify_maps`):

       python rotating-ellipse.py -o grid.fci.nc

2. Plot the grid (cell types in several y planes, forward field-line map):

       python plot_grid.py grid.fci.nc           # show the figures
       python plot_grid.py grid.fci.nc --save    # write grid_planes.png, grid_forward_map.png

   Options: `--slices 0 3` selects the y planes to plot.

## Run the test

In the build directory, under `tests/integrated/test-pardiv-fci/`:

    make test_pardiv_fci
    python rotating-ellipse.py -o grid.fci.nc
    ./test_pardiv_fci

With the current skeleton backend it stops with "fci backend selected".
