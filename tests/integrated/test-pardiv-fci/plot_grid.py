#!/usr/bin/env python
"""Plot the FCI grid written by rotating-ellipse.py

Figure 1: the poloidal (R, Z) grid in several y planes, with the cell types
          (interior, x boundary, and interior cells whose forward field line
          hits the wall, i.e. that have a forward boundary leg / virtual cell).
Figure 2: forward field-line map from one plane to the next: where the field
          line from each cell centre lands in the next plane.

Usage:
    python rotating-ellipse.py -o grid.fci.nc
    python plot_grid.py grid.fci.nc [--save]
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset

parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument("gridfile", nargs="?", default="grid.fci.nc")
parser.add_argument(
    "--slices", type=int, nargs="+", default=None, help="y indices to plot"
)
parser.add_argument("--save", action="store_true", help="save PNGs instead of showing")
args = parser.parse_args()

with Dataset(args.gridfile) as d:
    R = np.asarray(d["R"][:])  # (x, y, z) cell centres
    Z = np.asarray(d["Z"][:])
    mxg = int(d["MXG"][:])
    fwd_R = np.asarray(d["forward_R"][:])  # end of the forward leg in plane y+1
    fwd_Z = np.asarray(d["forward_Z"][:])
    fwd_bnd = np.asarray(d["forward_leg_boundary_number"][:])  # >= 0: hits wall
    bwd_bnd = np.asarray(d["backward_leg_boundary_number"][:])
    total_cells = int(d["total_cells"][:])
    n_fwd_legs = int(d["n_forward_legs"][:])
    n_bwd_legs = int(d["n_backward_legs"][:])

nx, ny, nz = R.shape
interior = np.zeros(R.shape, dtype=bool)
interior[mxg : nx - mxg] = True

n_int = interior.sum()
n_xb = (~interior).sum()
print(f"Grid {args.gridfile}: nx={nx} (incl. {2 * mxg} x-boundary), ny={ny}, nz={nz}")
print(f"  interior cells:   {n_int}")
print(f"  x-boundary cells: {n_xb}")
print(
    f"  virtual cells:    {total_cells - n_int - n_xb} "
    f"(forward {(fwd_bnd >= 0).sum()}, backward {(bwd_bnd >= 0).sum()})"
)
print(f"  forward legs: {n_fwd_legs}, backward legs: {n_bwd_legs}")


def draw_mesh(ax, j, **kw):
    """Grid lines through the cell centres of plane j (z is periodic)"""
    Rc = np.concatenate([R[:, j, :], R[:, j, :1]], axis=1)
    Zc = np.concatenate([Z[:, j, :], Z[:, j, :1]], axis=1)
    ax.plot(Rc, Zc, **kw)  # lines of constant x
    ax.plot(Rc.T, Zc.T, **kw)  # lines of constant z


# Figure 1: poloidal planes and cell types
slices = args.slices if args.slices is not None else list(range(0, ny, max(1, ny // 4)))
fig1, axs = plt.subplots(1, len(slices), figsize=(4 * len(slices), 4.4), squeeze=False)
for ax, j in zip(axs[0], slices):
    draw_mesh(ax, j, color="0.8", lw=0.5)
    xb = ~interior[:, j, :]
    hits = interior[:, j, :] & ((fwd_bnd[:, j, :] >= 0) | (bwd_bnd[:, j, :] >= 0))
    rest = interior[:, j, :] & ~hits
    ax.scatter(R[:, j, :][rest], Z[:, j, :][rest], s=6, c="C0", label="interior")
    ax.scatter(R[:, j, :][hits], Z[:, j, :][hits], s=12, c="C3",
               label="interior, field line hits wall")
    ax.scatter(R[:, j, :][xb], Z[:, j, :][xb], s=6, c="C2", marker="s",
               label="x boundary")
    ax.set_title(f"y index {j}")
    ax.set_aspect("equal")
    ax.set_xlabel("R [m]")
axs[0, 0].set_ylabel("Z [m]")
axs[0, 0].legend(loc="lower left", fontsize=7)
fig1.suptitle("Poloidal planes and cell types")
fig1.tight_layout()

# Figure 2: forward map from plane j to plane j+1
j = slices[0]
jn = (j + 1) % ny
fig2, ax = plt.subplots(figsize=(6, 6))
draw_mesh(ax, jn, color="0.75", lw=0.5)
ax.plot([], [], color="0.75", lw=0.5, label=f"grid of plane {jn}")
inner = interior[:, j, :] & (fwd_bnd[:, j, :] < 0)
wall = interior[:, j, :] & (fwd_bnd[:, j, :] >= 0)
ax.scatter(R[:, j, :][interior[:, j, :]], Z[:, j, :][interior[:, j, :]], s=5, c="C0",
           label=f"cell centres, plane {j}")
ax.scatter(fwd_R[:, j, :][inner], fwd_Z[:, j, :][inner], s=8, c="C1", marker="x",
           label=f"field line lands in plane {jn}")
ax.scatter(fwd_R[:, j, :][wall], fwd_Z[:, j, :][wall], s=12, c="C3", marker="x",
           label="field line hits the wall (boundary leg)")
ax.set_aspect("equal")
ax.set_xlabel("R [m]")
ax.set_ylabel("Z [m]")
ax.set_title(f"Forward map, plane {j} → {jn}")
ax.legend(loc="lower left", fontsize=7)
fig2.tight_layout()

if args.save:
    fig1.savefig("grid_planes.png", dpi=150)
    fig2.savefig("grid_forward_map.png", dpi=150)
    print("Saved grid_planes.png, grid_forward_map.png")
else:
    plt.show()
