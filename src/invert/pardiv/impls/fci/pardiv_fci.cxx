/************************************************************************
 *
 * FIXME
 *
 * Known issues:
 * ------------
 *
 * FIXME
 *
 *
 **************************************************************************
 * Copyright 2010 - 2022 BOUT++ contributors
 *
 * Contact: Ben Dudson, dudson2@llnl.gov
 *
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 ************************************************************************/

#include "pardiv_fci.hxx"
#include "bout/build_defines.hxx"

#if BOUT_USE_METRIC_3D and BOUT_HAS_PETSC

#include <bout/boutexception.hxx>
#include <bout/constants.hxx>
#include <bout/cyclic_reduction.hxx>
#include <bout/derivs.hxx>
#include <bout/fft.hxx>
#include <bout/globals.hxx>
#include <bout/surfaceiter.hxx>
#include <bout/utils.hxx>

#include <cmath>

InvertParDivFci::InvertParDivFci(Options* opt, CELL_LOC location, Mesh* mesh_in)
    : InvertParDiv(opt, location, mesh_in), ops(localmesh),
      // TODO: Neumann chosen, not validated against sheath BCs
      parallel(ops.getParallel()) {
    if (!localmesh->isFci()) {
        throw BoutException("Mesh is not FCI");
    }

  A = withParallelSlices(Field3D(1.0, localmesh));
  B = withParallelSlices(Field3D(0.0, localmesh));

}

Field3D InvertParDivFci::withParallelSlices(const Field3D& f) const {
  Field3D result = copy(f); // own data, so we don't modify the caller's field
  localmesh->communicate(result);
  result.applyParallelBoundary("parallel_neumann_o1");
  return result;
}

PetscCellOperator InvertParDivFci::assemble() const {
  // TODO: rebuilds the whole matrix every call; reuse structure, update values
  const auto& p = parallel;

  // B on the forward / backward field-line segments, as diagonal matrices
  const auto B_plus = ops.diagonal(p.Interp_plus(B));
  const auto B_minus = ops.diagonal(p.Interp_minus(B));

  // Div_par(B Grad_par), averaged over the two sides (as in Div_par_K_Grad_par)
  auto L = (p.Div_minus * B_plus * p.Grad_plus + p.Div_plus * B_minus * p.Grad_minus)
           * 0.5;

  return ops.diagonal(A) + L;
}

Field3D InvertParDivFci::solve(const Field3D& f) {
  const PetscCellOperator M = assemble();

  PetscInt rows = 0, cols = 0;
  BOUT_DO_PETSC(MatGetSize(M.raw(), &rows, &cols));
  output_info.write("InvertParDivFci: assembled {:d} x {:d} matrix\n", rows, cols);

  throw BoutException("fci backend: matrix assembled, solve not implemented");
}

#endif // BOUT_USE_METRIC_3D and BOUT_HAS_PETSC
