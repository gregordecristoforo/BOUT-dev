/************************************************************************
 * FCI backend for InvertParDiv: solves (A + Div_par(B Grad_par)) x = b
 * with PETSc, using the support-operator matrices from PetscOperators.
 *
 * Derivation and notation: https://md.sigma2.no/tk9wifEUQXeUvgg177gIvQ?view
 * cited below as [notes §N].
 *
 * Known issues:
 * ------------
 *
 * - Solve not implemented yet [notes §4.6, §4.7]
 * - Boundary rows not set explicitly yet [notes §3.7, §4.4, §4.5]
 * - Only considered 1 rank
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

#ifndef INV_PARDIV_FCI_H
#define INV_PARDIV_FCI_H

#include "bout/build_defines.hxx"
#include "bout/invert_pardiv.hxx"

#if not BOUT_USE_METRIC_3D or not BOUT_HAS_PETSC

namespace {
RegisterUnavailableInvertParDiv registerinvertpardivfci{
    PARDIVFCI, "BOUT++ was configured without 3D metrics or PETSc"};
}

#else

#include <bout/dcomplex.hxx>
#include <bout/globals.hxx>
#include <bout/utils.hxx>
#include <bout/petsc_operators.hxx>

// TODO: KSP member + options from [pardiv] [notes §4.7]
/// Solves M x = b with M = diag(A) + Div_par(B Grad_par) on an FCI grid [notes §4]
class InvertParDivFci : public InvertParDiv {
public:
  explicit InvertParDivFci(Options* opt, CELL_LOC location = CELL_CENTRE,
                             Mesh* mesh_in = bout::globals::mesh);
  using InvertParDiv::solve;
  Field3D solve(const Field3D& f) override;

  // A, B stored as Field3D (no z-average), with parallel slices [notes §4.1, §4.2]
  using InvertParDiv::setCoefA;
  void setCoefA(const Field2D& f) override { setCoefA(Field3D(f)); }
  void setCoefA(const Field3D& f) override {
    ASSERT1(localmesh == f.getMesh());
    ASSERT1(location == f.getLocation());
    A = withParallelSlices(f);
  }

  using InvertParDiv::setCoefB;
  void setCoefB(const Field2D& f) override { setCoefB(Field3D(f)); }
  void setCoefB(const Field3D& f) override {
    ASSERT1(localmesh == f.getMesh());
    ASSERT1(location == f.getLocation());
    B = withParallelSlices(f);
  }

private:
  Field3D A{1.0}, B{0.0};

  // TODO: >1 rank needs PetscOperators column-layout fix (54a68824a)
  PetscOperators ops;  // reading FCI maps/weights from the grid: F, Bk, I±, w± [notes §3.2]
  PetscOperators::Parallel parallel; // Grad±, Div±, Interp± [notes §3.3, §3.4]

  /// Copy of f with communicated guard cells and parallel (yup/ydown) slices
  /// (Neumann copy into virtual cells) [notes §4.2]
  Field3D withParallelSlices(const Field3D& f) const;

  /// Build A + Div_par(B Grad_par) as a PETSc matrix [notes §3.5, §4.3]
  PetscCellOperator assemble() const;

};

namespace {
RegisterInvertParDiv<InvertParDivFci> registerinvertpardivfci{PARDIVFCI};
}

#endif // not BOUT_USE_METRIC_3D or not BOUT_HAS_PETSC

#endif // INV_PARDIV_FCI_H
