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
#include <bout/globals.hxx>
#include <bout/region.hxx>

InvertParDivFci::InvertParDivFci(Options* opt, CELL_LOC location, Mesh* mesh_in)
    : InvertParDiv(opt, location, mesh_in),
      lib(opt == nullptr ? &(Options::root()["pardiv"]) : opt) {
  if (!localmesh->isFci()) {
    throw BoutException("Mesh is not FCI");
  }
  A = Field3D{1.0, localmesh};
  B = Field3D{0.0, localmesh};

  // Discrete operators from the grid file (Grad±, Div±, Interp±, ...)
  ops.emplace(localmesh);
  par = ops->getParallel();

  // N: identity on interior and x-boundary rows,
  //    (N x)_virt = x_interior on virtual rows
  neumann = makeNeumannOperator(ops->getCellMapping(), BoundaryDirection::Both);

  // Default: GMRES + block Jacobi (works in serial and parallel).
  // Override in the input file, e.g. [pardiv:petsc] ksp_type = preonly, pc_type = lu
  BOUT_DO_PETSC(KSPCreate(BoutComm::get(), &ksp));
  BOUT_DO_PETSC(KSPSetType(ksp, KSPGMRES));
  PC pc;
  BOUT_DO_PETSC(KSPGetPC(ksp, &pc));
  BOUT_DO_PETSC(PCSetType(pc, PCBJACOBI));
  lib.setOptionsFromInputFile(ksp);
}

InvertParDivFci::~InvertParDivFci() {
  if (ksp != nullptr) {
    KSPDestroy(&ksp);
  }
}

Field3D InvertParDivFci::withSlices(const Field3D& f, BoutReal virt) const {
  Field3D result = f;
  result.clearParallelSlices();
  result.splitParallelSlices();
  result.yup() = virt;
  result.ydown() = virt;
  return result;
}

void InvertParDivFci::assemble() {
  // Coefficients with values on virtual cells
  // A: irrelevant there (rows are replaced) -> 0
  // B: Neumann copy of the adjacent interior cell, needed for Interp±(B)
  const Field3D A_c = withSlices(A, 0.0);
  Field3D B_c = withSlices(B, 0.0);
  B_c.applyParallelBoundary("parallel_neumann_o1");

  // M_int = diag(A) + ½ [ Div− diag(Interp+ B) Grad+ + Div+ diag(Interp− B) Grad− ]
  auto L_plus = par.Div_minus * ops->diagonal(par.Interp_plus(B_c)) * par.Grad_plus;
  auto L_minus = par.Div_plus * ops->diagonal(par.Interp_minus(B_c)) * par.Grad_minus;
  auto M_int = ops->diagonal(A_c) + (L_plus + L_minus) * 0.5;

  // Row masks (0/1 on cell space C)
  Field3D int_main{0.0, localmesh};
  BOUT_FOR(i, int_main.getRegion("RGN_NOBNDRY")) { int_main[i] = 1.0; }

  const Field3D m_int = withSlices(int_main, 0.0);                 // interior cells
  const Field3D m_xb = withSlices(1.0 - int_main, 0.0);            // x-boundary cells
  const Field3D m_virt = withSlices(Field3D{0.0, localmesh}, 1.0); // virtual cells
  const auto I = ops->diagonal(withSlices(Field3D{1.0, localmesh}, 1.0));

  // M = diag(m_int) M_int + diag(m_virt) (I − N) + diag(m_xb)
  M = ops->diagonal(m_int) * M_int + ops->diagonal(m_virt) * (I - *neumann)
      + ops->diagonal(m_xb);

  BOUT_DO_PETSC(KSPSetOperators(ksp, M->raw(), M->raw()));
}

Field3D InvertParDivFci::solve(const Field3D& f) {
  ASSERT1(localmesh == f.getMesh());
  ASSERT1(location == f.getLocation());

  if (dirty) {
    assemble();
    dirty = false;
  }

  // b = v on interior and x-boundary cells, 0 on virtual cells
  const auto b = makePetscCellVector(ops->getCellMapping(), withSlices(f, 0.0));
  PetscCellVector x(ops->getCellMapping());

  BOUT_DO_PETSC(KSPSolve(ksp, b.raw(), x.raw()));

  KSPConvergedReason reason;
  BOUT_DO_PETSC(KSPGetConvergedReason(ksp, &reason));
  if (reason < 0) {
    throw BoutException("InvertParDivFci: KSP diverged, reason {:d}",
                        static_cast<int>(reason));
  }

  // Back to a Field3D; virtual-cell values are dropped
  Field3D result = toField3D(x, f);
  result.clearParallelSlices();
  return result;
}

#endif // BOUT_USE_METRIC_3D and BOUT_HAS_PETSC
