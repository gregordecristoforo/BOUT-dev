/************************************************************************
 * FIXME
 *
 * Known issues:
 * ------------
 *
 * FIXME
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

#if not BOUT_USE_METRIC_3D

namespace {
RegisterUnavailableInvertParDiv registerinvertpardivfci{
    PARDIVCYCLIC, "BOUT++ was configured without 3D metrics"};
}

#else

#include <bout/dcomplex.hxx>
#include <bout/globals.hxx>
#include <bout/utils.hxx>

class InvertParDivFci : public InvertParDiv {
public:
  explicit InvertParDivFci(Options* opt, CELL_LOC location = CELL_CENTRE,
                             Mesh* mesh_in = bout::globals::mesh);

  using InvertParDiv::solve;
  Field3D solve(const Field3D& f) override;

  using InvertParDiv::setCoefA;
  void setCoefA(const Field2D& f) override {
    ASSERT1(localmesh == f.getMesh());
    ASSERT1(location == f.getLocation());
    abort();
  }
  using InvertParDiv::setCoefB;
  void setCoefB(const Field2D& f) override {
    ASSERT1(localmesh == f.getMesh());
    ASSERT1(location == f.getLocation());
    abort();
  }

private:
  Field2D A{1.0}, B{0.0};

  int nsys;
};

namespace {
RegisterInvertParDiv<InvertParDivFci> registerinvertpardivfci{PARDIVFCI};
}

#endif // !BOUT_USE_METRIC_3D

#endif // INV_PARDIV_FCI_H
