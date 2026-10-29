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

#include "pardiv_petsc.hxx"
#include "bout/build_defines.hxx"

#if BOUT_USE_METRIC_3D

#include <bout/boutexception.hxx>
#include <bout/constants.hxx>
#include <bout/cyclic_reduction.hxx>
#include <bout/derivs.hxx>
#include <bout/fft.hxx>
#include <bout/globals.hxx>
#include <bout/surfaceiter.hxx>
#include <bout/utils.hxx>

#include <cmath>

InvertParDivPetsc::InvertParDivPetsc(Options* opt, CELL_LOC location, Mesh* mesh_in)
    : InvertParDiv(opt, location, mesh_in) {
    abort();
}

Field3D InvertParDivPetsc::solve(const Field3D& f) {
    abort();
}

#endif // BOUT_USE_METRIC_3D
