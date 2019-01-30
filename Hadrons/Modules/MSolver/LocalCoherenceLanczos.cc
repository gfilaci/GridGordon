/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MSolver/LocalCoherenceLanczos.cc

Copyright (C) 2015-2018

Author: Antonin Portelli <antonin.portelli@me.com>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */
#include <Hadrons/Modules/MSolver/LocalCoherenceLanczos.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MSolver;

template class Grid::Hadrons::MSolver::TLocalCoherenceLanczos<FIMPL,HADRONS_DEFAULT_LANCZOS_NBASIS>;
template class Grid::Hadrons::MSolver::TLocalCoherenceLanczos<ZFIMPL,HADRONS_DEFAULT_LANCZOS_NBASIS>;
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
template class Grid::Hadrons::MSolver::TLocalCoherenceLanczos<FIMPL,HADRONS_DEFAULT_LANCZOS_NBASIS, FIMPLF>;
template class Grid::Hadrons::MSolver::TLocalCoherenceLanczos<ZFIMPL,HADRONS_DEFAULT_LANCZOS_NBASIS, ZFIMPLF>;
#endif