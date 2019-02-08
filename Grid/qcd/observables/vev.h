/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/modules/vev.h

Copyright (C) 2018

Author: Gianluca Filaci <g.filaci@ed.ac.uk>

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */

#ifndef HMC_VEV_H
#define HMC_VEV_H

namespace Grid {
namespace QCD {

template <class Impl>
class  VevLogger : public HmcObservable<typename Impl::Field> {
 public:

  // necessary for HmcObservable compatibility
  typedef typename Impl::Field Field;
  
  void TrajectoryComplete(int traj, typename Impl::Field &U,
                          GridSerialRNG &sRNG,
                          GridParallelRNG &pRNG) {

    double e = sum(trace(U)) / (double)U._grid->gSites();

    int def_prec = std::cout.precision();

    std::cout << GridLogMessage
        << std::setprecision(std::numeric_limits<Real>::digits10 + 1)
        << "Vev: [ " << traj << " ] "<< e << std::endl;

    std::cout.precision(def_prec);

  }
};

}  // namespace QCD
}  // namespace Grid

#endif  // HMC_VEV_H
