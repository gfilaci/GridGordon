/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/modules/CPNenergy.h

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

#ifndef HMC_CPNENERGY_H
#define HMC_CPNENERGY_H

namespace Grid {
namespace QCD {

struct CPNEnergyParameters : Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(CPNEnergyParameters,
                                    double, beta)
    
    CPNEnergyParameters(double _a = 0.0):beta(_a){}
};
    
template <class Impl>
class CPNEnergyLogger : public HmcObservable<typename Impl::Field> {
    CPNEnergyParameters Pars;
public:
    
    // necessary for HmcObservable compatibility
    typedef typename Impl::Field Field;
    
    CPNEnergyLogger(double _a):Pars(_a){}
    
    CPNEnergyLogger(CPNEnergyParameters P):Pars(P){}
    
    void TrajectoryComplete(int traj, typename Impl::Field &U,
                            GridSerialRNG &sRNG,
                            GridParallelRNG &pRNG) {
        
        typename Impl::ZField zshifted(U._grid), ztmp(U._grid);
        decltype(peekSpin(U,0)) Umu(U._grid);
        
        // the definition of the CPN energy density can be found in
        // https://arxiv.org/pdf/1706.04443.pdf
        // or
        // https://arxiv.org/pdf/1508.07270.pdf
        
        typename Impl::ZField z = CPNObs<Impl>::extractZField(U);
        RealD res = 0;
        for (int mu=0; mu<Nd; mu++) {
            zshifted = Cshift(z,mu,1);
            Umu = peekSpin(U,mu);
            ztmp = conjugate(Umu) * z;
            res += real(innerProduct(zshifted,ztmp));
        }
        
        res = (double)Nd - res/(double)U._grid->gSites();
        
        int def_prec = std::cout.precision();
        
        std::cout << GridLogMessage
        << std::setprecision(std::numeric_limits<Real>::digits10 + 1)
        << "CPNEnergy : [ " << traj << " ] "<< res << std::endl;
        
        std::cout.precision(def_prec);
        
    }
};

}  // namespace QCD
}  // namespace Grid

#endif  // HMC_CPNENERGY_H
