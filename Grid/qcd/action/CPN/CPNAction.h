/*************************************************************************************

  Grid physics library, www.github.com/paboyle/Grid

  Source file: ./lib/qcd/action/CPN/CPNAction.h

  Copyright (C) 2015

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

#ifndef CPN_ACTION_H
#define CPN_ACTION_H

namespace Grid {
  // FIXME drop the QCD namespace everywhere here

template <class Impl>
class CPNAction : public QCD::Action<typename Impl::Field> {
 public:
    INHERIT_FIELD_TYPES(Impl);

 private:
    RealD beta;
    RealD factor;
    
 public:
    CPNAction(RealD beta_) : beta(beta_) {
        factor = Impl::NCPN * beta_;
    }

    virtual std::string LogParameters() {
      std::stringstream sstream;
      sstream << GridLogMessage << "[CPNAction] beta : " << beta << std::endl;
      return sstream.str();
    }
    virtual std::string action_name() {return "CPNAction";}

    virtual void refresh(const Field &U, GridParallelRNG &pRNG) {}  // noop as no pseudoferms

    virtual RealD S(const Field &p) {
        typename Impl::ZField zshifted(p._grid);
        decltype(QCD::peekSpin(p,0)) Umu(p._grid);
        
        typename Impl::ZField z = CPNObs<Impl>::extractZField(p);
        RealD res = 0;
        for (int mu=0; mu<QCD::Nd; mu++) {
            zshifted = Cshift(z,mu,1);
            Umu = QCD::peekSpin(p,mu);
            z = conjugate(Umu) * z;
            res += real(innerProduct(zshifted,z));
        }
        return (-2.*factor)*res;
    };

    virtual void deriv(const Field &p,
                       Field &force) {
        typename Impl::ZField zshifted(p._grid), Fz(p._grid);
        typename Impl::Gauge Fg(p._grid);
        decltype(QCD::peekSpin(p,0)) Umu(p._grid);
        
        typename Impl::ZField z = CPNObs<Impl>::extractZField(p);
        
        // FIXME: inefficent, both loops can be reduced to one
        
        Fz = zero;
        for (int mu=0; mu<QCD::Nd; mu++) {
            Umu = QCD::peekSpin(p,mu);
            zshifted = Cshift(z,mu,1);
            Fz += Umu * zshifted;
            
            zshifted = Cshift(z,mu,-1);
            Umu = Cshift(Umu,mu,-1);
            Fz += conjugate(Umu) * zshifted;
        }
        Fz = (-factor)*conjugate(Fz);
        
        for (int mu=0; mu<QCD::Nd; mu++) {
            zshifted = Cshift(z,mu,1);
            Umu = QCD::peekSpin(p,mu);
            zshifted = Umu * zshifted;
            
            // perform scalar product for each lattice site
            decltype(QCD::peekSpin(p,0)) sp(p._grid);
            sp = zero;
            for (int i=0; i<Impl::NCPN; i++) {
                auto tmpsp1 = conjugate(QCD::peekSpin(z,i));
                auto tmpsp2 = QCD::peekSpin(zshifted,i);
                sp += tmpsp1*tmpsp2;
            }
            
            Umu = (2.*factor)*imag(sp);
            QCD::pokeLorentz(Fg,Umu,mu);
        }
        force = CPNObs<Impl>::loadGaugeZ(Fg, Fz);
    }
};



}  // namespace Grid

#endif // CPN_ACTION_H
