/*************************************************************************************

  Grid physics library, www.github.com/paboyle/Grid

  Source file: ./lib/qcd/action/gauge/shGordonAction.h

  Copyright (C) 2018

  Author: Guido Cossu <guido.cossu@ed.ac.uk>

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

#ifndef SHGORDON_ACTION_H
#define SHGORDON_ACTION_H

namespace Grid {

    /*
     The simulated action is
     
     S = 1/(16 pi) (Dmu phi)^2 + 2 mu cosh(b phi) + eta^2/(16 pi) phi^2
     
     With the redefinitions
     
     psi  = phi / sqrt(8 pi)
     g    = sqrt(8 pi) b
     m0^2 = 16 pi mu b
     
     it becomes
     
     S = 1/2 (Dmu psi)^2 + m0^2/g^2 cosh(g psi) + eta^2/2 psi^2 =
       = 1/2 (Dmu psi)^2 + (m0^2+eta^2)/2 psi^2 + ...
     
     so that the bare mass is M0^2 = m0^2+eta^2
     */
    
template <class Impl>
class shGordonAction : public QCD::Action<typename Impl::Field> {
 public:
    INHERIT_FIELD_TYPES(Impl);

 private:
    RealD mu_param;
    RealD b_param;
    RealD eta_param;
    const RealD inveightpi = 1. / 8. / M_PI;

 public:
    shGordonAction(RealD mu_, RealD b_, RealD eta_=0) : mu_param(mu_), b_param(b_), eta_param(eta_) {}

    virtual std::string LogParameters() {
      std::stringstream sstream;
      sstream << GridLogMessage << "[shGordonAction] mu     : " << mu_param  << std::endl;
      sstream << GridLogMessage << "[shGordonAction] b      : " << b_param   << std::endl;
      sstream << GridLogMessage << "[shGordonAction] eta    : " << eta_param << std::endl;
      return sstream.str();
    }
    virtual std::string action_name() {return "shGordonAction";}

    virtual void refresh(const Field &U, GridParallelRNG &pRNG) {}  // noop as no pseudoferms

    virtual RealD S(const Field &phi) {
      return (0.5*eta_param*eta_param + QCD::Nd) * inveightpi * ScalarObs<Impl>::sumphisquared(phi) + inveightpi * ScalarObs<Impl>::sumphider(phi) + mu_param * sum(trace(exp(b_param*phi) + exp(-b_param*phi)))  ;
    };

    virtual void deriv(const Field &phi,
                       Field &force) {
        //std::cout << GridLogDebug << "Force total before :\n" << force << std::endl;
        Field tmp(phi._grid);
        tmp = (eta_param*eta_param + 2.0*QCD::Nd)*phi;
        for (int mu = 0; mu < QCD::Nd; mu++) tmp -= Cshift(phi, mu, 1) + Cshift(phi, mu, -1);
        tmp = inveightpi * tmp;

        //std::cout << GridLogDebug << "Phi norm : " << norm2(phi) << std::endl;
        force += tmp + b_param * mu_param * (exp(b_param*phi) - exp(-b_param*phi));
        //std::cout << GridLogDebug << "Force tmp :\n" << tmp << std::endl;
        //std::cout << GridLogDebug << "Force total after :\n" << force << std::endl;
    }
};



}  // namespace Grid

#endif // SHGORDON_ACTION_H
