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
     Using action
     S = 1/(16 pi) (Dmu phi)^2 + 2 mu cosh(bphi)
     */

template <class Impl>
class shGordonAction : public QCD::Action<typename Impl::Field> {
 public:
    INHERIT_FIELD_TYPES(Impl);

 private:
    RealD mu_param;
    RealD b_param;
    const RealD inveightpi = 1. / 8. / M_PI;
    const RealD inv = 1. / 12. ;
    const RealD inv_t = 1. / 3. ;

 public:
    shGordonAction(RealD mu_, RealD b_) : mu_param(mu_), b_param(b_) {}

    virtual std::string LogParameters() {
      std::stringstream sstream;
      sstream << GridLogMessage << "[shGordonAction] mu     : " << mu_param << std::endl;
      sstream << GridLogMessage << "[shGordonAction] b      : " << b_param  << std::endl;
      return sstream.str();
    }
    virtual std::string action_name() {return "shGordonAction";}

    virtual void refresh(const Field &U, GridParallelRNG &pRNG) {}  // noop as no pseudoferms

    virtual RealD S(const Field &phi) {
//      return QCD::Nd * inveightpi * ScalarObs<Impl>::sumphisquared(phi) + inveightpi * ScalarObs<Impl>::sumphider(phi) + mu_param * sum(trace(exp(b_param*phi) + exp(-b_param*phi)))   ;
      return QCD::Nd * inveightpi * ScalarObs<Impl>::sumphisquared(phi) + inveightpi * ScalarObs<Impl>::sumphider(phi) + mu_param * b_param*b_param*ScalarObs<Impl>::sumphisquared(phi)
              + inv* mu_param*b_param*b_param*b_param*b_param*sum(trace(phi*phi*phi*phi)) ;

    };

    virtual void deriv(const Field &phi,
                       Field &force) {
        //std::cout << GridLogDebug << "Force total before :\n" << force << std::endl;
        Field tmp(phi._grid);
        tmp = 2.0 * QCD::Nd*phi;
        for (int mu = 0; mu < QCD::Nd; mu++) tmp -= Cshift(phi, mu, 1) + Cshift(phi, mu, -1);
        tmp = inveightpi * tmp;

        //std::cout << GridLogDebug << "Phi norm : " << norm2(phi) << std::endl;
        //force += tmp + b_param * mu_param * (exp(b_param*phi) - exp(-b_param*phi));
        force += tmp + 2.*b_param* b_param * mu_param * phi + inv_t*mu_param*b_param*b_param*b_param*b_param*phi*phi*phi ;
        //std::cout << GridLogDebug << "Force tmp :\n" << tmp << std::endl;
        //std::cout << GridLogDebug << "Force total after :\n" << force << std::endl;
    }
};



}  // namespace Grid

#endif // SHGORDON_ACTION_H
