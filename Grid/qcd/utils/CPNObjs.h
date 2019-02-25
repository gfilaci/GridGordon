/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/qcd/utils/CPNObjs.h

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
#ifndef CPN_OBJS_H
#define CPN_OBJS_H
namespace Grid {

  // FIXME drop the QCD namespace in Nd
  

// Scalar field obs
template <class Impl>
class CPNObs {
 public:
  
  //////////////////////////////////////////////////
  // extract Gauge and Z fields
  //////////////////////////////////////////////////
  
  static typename Impl::Gauge extractGauge(const typename Impl::Field &f) {
      typename Impl::Gauge gaugeout(f._grid);
      for(int i=0; i<QCD::Nd; i++){
          auto tmp = QCD::peekSpin(f,i);
          QCD::pokeLorentz(gaugeout,tmp,i);
      }
      return gaugeout;
  }
  
  static typename Impl::ZField extractZField(const typename Impl::Field &f) {
      typename Impl::ZField zfieldout(f._grid);
      for(int i=QCD::Nd; i<Impl::Nfull; i++){
          auto tmp = QCD::peekSpin(f,i);
          QCD::pokeSpin(zfieldout,tmp,i-QCD::Nd);
      }
      return zfieldout;
  }

  //////////////////////////////////////////////////
  // load Gauge and Z fields
  //////////////////////////////////////////////////
    
  static typename Impl::Field loadGaugeZ(const typename Impl::Gauge &g, const typename Impl::ZField &z) {
      typename Impl::Field fieldout(g._grid);
      for(int i=0; i<QCD::Nd; i++){
          auto tmp = QCD::peekLorentz(g,i);
          QCD::pokeSpin(fieldout,tmp,i);
      }
      for(int i=QCD::Nd; i<Impl::Nfull; i++){
          auto tmp = QCD::peekSpin(z,i-QCD::Nd);
          QCD::pokeSpin(fieldout,tmp,i);
      }
      return fieldout;
  }

  //////////////////////////////////////////////////
  // project on CPN group
  //////////////////////////////////////////////////
    
  static typename Impl::Field ProjectOnCPN(const typename Impl::Field &f) {
      
      typename Impl::Field fieldout(f._grid);
      decltype(QCD::peekSpin(f,0)) tmp(f._grid), normtmp(f._grid);
      
      // project on U(1)
      for(int i=0; i<QCD::Nd; i++){
          tmp = QCD::peekSpin(f,i);
          normtmp = conjugate(tmp)*tmp;
          tmp = tmp / sqrt(normtmp);
          QCD::pokeSpin(fieldout,tmp,i);
      }
      
      // project on CP^N
      normtmp = zero;
      for(int i=QCD::Nd; i<Impl::Nfull; i++){
          tmp = QCD::peekSpin(f,i);
          normtmp = normtmp + conjugate(tmp)*tmp;
      }
      auto sqrtnorm = sqrt(normtmp);
      for(int i=QCD::Nd; i<Impl::Nfull; i++){
          tmp = QCD::peekSpin(f,i);
          tmp = tmp / sqrtnorm;
          QCD::pokeSpin(fieldout,tmp,i);
      }
      return fieldout;
  }

  //////////////////////////////////////////////////
  // project on space orthogonal to CPN
  //////////////////////////////////////////////////
  template <class T>
  static T ProjectOrthogonalCPN(const T &v, const T &z) {
      return v - outerProduct(z,z) * v;
  }
  static typename Impl::Field ProjectOrthogonalCPN(const typename Impl::Field &v, const typename Impl::Field &f) {
      auto U = CPNObs<Impl>::extractGauge(v);
      auto zv = CPNObs<Impl>::extractZField(v);
      auto zf = CPNObs<Impl>::extractZField(f);
      zv = ProjectOrthogonalCPN(zv,zf);
      return CPNObs<Impl>::loadGaugeZ(U, zv);
  }

};


}

#endif
