/*************************************************************************************
 
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: ./lib/qcd/action/CPN/CPNImpl.h
 
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

#ifndef CPN_IMPL
#define CPN_IMPL


namespace Grid {
  //namespace QCD {

template <int N, class S>
class CPNImplTypes {
 public:
    typedef S Simd;
    
    static const int Nfull = N + QCD::Nd;
    static const int NCPN = N;
    
    template <typename vtype>
    using iImplField = iScalar<iVector <iScalar<vtype>, Nfull > >;

    typedef iImplField<Simd> SiteField;
    typedef SiteField        SitePropagator;
    typedef SiteField        SiteComplex;

    typedef Lattice<SiteField> Field;
    typedef Field              ComplexField;
    typedef Field              FermionField;
    typedef Field              PropagatorField;

    template <typename vtype>
    using iImplGauge = iVector<iScalar<iScalar<vtype> >, QCD::Nd >;
    typedef iImplGauge<Simd> SiteGauge;
    typedef Lattice<SiteGauge> Gauge;
    
    template <typename vtype>
    using iImplZField = iScalar<iVector<iScalar<vtype>, N> >;
    typedef iImplZField<Simd> SiteZField;
    typedef Lattice<SiteZField> ZField;
    
    static inline void generate_momenta(Field& P, GridParallelRNG& pRNG){
        gaussian(pRNG, P);
        auto Pg = CPNObs<CPNImplTypes>::extractGauge(P);
        auto Pz = CPNObs<CPNImplTypes>::extractZField(P);
        // momenta conjugated to the gauge field must be real
        Pg += conjugate(Pg);
        Pg = 0.5 * Pg;
        P = CPNObs<CPNImplTypes>::loadGaugeZ(Pg,Pz);
    }

    static inline void project_momenta(Field& P, const Field& U){
        P = CPNObs<CPNImplTypes>::ProjectOrthogonalCPN(P,U);
    }

    static inline Field projectForce(Field& P){return P;}

    static inline void update_field(Field& P, Field& U, double ep) {

        decltype(QCD::peekSpin(U,0)) Ptmp(U._grid), Pcos(U._grid), Psin(U._grid), Pabs(U._grid), Utmp(U._grid), Uold(U._grid);
        
        Complex im(0.,1.);
        for(int i=0; i<QCD::Nd; i++){
            Ptmp = QCD::peekSpin(P,i);
            Utmp = QCD::peekSpin(U,i);
            Utmp = exp(im*ep*Ptmp) * Utmp;
            QCD::pokeSpin(U,Utmp,i);
        }
        
        Pabs = zero;
        for(int i=QCD::Nd; i<Nfull; i++){
            Ptmp = QCD::peekSpin(P,i);
            Pabs += conjugate(Ptmp)*Ptmp;
        }
        Pabs = sqrt(Pabs);
        Psin = sin(ep*Pabs);
        Pcos = cos(ep*Pabs);
        
        for(int i=QCD::Nd; i<Nfull; i++){
            Ptmp = QCD::peekSpin(P,i);
            Uold = QCD::peekSpin(U,i);
            
            // update
            Utmp =  Pcos * Uold        + Psin * Ptmp / Pabs;
            Ptmp = -Pabs * Psin * Uold + Pcos * Ptmp;
            
            QCD::pokeSpin(U,Utmp,i);
            QCD::pokeSpin(P,Ptmp,i);
        }
        
        // UPDATE CHECK
//        auto Ug = CPNObs<CPNImplTypes>::extractGauge(U);
//        auto Uz = CPNObs<CPNImplTypes>::extractZField(U);
//        auto Pg = CPNObs<CPNImplTypes>::extractGauge(P);
//        auto Pz = CPNObs<CPNImplTypes>::extractZField(P);
//        std::cout << "U belongs to U(1)     : " << norm2(Ug)/(double)QCD::Nd/(double)U._grid->gSites() - 1. << std::endl;
//        std::cout << "z belongs to CPN      : " << norm2(Uz)/(double)U._grid->gSites() - 1. << std::endl;
//        Pg = 0.5 * ( Pg - conjugate(Pg) );
//        std::cout << "PU is real            : " << norm2(Pg)/(double)QCD::Nd/(double)U._grid->gSites() << std::endl;
//        std::cout << "Pz is orthogonal to z : " << innerProduct(Uz,Pz)/(double)U._grid->gSites() << std::endl;
    }

    static inline RealD FieldSquareNorm(Field& U) {
        RealD res = 0;
        decltype(QCD::peekSpin(U,0)) tmp(U._grid);
        for(int i=0; i<Nfull; i++){
            tmp = QCD::peekSpin(U,i);
            res += norm2(tmp);
        }
      return -0.5*res;
    }

    static inline void HotConfiguration(GridParallelRNG &pRNG, Field &U) {
        random(pRNG, U);
        U = CPNObs<CPNImplTypes>::ProjectOnCPN(U);
    }

    static inline void TepidConfiguration(GridParallelRNG &pRNG, Field &U) {
      random(pRNG, U);
      U *= 0.01;
      U = CPNObs<CPNImplTypes>::ProjectOnCPN(U);
    }

    static inline void ColdConfiguration(GridParallelRNG &pRNG, Field &U) {
      decltype(QCD::peekSpin(U,0)) tmp(U._grid);
      tmp = 1.;
      for(int i=0; i<QCD::Nd+1; i++){
          QCD::pokeSpin(U,tmp,i);
      }
      tmp = zero;
      for(int i=QCD::Nd+1; i<Nfull; i++){
          QCD::pokeSpin(U,tmp,i);
      }
    }

    static void MomentumSpacePropagator(Field &out, RealD m)
    {
        std::cout << GridLogError << "MomentumSpacePropagator not implemented for CPNImplTypes" << std::endl;
        assert(0);
    }

    static void FreePropagator(const Field &in, Field &out,
                               const Field &momKernel)
    {
        std::cout << GridLogError << "FreePropagator not implemented for CPNImplTypes" << std::endl;
        assert(0);
    }

    static void FreePropagator(const Field &in, Field &out, RealD m)
    {
        std::cout << GridLogError << "FreePropagator not implemented for CPNImplTypes" << std::endl;
        assert(0);
    }

  };

  template<int N> using CPNImplR = CPNImplTypes<N, vComplex>;
  template<int N> using CPNImplF = CPNImplTypes<N, vComplexF>;
  template<int N> using CPNImplD = CPNImplTypes<N, vComplexD>;
  
}

#endif
