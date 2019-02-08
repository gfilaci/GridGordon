/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/modules/two_point.h

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

#ifndef HMC_TWO_POINT_H
#define HMC_TWO_POINT_H

namespace Grid {
namespace QCD {

struct TwoPointParameters : Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(TwoPointParameters,
    std::string, T)
    
    std::vector<int> Tvec;
    std::vector<std::string> Tnames;
    
    TwoPointParameters(int N = 0) {
        Tvec.resize(N);
        for(int i=0; i<N; i++) Tvec[i] = i;
        set_names();
    }
    
    TwoPointParameters(std::vector<int> Tvec_): Tvec(Tvec_) { set_names(); }
    
    template < class ReaderClass >
    TwoPointParameters(Reader<ReaderClass>& Reader){
        read(Reader, "TwoPoint", *this);
        Tvec = strToVec<int>(T);
        set_names();
    }
    
    void set_names(){
        Tnames.resize(Tvec.size());
        for(int i=0; i<Tvec.size(); i++){
            std::ostringstream stream;
            stream << Tvec[i];
            Tnames[i] = stream.str();
        }
    }
};

template <class Impl>
class TwoPointLogger : public HmcObservable<typename Impl::Field> {
  TwoPointParameters Pars;
 public:

  // necessary for HmcObservable compatibility
  typedef typename Impl::Field Field;
  typedef typename Impl::SiteField::scalar_object SiteField;
  
  TwoPointLogger(double _a):Pars(_a){}

  TwoPointLogger(TwoPointParameters P):Pars(P){}

  void TrajectoryComplete(int traj, typename Impl::Field &U,
                          GridSerialRNG &sRNG,
                          GridParallelRNG &pRNG) {
    
    // set which direction is the temporal direction
    const int tmpdir = 0;
    int Ltmp = U._grid->_fdimensions[tmpdir];
      
    // sum over all directions except for the temporal one
    std::vector<SiteField> sumresult;
    sliceSum(U, sumresult, tmpdir);
    
    int def_prec = std::cout.precision();
    
    for(int i=0; i<Pars.Tvec.size(); i++){
        
        SiteField corr = 0.;
        for(int j=0; j<Ltmp; j++) corr += sumresult[j] * sumresult[(j+Pars.Tvec[i])%Ltmp];
        double obs = trace(corr) / (double)U._grid->gSites();
        
        std::cout << GridLogMessage
        << std::setprecision(std::numeric_limits<Real>::digits10 + 1)
        << "TwoPoint" << Pars.Tnames[i] << ": [ " << traj << " ] "<< obs << std::endl;
        
        std::cout.precision(def_prec);
    }
  }

};

}  // namespace QCD
}  // namespace Grid

#endif  // HMC_TWO_POINT_H
