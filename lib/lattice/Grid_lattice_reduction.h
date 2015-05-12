#ifndef GRID_LATTICE_REDUCTION_H
#define GRID_LATTICE_REDUCTION_H

namespace Grid {
#ifdef GRID_WARN_SUBOPTIMAL
#warning "Optimisation alert all these reduction loops are NOT threaded "
#endif     

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Deterministic Reduction operations
    ////////////////////////////////////////////////////////////////////////////////////////////////////
  template<class vobj> inline RealD norm2(const Lattice<vobj> &arg){
    ComplexD nrm = innerProduct(arg,arg);
    return real(nrm); 
  }

    template<class vobj>
    inline ComplexD innerProduct(const Lattice<vobj> &left,const Lattice<vobj> &right) 
    {
      typedef typename vobj::scalar_type scalar_type;
      typedef typename vobj::vector_type vector_type;
      vector_type vnrm;
      scalar_type  nrm;

      GridBase *grid = left._grid;

      std::vector<vector_type,alignedAllocator<vector_type> > sumarray(grid->SumArraySize());
      for(int i=0;i<grid->SumArraySize();i++){
	sumarray[i]=zero;
      }

PARALLEL_FOR_LOOP
      for(int thr=0;thr<grid->SumArraySize();thr++){

	int nwork, mywork, myoff;
	GridThread::GetWork(left._grid->oSites(),thr,mywork,myoff);
	
	decltype(innerProduct(left._odata[0],right._odata[0])) vnrm=zero; // private to thread; sub summation
        for(int ss=myoff;ss<mywork+myoff; ss++){
	  vnrm = vnrm + innerProduct(left._odata[ss],right._odata[ss]);
	}
	sumarray[thr]=TensorRemove(vnrm) ;
      }
    
      vector_type vvnrm; vvnrm=zero;  // sum across threads
      for(int i=0;i<grid->SumArraySize();i++){
	vvnrm = vvnrm+sumarray[i];
      } 
      nrm = Reduce(vvnrm);// sum across simd
      right._grid->GlobalSum(nrm);
      return nrm;
    }

    template<class vobj>
    inline typename vobj::scalar_object sum(const Lattice<vobj> &arg){

      GridBase *grid=arg._grid;
      int Nsimd = grid->Nsimd();

      std::vector<vobj,alignedAllocator<vobj> > sumarray(grid->SumArraySize());
      for(int i=0;i<grid->SumArraySize();i++){
	sumarray[i]=zero;
      }

PARALLEL_FOR_LOOP
      for(int thr=0;thr<grid->SumArraySize();thr++){
	int nwork, mywork, myoff;
	GridThread::GetWork(grid->oSites(),thr,mywork,myoff);

	vobj vvsum=zero;
        for(int ss=myoff;ss<mywork+myoff; ss++){
	  vvsum = vvsum + arg._odata[ss];
	}
	sumarray[thr]=vvsum;
      }

      vobj vsum=zero;  // sum across threads
      for(int i=0;i<grid->SumArraySize();i++){
	vsum = vsum+sumarray[i];
      } 

      typedef typename vobj::scalar_object sobj;
      sobj ssum=zero;

      std::vector<sobj>               buf(Nsimd);
      extract(vsum,buf);

      for(int i=0;i<Nsimd;i++) ssum = ssum + buf[i];
      arg._grid->GlobalSum(ssum);

      return ssum;
    }



template<class vobj> inline void sliceSum(const Lattice<vobj> &Data,std::vector<typename vobj::scalar_object> &result,int orthogdim)
{
  typedef typename vobj::scalar_object sobj;
  GridBase  *grid = Data._grid;
  assert(grid!=NULL);

  // FIXME
  std::cout<<"WARNING ! SliceSum is unthreaded "<<grid->SumArraySize()<<" threads "<<std::endl;

  const int    Nd = grid->_ndimension;
  const int Nsimd = grid->Nsimd();

  assert(orthogdim >= 0);
  assert(orthogdim < Nd);

  int fd=grid->_fdimensions[orthogdim];
  int ld=grid->_ldimensions[orthogdim];
  int rd=grid->_rdimensions[orthogdim];

  std::vector<vobj,alignedAllocator<vobj> > lvSum(rd); // will locally sum vectors first
  std::vector<sobj> lsSum(ld,zero); // sum across these down to scalars
  std::vector<sobj> extracted(Nsimd);     // splitting the SIMD

  result.resize(fd); // And then global sum to return the same vector to every node for IO to file
  for(int r=0;r<rd;r++){
    lvSum[r]=zero;
  }

  std::vector<int>  coor(Nd);  

  // sum over reduced dimension planes, breaking out orthog dir

  for(int ss=0;ss<grid->oSites();ss++){
    GridBase::CoorFromIndex(coor,ss,grid->_rdimensions);
    int r = coor[orthogdim];
    lvSum[r]=lvSum[r]+Data._odata[ss];
  }  

  // Sum across simd lanes in the plane, breaking out orthog dir.
  std::vector<int> icoor(Nd);

  for(int rt=0;rt<rd;rt++){

    extract(lvSum[rt],extracted);

    for(int idx=0;idx<Nsimd;idx++){

      grid->iCoorFromIindex(icoor,idx);

      int ldx =rt+icoor[orthogdim]*rd;

      lsSum[ldx]=lsSum[ldx]+extracted[idx];

    }
  }
  
  // sum over nodes.
  sobj gsum;
  for(int t=0;t<fd;t++){
    int pt = t/ld; // processor plane
    int lt = t%ld;
    if ( pt == grid->_processor_coor[orthogdim] ) {
      gsum=lsSum[lt];
    } else {
      gsum=zero;
    }

    grid->GlobalSum(gsum);

    result[t]=gsum;
  }

}


}
#endif

