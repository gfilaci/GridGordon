#ifndef CPN_IMPL
#define CPN_IMPL


namespace Grid {
  //namespace QCD {

template <class S>
class CPNImplTypes {
 public:
    typedef S Simd;

    template <typename vtype>
    using iImplField = iScalar<iScalar<iScalar<vtype> > >;

    typedef iImplField<Simd> SiteField;
    typedef SiteField        SitePropagator;
    typedef SiteField        SiteComplex;

    typedef Lattice<SiteField> Field;
    typedef Field              ComplexField;
    typedef Field              FermionField;
    typedef Field              PropagatorField;

    static inline void generate_momenta(Field& P, GridParallelRNG& pRNG){
      gaussian(pRNG, P);
    }

    static inline Field projectForce(Field& P){return P;}

    static inline void update_field(Field& P, Field& U, double ep) {
      //std::cout << GridLogDebug << "P:\n" << P << std::endl;
      U += P*ep;
      //std::cout << GridLogDebug << "U:\n" << U << std::endl;
    }

    static inline RealD FieldSquareNorm(Field& U) {
      return (- sum(trace(U*U))/2.0);
    }

    static inline void HotConfiguration(GridParallelRNG &pRNG, Field &U) {
     random(pRNG, U);
    }

    static inline void TepidConfiguration(GridParallelRNG &pRNG, Field &U) {
      random(pRNG, U);
      U *= 0.01;
    }

    static inline void ColdConfiguration(GridParallelRNG &pRNG, Field &U) {
      U = 0.0;
      //std::cout << GridLogDebug << "Initial U:\n" << U << std::endl;
    }

    static void MomentumSpacePropagator(Field &out, RealD m)
    {
      GridBase           *grid = out._grid;
      Field              kmu(grid), one(grid);
      const unsigned int nd    = grid->_ndimension;
      std::vector<int>   &l    = grid->_fdimensions;

      one = Complex(1.0,0.0);
      out = m*m;
      for(int mu = 0; mu < nd; mu++)
      {
        Real twoPiL = M_PI*2./l[mu];

        LatticeCoordinate(kmu,mu);
        kmu = 2.*sin(.5*twoPiL*kmu);
        out = out + kmu*kmu;
      }
      out = one/out;
    }

    static void FreePropagator(const Field &in, Field &out,
                               const Field &momKernel)
    {
      FFT   fft((GridCartesian *)in._grid);
      Field inFT(in._grid);

      fft.FFT_all_dim(inFT, in, FFT::forward);
      inFT = inFT*momKernel;
      fft.FFT_all_dim(out, inFT, FFT::backward);
    }

    static void FreePropagator(const Field &in, Field &out, RealD m)
    {
      Field momKernel(in._grid);

      MomentumSpacePropagator(momKernel, m);
      FreePropagator(in, out, momKernel);
    }

  };


  typedef CPNImplTypes<vComplex> CPNImplR;
  typedef CPNImplTypes<vComplexF> CPNImplF;
  typedef CPNImplTypes<vComplexD> CPNImplD;

}

#endif
