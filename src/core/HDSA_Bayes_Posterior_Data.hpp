#ifndef HDSA_BAYES_POSTERIOR_DATA_HPP
#define HDSA_BAYES_POSTERIOR_DATA_HPP

namespace HDSA
{

  template <class RealT>
  class Bayes_Posterior_Data{

  public:
    RealT alpha;
    int N;
    HDSA::Ptr<HDSA::MultiVector<RealT> > Z;
    HDSA::Ptr<HDSA::MultiVector<RealT> > Y;
    HDSA::Ptr<HDSA::MultiVector<RealT> > Gamma_inv_Z;

    Bayes_Posterior_Data(void) 
    { }

    ~Bayes_Posterior_Data(void)
    { }
  
  };

}

#endif
