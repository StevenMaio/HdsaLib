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
    HDSA::Ptr<HDSA::MultiVector<RealT> > Mz_inv_Gamma_inv_Z;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > G;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > g_vecs;
    HDSA::Ptr<HDSA::Vector<RealT> > Lambda;
    HDSA::Ptr<HDSA::MultiVector<RealT> > u_ell;
    std::vector<HDSA::Ptr<HDSA::MultiVector<RealT> > > u_i_ell;
    std::vector<HDSA::Ptr<HDSA::MultiVector<RealT> > > u_hat;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > a_ell;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > b_i_ell;

    Bayes_Posterior_Data(void) 
    { }

    ~Bayes_Posterior_Data(void)
    { }
  
  };

}

#endif
