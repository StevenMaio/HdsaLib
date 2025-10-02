#ifndef HDSA_MD_POSTERIOR_VECTORS_HPP
#define HDSA_MD_POSTERIOR_VECTORS_HPP

namespace HDSA
{

  template <class RealT>
  class MD_Posterior_Vectors
  {

  public:
    int num_samples;
    HDSA::Ptr<HDSA::Vector<RealT>> mean;
    HDSA::Ptr<HDSA::MultiVector<RealT>> samples;

  public:
    MD_Posterior_Vectors(const int num_samples_in, const HDSA::Vector<RealT> &vec) : num_samples(num_samples_in)
    {
      mean = vec.Clone();
      samples = HDSA::makePtr<HDSA::MultiVector<RealT>>(num_samples_in, vec);
    }

    virtual ~MD_Posterior_Vectors()
    {
    }
  };

}

#endif
