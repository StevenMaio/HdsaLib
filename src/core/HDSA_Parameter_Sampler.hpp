#ifndef HDSA_PARAMETER_SAMPLER_HPP
#define HDSA_PARAMETER_SAMPLER_HPP

// This class is where the distribution of the uncertain parameters is defined and sampled from.

namespace HDSA
{

  template <class RealT>
  class Parameter_Sampler{

  private:

    
  public:
    Parameter_Sampler() {}
    
    virtual ~Parameter_Sampler() {}
    
    virtual std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > Draw_Samples(int num_samp) = 0;
    
  }; // Parameter_Sampler
  
}

#endif

