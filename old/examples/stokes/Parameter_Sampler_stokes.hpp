#ifndef PARAMETER_SAMPLER_STOKES_HPP
#define PARAMETER_SAMPLER_STOKES_HPP

template <class RealT>
class Parameter_Sampler_stokes: public HDSA::Parameter_Sampler<RealT>{

private:
  int theta_dim_;

public:
  Parameter_Sampler_stokes(void): HDSA::Parameter_Sampler<RealT>() {
    theta_dim_ = 1;
  }

  ~Parameter_Sampler_stokes() {}

  std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > Draw_Samples(int num_samp)
  {
    std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > samples;
    samples.resize(num_samp);
    for(int k = 0; k < num_samp; k++)
      {
	samples[k] = HDSA::makePtr<Std_Vector<RealT> >(theta_dim_);
      }
    return samples;
  }

}; // Parameter_Sampler_stokes

#endif

