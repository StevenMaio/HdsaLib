#ifndef PARAMETER_SAMPLER_SHALLOW_ICE_HPP
#define PARAMETER_SAMPLER_SHALLOW_ICE_HPP

template <class RealT>
class Parameter_Sampler_shallow_ice: public HDSA::Parameter_Sampler<RealT>{

private:
  int theta_dim_;

public:
  Parameter_Sampler_shallow_ice(const HDSA::Ptr<HDSA::ParameterList> & parlist): HDSA::Parameter_Sampler<RealT>() {
    int L_s = parlist->sublist("Problem").get("Number of Uncertain Basis Functions in Space", 10);
    int L_t = parlist->sublist("Problem").get("Number of Uncertain Basis Functions in Time", 10);
    theta_dim_ = 1 + (L_s+1)*(L_s+1) + (L_s+1)*(L_s+1)*(L_t+1);
  }

  ~Parameter_Sampler_shallow_ice() {}

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

}; // Parameter_Sampler_shallow_ice

#endif

