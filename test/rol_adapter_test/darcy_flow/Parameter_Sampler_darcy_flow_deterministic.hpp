#ifndef PARAMETER_SAMPLER_DARCY_FLOW_DETERMINISTIC_HPP
#define PARAMETER_SAMPLER_DARCY_FLOW_DETERMINISTIC_HPP

template <class RealT>
class Parameter_Sampler_darcy_flow_deterministic: public HDSA::Parameter_Sampler<RealT>{

private:
  int theta_dim_;

public:
  Parameter_Sampler_darcy_flow_deterministic(const HDSA::Ptr<HDSA::ParameterList> & parlist): HDSA::Parameter_Sampler<RealT>() {
    int L = parlist->sublist("Problem").get("Number of Uncertainty Basis Function", 10);    
    theta_dim_ = (L+1)*(L+1);
  }

  ~Parameter_Sampler_darcy_flow_deterministic() {}

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

}; // Parameter_Sampler_darcy_flow_deterministic

#endif

