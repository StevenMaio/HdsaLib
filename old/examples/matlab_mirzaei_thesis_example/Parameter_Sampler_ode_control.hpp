template <class RealT>
class Parameter_Sampler_ode_control: public HDSA::Parameter_Sampler<RealT>{

private:
  int theta_dim_;

public:
  Parameter_Sampler_ode_control(const int & theta_dim): HDSA::Parameter_Sampler<RealT>() {
    theta_dim_ = theta_dim;
  }

  ~Parameter_Sampler_ode_control() {}

  std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > Draw_Samples(int num_samp)
  {
    std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > samples;
    samples.resize(num_samp);
    for(int k = 0; k < num_samp; k++)
      {
	samples[k] = HDSA::makePtr<Std_Vector<RealT> >(theta_dim_);
	samples[k]->zero();
      }
    return samples;
  }

}; // Parameter_Sampler_synthetic_test
