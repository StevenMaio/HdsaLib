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
	samples[0]->Replace_Element(0,1.0); samples[0]->Replace_Element(1,2.0); samples[0]->Replace_Element(2,2.0); samples[0]->Replace_Element(3,1.0);
      }
    return samples;
  }

}; // Parameter_Sampler_synthetic_test
