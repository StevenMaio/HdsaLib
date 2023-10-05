#ifndef PARAMETER_SAMPLER_BAYES_HDSA_MODEL_ERROR_TEST_HPP
#define PARAMETER_SAMPLER_BAYES_HDSA_MODEL_ERROR_TEST_HPP

template <class RealT>
class Parameter_Sampler_bayes_hdsa_model_error_test: public HDSA::Parameter_Sampler<RealT>{

private:
  int theta_dim_;

public:
  Parameter_Sampler_bayes_hdsa_model_error_test(void): HDSA::Parameter_Sampler<RealT>() {
    theta_dim_ = 1;
  }

  ~Parameter_Sampler_bayes_hdsa_model_error_test() {}

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

}; // Parameter_Sampler_bayes_hdsa_model_error_test

#endif

