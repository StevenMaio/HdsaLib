#ifndef PARAMETER_SAMPLER_BRANDES_GRIESSE_HPP
#define PARAMETER_SAMPLER_BRANDES_GRIESSE_HPP

template <class RealT>
class Parameter_Sampler_brandes_griesse: public HDSA::Parameter_Sampler<RealT>{

private:
  int theta_dim_;

public:
  Parameter_Sampler_brandes_griesse(void): HDSA::Parameter_Sampler<RealT>() {
    theta_dim_ = 20;
  }

  ~Parameter_Sampler_brandes_griesse() {}

  std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > Draw_Samples(int num_samp)
  {
    std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > samples;
    samples.resize(num_samp);
    for(int k = 0; k < num_samp; k++)
      {
	samples[k] = HDSA::makePtr<Std_Vector<RealT> >(theta_dim_);
	samples[k]->Replace_Element(0,2.0);
	samples[k]->Replace_Element(1,1.0);
	for(int i = 2; i < theta_dim_; i++)
	  {
	    samples[k]->Replace_Element(i,0.0);
	  }

      }
    return samples;
  }

}; // Parameter_Sampler_brandes_griesse

#endif

