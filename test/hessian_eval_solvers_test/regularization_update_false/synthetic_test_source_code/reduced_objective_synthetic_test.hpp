#ifndef HDSA_SYNTHETIC_TEST_RS_OBJECTIVE_HPP
#define HDSA_SYNTHETIC_TEST_RS_OBJECTIVE_HPP

template <class RealT>
class synthetic_test_RS_Objective: public HDSA::RS_Objective<RealT> {

private:
  std::vector<RealT> a_;
  int dim_;

public:

  synthetic_test_RS_Objective(std::vector<RealT> & a): a_(a)
  { 
    dim_ = a.size();
  }

  ~synthetic_test_RS_Objective()
  { }

  // evaluate objective function
  RealT value(const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    RealT val = 0.0;
    for(int k = 0; k < dim_; k++)
      {
	val += std::pow( z(k)*a_[k]-theta(k), 2);
      }
    for(int k = dim_; k < z.dimension(); k++)
      {
	val += std::pow( z(k), 2);
      }
    return val;
  }

  // evaluate the gradient with respect to z
  void gradient_z(HDSA::Vector<RealT> & grad, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true) 
  {
    for(int k = 0; k < dim_; k++)
      {
	RealT val = 2.0*(z(k)*a_[k]-theta(k))*a_[k];
	grad.Replace_Element(k,val);
      }
    for(int k = dim_; k < z.dimension(); k++)
      {
	grad.Replace_Element(k,2.0*z(k));
      }
  } 
 
};


#endif
