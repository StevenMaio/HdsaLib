#ifndef HDSA_SYNTHETIC_TEST_FS_OBJECTIVE_HPP
#define HDSA_SYNTHETIC_TEST_FS_OBJECTIVE_HPP

template <class RealT>
class synthetic_test_FS_Objective: public HDSA::FS_Objective<RealT> {

private:
  std::vector<RealT> a_;

public:

  synthetic_test_FS_Objective(std::vector<RealT> & a): a_(a)
  { }

  ~synthetic_test_FS_Objective()
  { }

  // evaluate objective function
  RealT value(const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    RealT val = 0.0;
    int dim = a_.size();
    for(int k = 0; k < dim; k++)
      {
	val += std::pow( u(k)-a_[k]*theta(k), 2);
      }
    for(int k = dim; k < u.dimension(); k++)
      {
	val += std::pow( u(k), 2);
      }
    return val;    
  }

  // evaluate the gradient with respect to u
  void gradient_u(HDSA::Vector<RealT> & grad, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    int dim = a_.size();
    for(int k = 0; k < dim; k++)
      {
	RealT val = 2.0*(u(k)-a_[k]*theta(k));
	grad.Replace_Element(k,val);
      }
    for(int k = dim; k < u.dimension(); k++)
      {
	grad.Replace_Element(k,2.0*u(k));
      }
  }

  // evaluate the gradient with respect to z
  void gradient_z(HDSA::Vector<RealT> & grad, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    grad.zero();
  }
  
};


#endif
