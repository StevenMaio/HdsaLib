#ifndef HDSA_SYNTHETIC_TEST_RS_OBJECTIVE_HPP
#define HDSA_SYNTHETIC_TEST_RS_OBJECTIVE_HPP

template <class RealT>
class synthetic_test_RS_Objective: public HDSA::RS_Objective<RealT> {

private:
  std::vector<RealT> a_;
  std::vector<RealT> b_;
  int dim_;

public:

  synthetic_test_RS_Objective(std::vector<RealT> & a, std::vector<RealT> & b): a_(a), b_(b)
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
	val += std::pow( z(k)-a_[k]*theta(k), 2) + b_[k]*std::pow(z(k),2);
      }
    for(int k = dim_; k < z.dimension(); k++)
      {
	val += std::pow( z(k), 2);
      }
    val = 0.5*val;
    return val;
  }

  // evaluate the gradient with respect to z
  void gradient_z(HDSA::Vector<RealT> & grad, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true) 
  {
    for(int k = 0; k < dim_; k++)
      {
	RealT val = (z(k)-a_[k]*theta(k)) + b_[k]*z(k);
	grad.Replace_Element(k,val);
      }
    for(int k = dim_; k < z.dimension(); k++)
      {
	grad.Replace_Element(k,z(k));
      }
  } 

  // evaluate the misfit z,z hessian vector product (for computed likelihood informed subspaces only)
  void Misfit_hessVec_z_z(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, 
			  const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr)
  { 
    hv.set(v);
  }

  // evaluate the regularization z,z hessian vector product (for computed likelihood informed subspaces only)
  void Regularization_hessVec_z_z(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, 
				  const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr)
  {
    for(int k = 0; k < dim_; k++)
      {
	hv.Replace_Element(k,b_[k]*v(k));
      }
    for(int k = dim_; k < v.dimension(); k++)
      {
        hv.Replace_Element(k,v(k));
      }
  } 

};


#endif
