#ifndef HDSA_ROL_B_TRANSPOSE_HPP
#define HDSA_ROL_B_TRANSPOSE_HPP

template <class RealT>
class ROL_B_Transpose {

public:

  ROL_B_Transpose()
  { }

  // evaluate the theta,z hessian vector product, i.e. -B^T in HDSA
  virtual void hessVec_theta_z(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, 
			       const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr) = 0;
 
};


#endif
