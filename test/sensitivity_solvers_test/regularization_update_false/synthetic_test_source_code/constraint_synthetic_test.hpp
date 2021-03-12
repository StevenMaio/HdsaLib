#ifndef HDSA_SYNTHETIC_TEST_CONSTRAINT_HPP
#define HDSA_SYNTHETIC_TEST_CONSTRAINT_HPP

template <class RealT>
class synthetic_test_Constraint: public HDSA::Constraint<RealT> {


public:

  synthetic_test_Constraint()
  { }

  ~synthetic_test_Constraint()
  { }

  // evaluate constraint residual
  void value(HDSA::Vector<RealT> & r, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    r.set(u);
    r.axpy(-1.0,z);
  }

  // evaluate the jacobian with respect to u vector product
  void jacobian_u(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, 
		  const HDSA::Vector<RealT> & theta, const bool & update = true) 
  {
    Jv.set(v);
  } 

  // evaluate the inverse jacobian with respect to u vector product
  void jacobian_u_inverse(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, 
			  const HDSA::Vector<RealT> & theta, const bool & update = true) 
  {
    Jv.set(v);
  } 

  // evaluate the jacobian with respect to z vector product
  void jacobian_z(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, 
		  const HDSA::Vector<RealT> & theta, const bool & update = true) 
  {
    Jv.set(v);
    Jv.scale(-1.0);
  } 

  // evaluate the adjoint jacobian with respect to u vector product
  void jacobian_u_adjoint(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, 
			  const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    Jv.set(v);
  } 

  // evaluate the inverse adjoint jacobian with respect to u vector product
  void jacobian_u_adjoint_inverse(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, 
				  const HDSA::Vector<RealT> & theta, const bool & update = true) 
  {
    Jv.set(v);
  }

  // evaluate the adjoint jacobian with respect to z vector product
  void jacobian_z_adjoint(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, 
			  const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    Jv.set(v);
    Jv.scale(-1.0);
  }

  // evaluate the adjoint hessian with respect to u,u vector product, i.e. to compute second derivatives of lambda^T*c(u,z,theta) with respect to u,u
  void hessian_u_u_adjoint(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & lambda, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, 
			   const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true) 
  {
    Jv.zero();
  }

  // evaluate the adjoint hessian with respect to u,z vector product, i.e. to compute second derivatives of lambda^T*c(u,z,theta) with respect to u,z
  void hessian_u_z_adjoint(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & lambda, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, 
			   const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    Jv.zero();
  } 

  // evaluate the adjoint hessian with respect to z,z vector product, i.e. to compute second derivatives of lambda^T*c(u,z,theta) with respect to z,z
  void hessian_z_z_adjoint(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & lambda, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, 
			   const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    Jv.zero();
  }     

  // evaluate the adjoint hessian with respect to z,u vector product, i.e. to compute second derivatives of lambda^T*c(u,z,theta) with respect to z,u
  void hessian_z_u_adjoint(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & lambda, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, 
			   const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true) 
  {
    Jv.zero();
  }  

};


#endif
