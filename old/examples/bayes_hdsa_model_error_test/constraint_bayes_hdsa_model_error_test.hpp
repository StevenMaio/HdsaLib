#ifndef HDSA_BAYES_HDSA_MODEL_ERROR_TEST_CONSTRAINT_HPP
#define HDSA_BAYES_HDSA_MODEL_ERROR_TEST_CONSTRAINT_HPP

template <class RealT>
class bayes_hdsa_model_error_test_Constraint: public HDSA::Constraint<RealT> {

private:

public:

  bayes_hdsa_model_error_test_Constraint()
  { }

  ~bayes_hdsa_model_error_test_Constraint()
  { }

  // evaluate constraint residual
  void value(HDSA::Vector<RealT> & r, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    HDSA::Ptr<HDSA::Vector<RealT> > tmp = u.Clone();
    for(int k = 0; k < u.dimension(); k++)
      {
	tmp->Replace_Element(k,std::pow(z(k),3.0));
      }
    r.set(u);
    r.axpy(-1.0,*tmp);
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
    for(int k = 0; k < u.dimension(); k++)
      {
	RealT val = -3.0*std::pow(z(k),2.0)*v(k);
	Jv.Replace_Element(k,val);
      }
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
     for(int k = 0; k < u.dimension(); k++)
      {
	RealT val = -3.0*std::pow(z(k),2.0)*v(k);
	Jv.Replace_Element(k,val);
      }
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
    for(int k = 0; k < z.dimension(); k++)
      {
	RealT val = -6.0*z(k)*lambda(k)*v(k);
	Jv.Replace_Element(k,val);
      }
  }     

  // evaluate the adjoint hessian with respect to z,u vector product, i.e. to compute second derivatives of lambda^T*c(u,z,theta) with respect to z,u
  void hessian_z_u_adjoint(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & lambda, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, 
			   const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true) 
  {
    Jv.zero();
  }  

};


#endif
