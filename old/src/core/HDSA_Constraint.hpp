#ifndef HDSA_CONSTRAINT_HPP
#define HDSA_CONSTRAINT_HPP

namespace HDSA
{

template <class RealT>
class Constraint {

public:

  Constraint() 
  { }

  virtual ~Constraint()
  { }

  // evaluate constraint residual
  virtual void value(HDSA::Vector<RealT> & r, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true) = 0;

  // evaluate the jacobian with respect to u vector product
  virtual void jacobian_u(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, 
			  const HDSA::Vector<RealT> & theta, const bool & update = true) = 0; 

  // evaluate the inverse jacobian with respect to u vector product
  virtual void jacobian_u_inverse(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, 
				  const HDSA::Vector<RealT> & theta, const bool & update = true) = 0; 

  // evaluate the jacobian with respect to z vector product
  virtual void jacobian_z(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, 
			  const HDSA::Vector<RealT> & theta, const bool & update = true) = 0;  

  // evaluate the jacobian with respect to theta vector product
  virtual void jacobian_theta(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, 
			      const HDSA::Vector<RealT> & theta, const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & residual_at_input = HDSA::nullPtr)
  {
    HDSA::Ptr<HDSA::Vector<RealT> > res;
    if(residual_at_input == HDSA::nullPtr)
      {
        res = Jv.Clone();
        value(*res,u,z,theta);
      }
    else
      {
        res = residual_at_input->Clone();
	res->set(*residual_at_input);
      }
    RealT h = 1.e-4;
    RealT scale = theta.norm()/v.norm();
    if(scale == 0.0)
      {
	scale = 1.0;
      }
    HDSA::Ptr<HDSA::Vector<RealT> > vec = theta.Clone();
    vec->set(theta);
    vec->axpy(scale*h,v);
    value(Jv,u,z,*vec);
    Jv.axpy(-1.0,*res);
    Jv.scale(1.0/(scale*h));
  }  

  // evaluate the adjoint jacobian with respect to u vector product
  virtual void jacobian_u_adjoint(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, 
				  const HDSA::Vector<RealT> & theta, const bool & update = true) = 0; 

  // evaluate the inverse adjoint jacobian with respect to u vector product
  virtual void jacobian_u_adjoint_inverse(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, 
					  const HDSA::Vector<RealT> & theta, const bool & update = true) = 0; 

  // evaluate the adjoint jacobian with respect to z vector product
  virtual void jacobian_z_adjoint(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, 
				  const HDSA::Vector<RealT> & theta, const bool & update = true) = 0; 

  // evaluate the adjoint hessian with respect to u,u vector product, i.e. to compute second derivatives of lambda^T*c(u,z,theta) with respect to u,u
  virtual void hessian_u_u_adjoint(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & lambda, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, 
				   const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true) = 0;

  // evaluate the adjoint hessian with respect to u,z vector product, i.e. to compute second derivatives of lambda^T*c(u,z,theta) with respect to u,z
  virtual void hessian_u_z_adjoint(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & lambda, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, 
				   const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true) = 0;   

  // evaluate the adjoint hessian with respect to u,z vector product, i.e. to compute second derivatives of lambda^T*c(u,z,theta) with respect to u,theta
  virtual void hessian_u_theta_adjoint(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & lambda, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, 
				       const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    HDSA::Ptr<HDSA::Vector<RealT> > Ju = Jv.Clone();
    jacobian_u_adjoint(*Ju,lambda,u,z,theta);
    RealT h = 1.e-4;
    RealT scale = theta.norm()/v.norm();
    if(scale == 0.0)
      {
	scale = 1.0;
      }
    HDSA::Ptr<HDSA::Vector<RealT> > vec = theta.Clone();
    vec->set(theta);
    vec->axpy(scale*h,v);
    jacobian_u_adjoint(Jv,lambda,u,z,*vec);
    Jv.axpy(-1.0,*Ju);
    Jv.scale(1.0/(scale*h));
  }   

  // evaluate the adjoint hessian with respect to z,z vector product, i.e. to compute second derivatives of lambda^T*c(u,z,theta) with respect to z,z
  virtual void hessian_z_z_adjoint(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & lambda, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, 
				   const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true) = 0;     

  // evaluate the adjoint hessian with respect to z,u vector product, i.e. to compute second derivatives of lambda^T*c(u,z,theta) with respect to z,u
  virtual void hessian_z_u_adjoint(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & lambda, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, 
				   const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true) = 0;    

  // evaluate the adjoint hessian with respect to u,z vector product, i.e. to compute second derivatives of lambda^T*c(u,z,theta) with respect to z,theta
  virtual void hessian_z_theta_adjoint(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & lambda, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, 
				       const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    HDSA::Ptr<HDSA::Vector<RealT> > Jz = Jv.Clone();
    jacobian_z_adjoint(*Jz,lambda,u,z,theta);
    RealT h = 1.e-4;
    RealT scale = theta.norm()/v.norm();
    if(scale == 0.0)
      {
	scale = 1.0;
      }
    HDSA::Ptr<HDSA::Vector<RealT> > vec = theta.Clone();
    vec->set(theta);
    vec->axpy(scale*h,v);
    jacobian_z_adjoint(Jv,lambda,u,z,*vec);
    Jv.axpy(-1.0,*Jz);
    Jv.scale(1.0/(scale*h));
  }   
 
};

}

#endif
