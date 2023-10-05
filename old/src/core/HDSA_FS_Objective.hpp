#ifndef HDSA_FS_OBJECTIVE_HPP
#define HDSA_FS_OBJECTIVE_HPP

namespace HDSA
{

template <class RealT>
class FS_Objective {

public:

  FS_Objective() 
  { }

  virtual ~FS_Objective()
  { }

  // evaluate objective function
  virtual RealT value(const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true) = 0;

  // evaluate the gradient with respect to u
  virtual void gradient_u(HDSA::Vector<RealT> & grad, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true) = 0; 

  // evaluate the gradient with respect to z
  virtual void gradient_z(HDSA::Vector<RealT> & grad, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true) = 0; 
  
  // evaluate the u,u hessian vector product
  virtual void hessVec_u_u(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, 
			   const HDSA::Vector<RealT> & theta, const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr)
  {
    RealT vnorm = v.norm();
    if(vnorm == 0.0)
      {
	hv.zero();
      }
    else
      {
	HDSA::Ptr<HDSA::Vector<RealT> > grad;
	if(grad_at_input == HDSA::nullPtr)
	  {
	    grad = hv.Clone();
	    gradient_u(*grad,u,z,theta);
	  }
	else
	  {
	    grad = grad_at_input->Clone();
	    grad->set(*grad_at_input);
	  }
	RealT h = 1.e-4;
	RealT scale = u.norm()/vnorm;
	if(scale == 0.0)
	  {
	    scale = 1.0;
	  }
	HDSA::Ptr<HDSA::Vector<RealT> > vec = u.Clone();
	vec->set(u);
	vec->axpy(scale*h,v);
	gradient_u(hv,*vec,z,theta);
	hv.axpy(-1.0,*grad);
	hv.scale(1.0/(scale*h));
      }
  }

  // evaluate the u,z hessian vector product
  virtual void hessVec_u_z(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z,
			   const HDSA::Vector<RealT> & theta, const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr)
  {
    RealT vnorm = v.norm();
    if(vnorm == 0.0)
      {
	hv.zero();
      }
    else
      {
	HDSA::Ptr<HDSA::Vector<RealT> > grad;
	if(grad_at_input == HDSA::nullPtr)
	  {
	    grad = hv.Clone();
	    gradient_u(*grad,u,z,theta);
	  }
	else
	  {
	    grad = grad_at_input->Clone();
	    grad->set(*grad_at_input);
	  }
	RealT h = 1.e-4;
	RealT scale = z.norm()/vnorm;
	if(scale == 0.0)
	  {
	    scale = 1.0;
	  }
	HDSA::Ptr<HDSA::Vector<RealT> > vec = z.Clone();
	vec->set(z);
	vec->axpy(scale*h,v);
	gradient_u(hv,u,*vec,theta);
	hv.axpy(-1.0,*grad);
	hv.scale(1.0/(scale*h));
      }
  }

  // evaluate the z,z hessian vector product
  virtual void hessVec_z_z(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z,
			   const HDSA::Vector<RealT> & theta, const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr)
  {
    RealT vnorm = v.norm();
    if(vnorm == 0.0)
      {
	hv.zero();
      }
    else
      {
	HDSA::Ptr<HDSA::Vector<RealT> > grad;
	if(grad_at_input == HDSA::nullPtr)
	  {
	    grad = hv.Clone();
	    gradient_z(*grad,u,z,theta);
	  }
	else
	  {
	    grad = grad_at_input->Clone();
	    grad->set(*grad_at_input);
	  }
	RealT h = 1.e-4;
	RealT scale = z.norm()/vnorm;
	if(scale == 0.0)
	  {
	    scale = 1.0;
	  }
	HDSA::Ptr<HDSA::Vector<RealT> > vec = z.Clone();
	vec->set(z);
	vec->axpy(scale*h,v);
	gradient_z(hv,u,*vec,theta);
	hv.axpy(-1.0,*grad);
	hv.scale(1.0/(scale*h));
      }
  }

  // evaluate the z,u hessian vector product
  virtual void hessVec_z_u(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z,
			   const HDSA::Vector<RealT> & theta, const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr)
  {
    RealT vnorm = v.norm();
    if(vnorm == 0.0)
      {
	hv.zero();
      }
    else
      {
	HDSA::Ptr<HDSA::Vector<RealT> > grad;
	if(grad_at_input == HDSA::nullPtr)
	  {
	    grad = hv.Clone();
	    gradient_z(*grad,u,z,theta);
	  }
	else
	  {
	    grad = grad_at_input->Clone();
	    grad->set(*grad_at_input);
	  }
	RealT h = 1.e-4;
	RealT scale = u.norm()/vnorm;
	if(scale == 0.0)
	  {
	    scale = 1.0;
	  }
	HDSA::Ptr<HDSA::Vector<RealT> > vec = u.Clone();
	vec->set(u);
	vec->axpy(scale*h,v);
	gradient_z(hv,*vec,z,theta);
	hv.axpy(-1.0,*grad);
	hv.scale(1.0/(scale*h));
      }
  }

  // evaluate the u,theta hessian vector product
  virtual void hessVec_u_theta(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, 
			       const HDSA::Vector<RealT> & theta, const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr)
  {
    RealT vnorm = v.norm();
    if(vnorm == 0.0)
      {
	hv.zero();
      }
    else
      {
	HDSA::Ptr<HDSA::Vector<RealT> > grad;
	if(grad_at_input == HDSA::nullPtr)
	  {
	    grad = hv.Clone();
	    gradient_u(*grad,u,z,theta);
	  }
	else
	  {
	    grad = grad_at_input->Clone();
	    grad->set(*grad_at_input);
	  }
	RealT h = 1.e-4;
	RealT scale = theta.norm()/vnorm;
	if(scale == 0.0)
	  {
	    scale = 1.0;
	  }
	HDSA::Ptr<HDSA::Vector<RealT> > vec = theta.Clone();
	vec->set(theta);
	vec->axpy(scale*h,v);
	gradient_u(hv,u,z,*vec);
	hv.axpy(-1.0,*grad);
	hv.scale(1.0/(scale*h));
      }
  }

  // evaluate the z,theta hessian vector product
  virtual void hessVec_z_theta(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, 
			       const HDSA::Vector<RealT> & theta, const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr)
  {
    RealT vnorm = v.norm();
    if(vnorm == 0.0)
      {
	hv.zero();
      }
    else
      {
	HDSA::Ptr<HDSA::Vector<RealT> > grad;
	if(grad_at_input == HDSA::nullPtr)
	  {
	    grad = hv.Clone();
	    gradient_z(*grad,u,z,theta);
	  }
	else
	  {
	    grad = grad_at_input->Clone();
	    grad->set(*grad_at_input);
	  }
	RealT h = 1.e-4;
	RealT scale = theta.norm()/vnorm;
	if(scale == 0.0)
	  {
	    scale = 1.0;
	  }
	HDSA::Ptr<HDSA::Vector<RealT> > vec = theta.Clone();
	vec->set(theta);
	vec->axpy(scale*h,v);
	gradient_z(hv,u,z,*vec);
	hv.axpy(-1.0,*grad);
	hv.scale(1.0/(scale*h));
      } 
  }

  // evaluate the theta,u hessian vector product
  virtual void hessVec_theta_u(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, 
			       const HDSA::Vector<RealT> & theta, const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr)
  {
    HDSA::Ptr<HDSA::Vector<RealT> > e = theta.Clone();
    HDSA::Ptr<HDSA::Vector<RealT> > He = u.Clone();
    int theta_dim = theta.dimension();
    for(int i = 0; i < theta_dim; i++)
      {
	e->basis(i);
	He->zero();
	hessVec_u_theta(*He, *e, u, z, theta, true, grad_at_input);
	hv.Replace_Element(i,He->dot(v));
      }
  }

  // evaluate the theta,z hessian vector product
  virtual void hessVec_theta_z(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, 
			       const HDSA::Vector<RealT> & theta, const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr)
  {
    HDSA::Ptr<HDSA::Vector<RealT> > e = theta.Clone();
    HDSA::Ptr<HDSA::Vector<RealT> > He = z.Clone();
    int theta_dim = theta.dimension();
    for(int i = 0; i < theta_dim; i++)
      {
	e->basis(i);
	He->zero();
	hessVec_z_theta(*He, *e, u, z, theta, true, grad_at_input);
	hv.Replace_Element(i,He->dot(v));
      }
  }
 
};

}

#endif
