#ifndef HDSA_PC_PSEUDO_TIME_CONTINUATION_HPP
#define HDSA_PC_PSEUDO_TIME_CONTINUATION_HPP

namespace HDSA
{

  template <class RealT>
  class PC_Pseudo_Time_Continuation{

  private:
    HDSA::Ptr<HDSA::Vector<RealT> > z_bar_;
    HDSA::Ptr<HDSA::Vector<RealT> > theta_bar_;
    HDSA::Ptr<HDSA::PC_Sensitivity_Operator_Interface<RealT> > sen_op_;
    bool use_bfgs_prec_;
    int num_bfgs_vecs_;
    std::vector<RealT> rho_;
    HDSA::Ptr<HDSA::MultiVector<RealT> > s_;
    HDSA::Ptr<HDSA::MultiVector<RealT> > y_;
    bool print_cg_output_;
    
  protected:

    // Overload this function if a better initialization is available
    virtual void Apply_Initial_Inverse_BFGS_Hessian(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const
    {
      z_out.set(z_in);
    }

    virtual void Apply_Inverse_Hessian(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta) const
    {
      RealT tol = 1.e-5;
      int max_iter = z_out.dimension();
      
      z_out.zeros();
      HDSA::Ptr<HDSA::Vector<RealT> > r = z_in.clone();
      r->set(z_in);
      HDSA::Ptr<HDSA::Vector<RealT> > v = z_in.clone();
      Apply_Inverse_BFGS_Hessian(*v,*r,num_bfgs_vecs_);
      HDSA::Ptr<HDSA::Vector<RealT> > p = z_in.clone();      
      p->set(*v);
      RealT scalar = r->dot(*p);
      RealT rel_tol = tol*z_in.norm();
      int iter = 0;
      HDSA::Ptr<HDSA::Vector<RealT> > w = z_in.clone();

      while( (std::sqrt(scalar) > rel_tol) && (r->norm() > rel_tol) && (iter < max_iter) )
	{
	  iter += 1;
	  sen_op_->Apply_Hessian(*w,*p,z,theta);
	  RealT alpha = scalar/(w->dot(*p));
	  z_out.axpy(alpha,*p);
	  r->axpy(-alpha,*w);
	  Apply_Inverse_BFGS_Hessian(*v,*r,num_bfgs_vecs_);
	  RealT scalar_old = scalar;
	  scalar = v->dot(*r);
	  p->scale(scalar/scalar_old);
	  p->plus(*v);
	}
      if(print_cg_output_)
	{
	  std::cout << "Total iterations = " << iter << std::endl;
	  std::cout << "Relative residual = " << std::sqrt(scalar)/z_in.norm() << std::endl;
	}
    }

  public:
    
    PC_Pseudo_Time_Continuation(const HDSA::Ptr<HDSA::Vector<RealT> > & z_bar, const HDSA::Ptr<HDSA::Vector<RealT> > & theta_bar, const HDSA::Ptr<HDSA::PC_Sensitivity_Operator_Interface<RealT> > & sen_op): 
      z_bar_(z_bar), theta_bar_(theta_bar), sen_op_(sen_op)
    {
      use_bfgs_prec_ = true;
      num_bfgs_vecs_ = 0;
      print_cg_output_ = false;
    }

    virtual ~PC_Pseudo_Time_Continuation()
    { }
    
    void Pseudo_Time_Continuation_Forward_Euler(HDSA::Vector<RealT> & z_star, HDSA::Vector<RealT> & grad_star, const HDSA::Vector<RealT> & theta_star, const int & N) 
    {
      HDSA::Ptr<HDSA::Vector<RealT> > z_current = z_star.clone();
      HDSA::Ptr<HDSA::Vector<RealT> > z_new = z_star.clone();
      HDSA::Ptr<HDSA::Vector<RealT> > z_tmp = z_star.clone();
      HDSA::Ptr<HDSA::Vector<RealT> > grad_current = z_star.clone();
      HDSA::Ptr<HDSA::Vector<RealT> > theta_current = theta_star.clone();

      if(use_bfgs_prec_)
	{
	  rho_ = std::vector<RealT>(N);
	  s_ = HDSA::makePtr<HDSA::MultiVector<RealT> >(N,z_star);
	  y_ = HDSA::makePtr<HDSA::MultiVector<RealT> >(N,z_star);
	}
      
      RealT dt = 1.0/static_cast<RealT>(N);
      HDSA::Ptr<HDSA::Vector<RealT> > d_theta = theta_star.clone();
      d_theta->set(theta_star);
      d_theta->axpy(-1.0,*theta_bar_);

      z_current->set(*z_bar_);
      theta_current->set(*theta_bar_);
      sen_op_->Gradient(*grad_current,*z_current,*theta_current);

      for(int k = 0; k < N; k++)
	{
	  sen_op_->Apply_B(*z_tmp,*d_theta,*z_current,*theta_current);
	  Apply_Inverse_Hessian(*z_new,*z_tmp,*z_current,*theta_current);
	  z_new->scale(-dt);
	  z_new->plus(*z_current);

	  if(use_bfgs_prec_)
	    {
	      (*s_)[num_bfgs_vecs_]->set(*z_new);
	      (*s_)[num_bfgs_vecs_]->axpy(-1.0,*z_current);

	      (*y_)[num_bfgs_vecs_]->set(*grad_current);
	      (*y_)[num_bfgs_vecs_]->axpy(dt,*z_tmp);
	      (*y_)[num_bfgs_vecs_]->scale(-1.0);
	    }
	  
	  z_current->set(*z_new);
	  theta_current->axpy(dt,*d_theta);
	  sen_op_->Gradient(*grad_current,*z_current,*theta_current);
	  if(use_bfgs_prec_)
	    {
	      (*y_)[num_bfgs_vecs_]->plus(*grad_current);
	      rho_[num_bfgs_vecs_] = 1.0/((*y_)[num_bfgs_vecs_]->dot(*(*s_)[num_bfgs_vecs_]));
	      num_bfgs_vecs_ += 1;
	    }
	}
      z_star.set(*z_current);
      grad_star.set(*grad_current);
    }

    void Pseudo_Time_Continuation_Modified_Euler(HDSA::Vector<RealT> & z_star, HDSA::Vector<RealT> & grad_star, const HDSA::Vector<RealT> & theta_star, const int & N) 
    {
      HDSA::Ptr<HDSA::Vector<RealT> > z_current = z_star.clone();
      HDSA::Ptr<HDSA::Vector<RealT> > z_new = z_star.clone();
      HDSA::Ptr<HDSA::Vector<RealT> > z_tmp = z_star.clone();
      HDSA::Ptr<HDSA::Vector<RealT> > z_store = z_star.clone();
      HDSA::Ptr<HDSA::Vector<RealT> > grad_current = z_star.clone();
      HDSA::Ptr<HDSA::Vector<RealT> > theta_current = theta_star.clone();

      if(use_bfgs_prec_)
	{
	  rho_ = std::vector<RealT>(N);
	  s_ = HDSA::makePtr<HDSA::MultiVector<RealT> >(N,z_star);
	  y_ = HDSA::makePtr<HDSA::MultiVector<RealT> >(N,z_star);
	}
      
      RealT dt = 1.0/static_cast<RealT>(N);
      HDSA::Ptr<HDSA::Vector<RealT> > d_theta = theta_star.clone();
      d_theta->set(theta_star);
      d_theta->axpy(-1.0,*theta_bar_);

      z_current->set(*z_bar_);
      theta_current->set(*theta_bar_);
      sen_op_->Gradient(*grad_current,*z_current,*theta_current);

      for(int k = 0; k < N; k++)
	{
	  sen_op_->Apply_B(*z_tmp,*d_theta,*z_current,*theta_current);
	  if(use_bfgs_prec_)
	    {
	      (*y_)[num_bfgs_vecs_]->set(*grad_current);
	      (*y_)[num_bfgs_vecs_]->axpy(dt,*z_tmp);
	      (*y_)[num_bfgs_vecs_]->scale(-1.0);
	    }
	  Apply_Inverse_Hessian(*z_new,*z_tmp,*z_current,*theta_current);
	  z_new->scale(-0.5*dt);
	  z_new->plus(*z_current);

	  z_store->set(*z_new);
	  z_new->zeros();
	  theta_current->axpy(0.5*dt,*d_theta);
	  sen_op_->Apply_B(*z_tmp,*d_theta,*z_store,*theta_current);
	  Apply_Inverse_Hessian(*z_new,*z_tmp,*z_store,*theta_current);
	  z_new->scale(-dt);
	  z_new->plus(*z_current);

	  if(use_bfgs_prec_)
	    {
	      (*s_)[num_bfgs_vecs_]->set(*z_new);
	      (*s_)[num_bfgs_vecs_]->axpy(-1.0,*z_current);
	    }
	  
	  z_current->set(*z_new);
	  theta_current->axpy(0.5*dt,*d_theta);
	  sen_op_->Gradient(*grad_current,*z_current,*theta_current);
	  if(use_bfgs_prec_)
	    {
	      (*y_)[num_bfgs_vecs_]->plus(*grad_current);
	      rho_[num_bfgs_vecs_] = 1.0/((*y_)[num_bfgs_vecs_]->dot(*(*s_)[num_bfgs_vecs_]));
	      num_bfgs_vecs_ += 1;
	    }
	}
      z_star.set(*z_current);
      grad_star.set(*grad_current);
    }

    void Apply_Inverse_BFGS_Hessian(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in, const int & num_vecs) const
    {
      if(!use_bfgs_prec_ || num_vecs==0)
	{
	  Apply_Initial_Inverse_BFGS_Hessian(z_out,z_in);
	}
      else
	{
	  RealT alpha = (*s_)[num_vecs-1]->dot(z_in);
	  HDSA::Ptr<HDSA::Vector<RealT> > tmp = z_out.clone();
	  tmp->set(z_in);
	  tmp->axpy(-rho_[num_vecs-1]*alpha,*(*y_)[num_vecs-1]);
	  HDSA::Ptr<HDSA::Vector<RealT> > tmp_out = z_out.clone();
	  if(num_vecs==1)
	    {
	      Apply_Initial_Inverse_BFGS_Hessian(*tmp_out,*tmp);
	    }
	  else
	    {
	      Apply_Inverse_BFGS_Hessian(*tmp_out,*tmp,num_vecs-1);
	    }
	  z_out.set(*tmp_out);
	  RealT coeff = rho_[num_vecs-1] * (alpha - (*y_)[num_vecs-1]->dot(*tmp_out));
	  z_out.axpy(coeff,*(*s_)[num_vecs-1]);
	}

    }
    
  };

}

#endif
