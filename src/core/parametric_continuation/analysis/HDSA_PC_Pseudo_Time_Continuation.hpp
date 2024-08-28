#ifndef HDSA_PC_PSEUDO_TIME_CONTINUATION_HPP
#define HDSA_PC_PSEUDO_TIME_CONTINUATION_HPP

namespace HDSA
{

  template <class RealT>
  class PC_Pseudo_Time_Continuation{

  private:
    HDSA::Ptr<HDSA::Vector<RealT> > z_bar_;
    HDSA::Ptr<HDSA::Vector<RealT> > theta_bar_;
    HDSA::Ptr<HDSA::PC_Sensitivity_Operator_Interface<RealT> > sen_op_interface_;
    HDSA::Ptr<HDSA::PC_Quasi_Newton_Preconditioner<RealT> > qn_prec_;
    bool use_qn_prec_;
    bool print_cg_output_;
    
  protected:

    virtual void Apply_Inverse_Hessian(std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > & P, std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > & W, HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta) const
    {
      RealT tol = 1.e-5;
      int max_iter = z_out.dimension();
      
      z_out.zeros();
      HDSA::Ptr<HDSA::Vector<RealT> > r = z_in.clone();
      r->set(z_in);
      HDSA::Ptr<HDSA::Vector<RealT> > v = z_in.clone();
      qn_prec_->Apply_Inverse_Hessian_Approximation(*v,*r);
      HDSA::Ptr<HDSA::Vector<RealT> > p = z_in.clone();      
      p->set(*v);
      RealT scalar = r->dot(*p);
      RealT rel_tol = tol*z_in.norm();
      int iter = 0;
      HDSA::Ptr<HDSA::Vector<RealT> > w = z_in.clone();
      
      while( (std::sqrt(scalar) > rel_tol) && (r->norm() > rel_tol) && (iter < max_iter) )
	{
	  iter += 1;
	  sen_op_interface_->Apply_Hessian(*w,*p,z,theta);

	  HDSA::Ptr<HDSA::Vector<RealT> > p_tmp = z_in.clone();  
	  p_tmp->set(*p);
	  HDSA::Ptr<HDSA::Vector<RealT> > w_tmp = z_in.clone();  
	  w_tmp->set(*w);
	  P.push_back(p_tmp);
	  W.push_back(w_tmp);

	  RealT alpha = scalar/(w->dot(*p));
	  z_out.axpy(alpha,*p);
	  r->axpy(-alpha,*w);
	  qn_prec_->Apply_Inverse_Hessian_Approximation(*v,*r);
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
    
    PC_Pseudo_Time_Continuation(const HDSA::Ptr<HDSA::Vector<RealT> > & z_bar, const HDSA::Ptr<HDSA::Vector<RealT> > & theta_bar,
				const HDSA::Ptr<HDSA::PC_Sensitivity_Operator_Interface<RealT> > & sen_op_interface, const HDSA::Ptr<HDSA::PC_Quasi_Newton_Preconditioner<RealT> > & qn_prec): 
      z_bar_(z_bar), theta_bar_(theta_bar), sen_op_interface_(sen_op_interface), qn_prec_(qn_prec)
    {
      use_qn_prec_ = true;
      print_cg_output_ = true;
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
      
      RealT dt = 1.0/static_cast<RealT>(N);
      HDSA::Ptr<HDSA::Vector<RealT> > d_theta = theta_star.clone();
      d_theta->set(theta_star);
      d_theta->axpy(-1.0,*theta_bar_);

      z_current->set(*z_bar_);
      theta_current->set(*theta_bar_);
      sen_op_interface_->Gradient(*grad_current,*z_current,*theta_current);

      HDSA::Ptr<HDSA::Vector<RealT> > s,y;
      if(use_qn_prec_)
	{
	  qn_prec_->Set_N(N);
	  s = z_new->clone();
	  y = z_new->clone();
	}
	  
      for(int k = 0; k < N; k++)
	{
	  sen_op_interface_->Apply_B(*z_tmp,*d_theta,*z_current,*theta_current);
	  std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > P, W;
	  Apply_Inverse_Hessian(P,W,*z_new,*z_tmp,*z_current,*theta_current);
	  z_new->scale(-dt);
	  z_new->plus(*z_current);

	  if(use_qn_prec_)
	    {
	      s->set(*z_new);
	      s->axpy(-1.0,*z_current);

	      y->set(*grad_current);
	      y->axpy(dt,*z_tmp);
	      y->scale(-1.0);
	    }
	  
	  z_current->set(*z_new);
	  theta_current->axpy(dt,*d_theta);
	  sen_op_interface_->Gradient(*grad_current,*z_current,*theta_current);
	  if(use_qn_prec_)
	    {
	      qn_prec_->Add_Block_Quasi_Newton_Data(P,W);

	      y->plus(*grad_current);
	      qn_prec_->Add_Parametric_Quasi_Newton_Data(*s,*y);
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
      
      RealT dt = 1.0/static_cast<RealT>(N);
      HDSA::Ptr<HDSA::Vector<RealT> > d_theta = theta_star.clone();
      d_theta->set(theta_star);
      d_theta->axpy(-1.0,*theta_bar_);

      z_current->set(*z_bar_);
      theta_current->set(*theta_bar_);
      sen_op_interface_->Gradient(*grad_current,*z_current,*theta_current);

      HDSA::Ptr<HDSA::Vector<RealT> > s,y;
      if(use_qn_prec_)
	{
	  qn_prec_->Set_N(2*N);
	  s = z_new->clone();
	  y = z_new->clone();
	}
      
      for(int k = 0; k < N; k++)
	{
	  z_store->set(*z_current);
	  
	  sen_op_interface_->Apply_B(*z_tmp,*d_theta,*z_current,*theta_current);
	  std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > P, W;
	  Apply_Inverse_Hessian(P,W,*z_new,*z_tmp,*z_current,*theta_current);
	  z_new->scale(-0.5*dt);
	  z_new->plus(*z_current);

	  if(use_qn_prec_)
	    {
	      s->set(*z_new);
	      s->axpy(-1.0,*z_current);
	      
	      y->set(*grad_current);
	      y->axpy(0.5*dt,*z_tmp);
	      y->scale(-1.0);
	    }


	  z_current->set(*z_new);
	  z_new->zeros();
	  theta_current->axpy(0.5*dt,*d_theta);
	  
	  sen_op_interface_->Gradient(*grad_current,*z_current,*theta_current);
	  if(use_qn_prec_)
	    {
	      qn_prec_->Add_Block_Quasi_Newton_Data(P,W);

	      y->plus(*grad_current);
	      qn_prec_->Add_Parametric_Quasi_Newton_Data(*s,*y);
	    }
	  
	  sen_op_interface_->Apply_B(*z_tmp,*d_theta,*z_current,*theta_current);
	  P.clear();
	  W.clear();
	  Apply_Inverse_Hessian(P,W,*z_new,*z_tmp,*z_current,*theta_current);
	  z_new->scale(-dt);
	  z_new->plus(*z_store);

	  if(use_qn_prec_)
	    {
	      s->set(*z_new);
	      s->axpy(-1.0,*z_current);

	      y->set(*grad_current);
	      y->axpy(0.5*dt,*z_tmp);
	      y->scale(-1.0);
	    }
	  
	  z_current->set(*z_new);
	  theta_current->axpy(0.5*dt,*d_theta);
	  sen_op_interface_->Gradient(*grad_current,*z_current,*theta_current);

	  if(use_qn_prec_)
	    {
	      qn_prec_->Add_Block_Quasi_Newton_Data(P,W);

	      y->plus(*grad_current);
	      qn_prec_->Add_Parametric_Quasi_Newton_Data(*s,*y);
	    }
	}
      z_star.set(*z_current);
      grad_star.set(*grad_current);
    }
    
  };

}

#endif
