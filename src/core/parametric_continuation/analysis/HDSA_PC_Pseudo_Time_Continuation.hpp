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
	bool print_cg_iter_;
	RealT cg_tol_;
	int max_cg_iter_;
    
  protected:

    virtual int Apply_Inverse_Hessian(std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > & P, std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > & W, HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta) const
    {

	  //Build_Hessian(z,theta);
      
      z_out.zeros();
      HDSA::Ptr<HDSA::Vector<RealT> > r = z_in.clone();
      r->set(z_in);
      HDSA::Ptr<HDSA::Vector<RealT> > v = z_in.clone();
      qn_prec_->Apply_Inverse_Hessian_Approximation(*v,*r);
      HDSA::Ptr<HDSA::Vector<RealT> > p = z_in.clone();      
      p->set(*v);
      RealT scalar = r->dot(*p);
      RealT rel_tol = cg_tol_*z_in.norm();
	  RealT r_norm = r->norm();
      int iter = 0;
      HDSA::Ptr<HDSA::Vector<RealT> > w = z_in.clone();
      
      while( (std::sqrt(scalar) > rel_tol) && (r_norm > rel_tol) && (iter < max_cg_iter_) )
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
		r_norm = r->norm();

		if(print_cg_iter_)
		{
			std::cout << "Iteration = " << iter << " with relative residual = " << r_norm/(rel_tol/cg_tol_) << std::endl;
		}
	  }

      if(print_cg_output_)
	  {
		std::cout << "Total iterations = " << iter << std::endl;
		std::cout << "Relative residual = " << r_norm/z_in.norm() << std::endl;
	  }

	  return iter;
    }

	void Build_Hessian(const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta) const
	{
		std::cout << "Beginnning Hessian construction" << std::endl;
		int n = z.dimension();
		std::vector<std::vector<RealT> > H;
		H.resize(n);
		for(int i = 0; i < n; i++)
		{
			H[i].resize(n);
		}

		HDSA::Ptr<HDSA::Vector<RealT> > z_out = z.clone();
		HDSA::ROL_Vector<RealT>& z_out_rol = dynamic_cast<HDSA::ROL_Vector<RealT>&>(*z_out);
		for(int j = 0; j < n; j++)
		{
			HDSA::Ptr<ROL::Vector<RealT> > z_in_rol = z_out_rol.rol_vec->basis(j);
			HDSA::Ptr<HDSA::Vector<RealT> > z_in = HDSA::makePtr<HDSA::ROL_Vector<RealT> >(z_in_rol);
			sen_op_interface_->Apply_Hessian(*z_out,*z_in,z,theta);
			for(int i = 0; i < n; i++)
			{
				HDSA::Ptr<ROL::Vector<RealT> > z_b_rol = z_out_rol.rol_vec->basis(i);
				HDSA::Ptr<HDSA::Vector<RealT> > z_b = HDSA::makePtr<HDSA::ROL_Vector<RealT> >(z_b_rol);
				H[i][j] = z_out->dot(*z_b);
			}
		}

	  std::string name = "H.txt";
	  std::ofstream fout;
      fout.open(name);
	  for(int i = 0; i < n; i++)
	  {
		for(int j = 0 ; j < n; j++)
		{
			fout << H[i][j] << " ";
		}
		fout << std::endl;
	  }
      fout.close();
	  std::cout << "Completed Hessian construction" << std::endl;
	}

  public:
    
    PC_Pseudo_Time_Continuation(const HDSA::Ptr<HDSA::Vector<RealT> > & z_bar, const HDSA::Ptr<HDSA::Vector<RealT> > & theta_bar,
								const HDSA::Ptr<HDSA::PC_Sensitivity_Operator_Interface<RealT> > & sen_op_interface, const HDSA::Ptr<HDSA::PC_Quasi_Newton_Preconditioner<RealT> > & qn_prec,
								const bool use_qn_prec = true, const bool print_cg_output = true, const bool print_cg_iter = false, const RealT cg_tol = 1.e-5, const int max_cg_iter = 100): 
      z_bar_(z_bar), theta_bar_(theta_bar), sen_op_interface_(sen_op_interface), qn_prec_(qn_prec), use_qn_prec_(use_qn_prec), print_cg_output_(print_cg_output), print_cg_iter_(print_cg_iter), cg_tol_(cg_tol), max_cg_iter_(max_cg_iter)
    { }

    virtual ~PC_Pseudo_Time_Continuation()
    { }
    
    void Pseudo_Time_Continuation_Forward_Euler(HDSA::Vector<RealT> & z_star, HDSA::Vector<RealT> & grad_star, const HDSA::Vector<RealT> & theta_star, const int & N) 
    {
      HDSA::Ptr<HDSA::Vector<RealT> > z_current = z_star.clone();
      HDSA::Ptr<HDSA::Vector<RealT> > z_new = z_star.clone();
      HDSA::Ptr<HDSA::Vector<RealT> > z_tmp = z_star.clone();
      HDSA::Ptr<HDSA::Vector<RealT> > grad_current = z_star.clone();
      HDSA::Ptr<HDSA::Vector<RealT> > theta_current = theta_star.clone();
      
	  int num_Hvecs = 0;
	  int num_Bvecs = 0;
	  int num_grads = 0;

      RealT dt = 1.0/static_cast<RealT>(N);
      HDSA::Ptr<HDSA::Vector<RealT> > d_theta = theta_star.clone();
      d_theta->set(theta_star);
      d_theta->axpy(-1.0,*theta_bar_);

      z_current->set(*z_bar_);
      theta_current->set(*theta_bar_);
      sen_op_interface_->Gradient(*grad_current,*z_current,*theta_current);
	  num_grads += 1;
      
      HDSA::Ptr<HDSA::Vector<RealT> > s,y;
      if(use_qn_prec_)
		{
			qn_prec_->Set_N(N);
			s = z_new->clone();
			y = z_new->clone();
		}
	  
      for(int k = 0; k < N; k++)
		{
		std::cout << "The gradient norm after step: " << k << " is " << grad_current->norm() << std::endl;
		std::cout << "Beginning B matvec at time step " << k+1 << std::endl;
		sen_op_interface_->Apply_B(*z_tmp,*d_theta,*z_current,*theta_current);
		num_Bvecs += 1;
		std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > P, W;
		std::cout << "Beginning inverse Hessian matvec at time step " << k+1 << std::endl;
		int iters = Apply_Inverse_Hessian(P,W,*z_new,*z_tmp,*z_current,*theta_current);
		num_Hvecs += iters;
		z_new->scale(-dt);
		z_new->plus(*z_current);

		if(use_qn_prec_ & (k < N-1))
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
		num_grads += 1;
		if(use_qn_prec_ & (k < N-1))
			{
			qn_prec_->Add_Block_Quasi_Newton_Data(P,W);

			y->plus(*grad_current);
			qn_prec_->Add_Parametric_Quasi_Newton_Data(*s,*y);
			}
		}

      z_star.set(*z_current);
      grad_star.set(*grad_current);

	  std::string name = "Forward_Euler_Cost_Report.txt";
	  std::ofstream fout;
      fout.open(name);
	  fout << "Number of B-vector products: " << num_Bvecs << std::endl;
	  fout << "Number of H-vector products: " << num_Hvecs << std::endl;
	  fout << "Number of gradient evaluations: " << num_grads << std::endl;
	  fout << "Number of vectors stored for preconditioner: " << qn_prec_->Get_Number_Vecs_Stored() << std::endl;
      fout.close();

    }

    void Pseudo_Time_Continuation_Modified_Euler(HDSA::Vector<RealT> & z_star, HDSA::Vector<RealT> & grad_star, const HDSA::Vector<RealT> & theta_star, const int & N) 
    {
      HDSA::Ptr<HDSA::Vector<RealT> > z_current = z_star.clone();
      HDSA::Ptr<HDSA::Vector<RealT> > z_new = z_star.clone();
      HDSA::Ptr<HDSA::Vector<RealT> > z_tmp = z_star.clone();
      HDSA::Ptr<HDSA::Vector<RealT> > z_store = z_star.clone();
      HDSA::Ptr<HDSA::Vector<RealT> > grad_current = z_star.clone();
      HDSA::Ptr<HDSA::Vector<RealT> > theta_current = theta_star.clone();

	  int num_Hvecs = 0;
	  int num_Bvecs = 0;
	  int num_grads = 0;
      
      RealT dt = 1.0/static_cast<RealT>(N);
      HDSA::Ptr<HDSA::Vector<RealT> > d_theta = theta_star.clone();
      d_theta->set(theta_star);
      d_theta->axpy(-1.0,*theta_bar_);

      z_current->set(*z_bar_);
      theta_current->set(*theta_bar_);
      sen_op_interface_->Gradient(*grad_current,*z_current,*theta_current);
	  num_grads += 1;

      HDSA::Ptr<HDSA::Vector<RealT> > s,y;
      if(use_qn_prec_)
		{
			qn_prec_->Set_N(2*N);
			s = z_new->clone();
			y = z_new->clone();
		}
      
      for(int k = 0; k < N; k++)
		{
			std::cout << "The gradient norm after step: " << k << " is " << grad_current->norm() << std::endl;
			
			z_store->set(*z_current);
			
			std::cout << "Beginning B matvec at time step " << k+1 << std::endl;
			sen_op_interface_->Apply_B(*z_tmp,*d_theta,*z_current,*theta_current);
			num_Bvecs += 1;
			std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > P, W;
			std::cout << "Beginning inverse Hessian matvec at time step " << k+1 << std::endl;
			int iters = Apply_Inverse_Hessian(P,W,*z_new,*z_tmp,*z_current,*theta_current);
			num_Hvecs += iters;
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
			num_grads += 1;
			if(use_qn_prec_)
			{
				qn_prec_->Add_Block_Quasi_Newton_Data(P,W);

				y->plus(*grad_current);
				qn_prec_->Add_Parametric_Quasi_Newton_Data(*s,*y);
			}
			
			std::cout << "Beginning B matvec at time step " << k+1 << " 1/2" << std::endl;
			sen_op_interface_->Apply_B(*z_tmp,*d_theta,*z_current,*theta_current);
			num_Bvecs += 1;
			P.clear();
			W.clear();
			std::cout << "Beginning inverse Hessian matvec at time step " << k+1 << " 1/2" << std::endl;
			iters = Apply_Inverse_Hessian(P,W,*z_new,*z_tmp,*z_current,*theta_current);
			num_Hvecs += iters;
			z_new->scale(-dt);
			z_new->plus(*z_store);

			if(use_qn_prec_ & (k < N-1))
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
			num_grads += 1;

			if(use_qn_prec_ & (k < N-1))
			{
				qn_prec_->Add_Block_Quasi_Newton_Data(P,W);
				y->plus(*grad_current);
				qn_prec_->Add_Parametric_Quasi_Newton_Data(*s,*y);
			}
		}
      z_star.set(*z_current);
      grad_star.set(*grad_current);

	  std::string name = "Modified_Euler_Cost_Report.txt";
	  std::ofstream fout;
      fout.open(name);
	  fout << "Number of B-vector products: " << num_Bvecs << std::endl;
	  fout << "Number of H-vector products: " << num_Hvecs << std::endl;
	  fout << "Number of gradient evaluations: " << num_grads << std::endl;
	  fout << "Number of vectors stored for preconditioner: " << qn_prec_->Get_Number_Vecs_Stored() << std::endl;
      fout.close();

    }
    
  };

}

#endif
