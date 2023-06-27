#ifndef HDSA_MD_INTERFACE_POISSON_1D_TEST_HPP
#define HDSA_MD_INTERFACE_POISSON_1D_TEST_HPP

template <class RealT>
class Model_Discrepancy_Interface_Poisson_1D_Test : public HDSA::Model_Discrepancy_Interface_Elliptic_Prior<RealT> {

private:
  int m_; // Mesh resolution
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > x_; // Mesh nodes on [0,1]
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > target_; // Target state solution in misfit
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > S_; // Stiffness matrix
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > M_; // Mass matrix
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > robin_bc_; // Robin boundary condition matrix
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > sol_op_lofi_; // Low-fidelity solution operator
  RealT diff_coeff_; // Diffusion coefficient
  RealT robin_coeff_; // Robin boundary condition coefficient
  RealT reg_coeff_; // Regularization coefficient
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > E_L_; // State elliptic operator
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > E_Gamma_; // Control elliptic

public:

  Model_Discrepancy_Interface_Poisson_1D_Test()
  {  
    m_ = 200;
    RealT h = 1.0/static_cast<RealT>(m_-1);
    x_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    target_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    for(int k = 0; k < m_; k++)
      {
	x_->Replace_Element(k,0,static_cast<RealT>(k)/static_cast<RealT>(m_-1));
	target_->Replace_Element(k,0,50.0-30.0*std::pow((*x_)(k,0)-0.5,2.0));
      }

    S_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);
    M_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);
    robin_bc_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);

    S_->Replace_Element(0,0,1.0/h);
    S_->Replace_Element(0,1,-1.0/h);
    for(int i = 1; i < m_-1; i++)
      {
	S_->Replace_Element(i,i,2.0/h);
	S_->Replace_Element(i,i-1,-1.0/h);
	S_->Replace_Element(i,i+1,-1.0/h);
      }
    S_->Replace_Element(m_-1,m_-2,-1.0/h);
    S_->Replace_Element(m_-1,m_-1,1.0/h);    

    M_->Replace_Element(0,0,(1.0/3.0)*h);
    M_->Replace_Element(0,1,(1.0/6.0)*h);
    for(int i = 1; i < m_-1; i++)
      {
	M_->Replace_Element(i,i,(2.0/3.0)*h);
	M_->Replace_Element(i,i-1,(1.0/6.0)*h);
	M_->Replace_Element(i,i+1,(1.0/6.0)*h);
      }
    M_->Replace_Element(m_-1,m_-2,(1.0/6.0)*h);
    M_->Replace_Element(m_-1,m_-1,(1.0/3.0)*h);

    robin_bc_->Replace_Element(0,0,1.0);
    robin_bc_->Replace_Element(m_-1,m_-1,1.0);

    diff_coeff_ = 1.0;
    robin_coeff_ = 2.0;
    reg_coeff_ = 10.0;

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > tmp_lofi = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);
    for(int i = 0; i < m_; i++)
      {
	for(int j = 0; j < m_; j++)
	  {
	    RealT val = diff_coeff_*(*S_)(i,j) + robin_coeff_*(*robin_bc_)(i,j);
	    tmp_lofi->Replace_Element(i,j,(1.e-2)*val);
	  }
      }
    sol_op_lofi_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);
    HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*tmp_lofi,*sol_op_lofi_,*M_);

    E_L_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);
    E_Gamma_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);
    for(int i = 0; i < m_; i++)
      {
	for(int j = 0; j < m_; j++)
	  {
	    RealT val = 10.0*( (5.e-2)*(*S_)(i,j) + (*M_)(i,j) );
	    E_L_->Replace_Element(i,j,val);
	    val = (1.e2)*( (1.e-3)*(*S_)(i,j) + (*M_)(i,j) );
	    E_Gamma_->Replace_Element(i,j,val);
	  }
      }

    int num_sing_vals = 200;
    int oversampling = 0;
    int num_subspace_iters = 1;
    HDSA::Ptr<HDSA::Vector<RealT> > u_vec = HDSA::makePtr<Std_Vector<RealT> >(m_);
    HDSA::Model_Discrepancy_Interface_Elliptic_Prior<RealT>::Compute_Elliptic_GSVD(num_sing_vals,oversampling,num_subspace_iters,*u_vec);
  }

  virtual ~Model_Discrepancy_Interface_Poisson_1D_Test()
  { }

  void Apply_u_Elliptic_Operator_Inverse(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const 
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > b = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    const Std_Vector<RealT>& u_in_std = dynamic_cast<const Std_Vector<RealT>&>(u_in);
    Std_Vector<RealT>& u_out_std = dynamic_cast<Std_Vector<RealT>&>(u_out);
    for(int k = 0; k < m_; k++)
      {
	b->Replace_Element(k,0,u_in_std(k));
      }

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > x = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*E_L_,*x,*b);

    for(int k = 0; k < m_; k++)
      {
	u_out_std.Replace_Element(k,(*x)(k,0));
      }
  }

  void Apply_u_Elliptic_Operator_Inverse_Transpose(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const 
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > b = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    const Std_Vector<RealT>& u_in_std = dynamic_cast<const Std_Vector<RealT>&>(u_in);
    Std_Vector<RealT>& u_out_std = dynamic_cast<Std_Vector<RealT>&>(u_out);
    for(int k = 0; k < m_; k++)
      {
	b->Replace_Element(k,0,u_in_std(k));
      }

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > x = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*E_L_,*x,*b);

    for(int k = 0; k < m_; k++)
      {
	u_out_std.Replace_Element(k,(*x)(k,0));
      }
  }

  void Apply_u_Mass_Mat(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const 
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > b = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    const Std_Vector<RealT>& u_in_std = dynamic_cast<const Std_Vector<RealT>&>(u_in);
    Std_Vector<RealT>& u_out_std = dynamic_cast<Std_Vector<RealT>&>(u_out);
    for(int k = 0; k < m_; k++)
      {
	b->Replace_Element(k,0,u_in_std(k));
      }

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > x = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    M_->Multiply(*x,*b);

    for(int k = 0; k < m_; k++)
      {
	u_out_std.Replace_Element(k,(*x)(k,0));
      }
  }

  void Apply_u_Mass_Mat_Inverse(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const 
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > b = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    const Std_Vector<RealT>& u_in_std = dynamic_cast<const Std_Vector<RealT>&>(u_in);
    Std_Vector<RealT>& u_out_std = dynamic_cast<Std_Vector<RealT>&>(u_out);
    for(int k = 0; k < m_; k++)
      {
	b->Replace_Element(k,0,u_in_std(k));
      }

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > x = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*M_,*x,*b);

    for(int k = 0; k < m_; k++)
      {
	u_out_std.Replace_Element(k,(*x)(k,0));
      }
  }

  void Apply_Gamma_Mat_Inverse(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const 
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > b = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    const Std_Vector<RealT>& z_in_std = dynamic_cast<const Std_Vector<RealT>&>(z_in);
    Std_Vector<RealT>& z_out_std = dynamic_cast<Std_Vector<RealT>&>(z_out);
    for(int k = 0; k < m_; k++)
      {
	b->Replace_Element(k,0,z_in_std(k));
      }

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > tmp1 = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*E_Gamma_,*tmp1,*b);
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > tmp2 = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    M_->Multiply(*tmp2,*tmp1);
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > tmp3 = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*E_Gamma_,*tmp3,*tmp2);

    for(int k = 0; k < m_; k++)
      {
	z_out_std.Replace_Element(k,(*tmp3)(k,0));
      }
  }

  void Apply_Solution_Operator_z_Jacobian_Transpose(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & u_in, const HDSA::Vector<RealT> & z) const 
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > b = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    const Std_Vector<RealT>& u_in_std = dynamic_cast<const Std_Vector<RealT>&>(u_in);
    for(int k = 0; k < m_; k++)
      {
	b->Replace_Element(k,0,u_in_std(k));
      }

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > tmp = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    sol_op_lofi_->Multiply(*tmp,*b,true);

    Std_Vector<RealT> z_out_std = dynamic_cast<const Std_Vector<RealT>&>(z_out);
    for(int k = 0; k < m_; k++)
      {
	z_out_std.Replace_Element(k,(*tmp)(k,0));
      }
  }

  void Apply_RS_Hessian(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in, const HDSA::Vector<RealT> & z) const 
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > b = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    const Std_Vector<RealT>& z_in_std = dynamic_cast<const Std_Vector<RealT>&>(z_in);
    Std_Vector<RealT>& z_out_std = dynamic_cast<Std_Vector<RealT>&>(z_out);
    for(int k = 0; k < m_; k++)
      {
	b->Replace_Element(k,0,z_in_std(k));
      }

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > tmp1 = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    sol_op_lofi_->Multiply(*tmp1,*b);
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > tmp2 = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    M_->Multiply(*tmp2,*tmp1);
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > tmp3 = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    sol_op_lofi_->Multiply(*tmp3,*tmp2,true);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > tmp4 = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    M_->Multiply(*tmp4,*b);

    for(int k = 0; k < m_; k++)
      {
	z_out_std.Replace_Element(k,(*tmp3)(k,0)+reg_coeff_*(*tmp4)(k,0));
      }					
  }

  void Misfit_Gradient(HDSA::Vector<RealT> & u_grad, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z) const 
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > v = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    const Std_Vector<RealT> u_std = dynamic_cast<const Std_Vector<RealT>&>(u);
    Std_Vector<RealT> u_grad_std = dynamic_cast<Std_Vector<RealT>&>(u_grad);
    for(int k = 0; k < m_; k++)
      {
	v->Replace_Element(k,0,u_std(k)-(*target_)(k,0));
      }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > grad = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    M_->Multiply(*grad,*v);
    for(int k = 0; k < m_; k++)
      {
	u_grad_std.Replace_Element(k,(*grad)(k,0));
      }
  }

  void Apply_Misfit_Hessian(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z) const
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > v = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    const Std_Vector<RealT> u_in_std = dynamic_cast<const Std_Vector<RealT>&>(u_in);
    Std_Vector<RealT> u_out_std = dynamic_cast<Std_Vector<RealT>&>(u_out);
    for(int k = 0; k < m_; k++)
      {
	v->Replace_Element(k,0,u_in_std(k));
      }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Hv = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    M_->Multiply(*Hv,*v);
    for(int k = 0; k < m_; k++)
      {
	u_out_std.Replace_Element(k,(*Hv)(k,0));
      }
  }

  HDSA::Ptr<HDSA::Vector<RealT> > Load_Optimal_u( ) const 
  {
    HDSA::Ptr<Std_Vector<RealT> > u_opt = HDSA::makePtr<Std_Vector<RealT> >(m_);

    RealT val = 0.0;
    // read in data
    std::ifstream in("u_opt.txt");           
    // read the elements in the file into a vector  
    // test file open   
    if (in) {   
      for(int i = 0; i < m_; i++)
	{
	  in >> val;
	  u_opt->Replace_Element(i,val);
	}   
    }
    else
      {
	std::cout << "Error loading the data from u_opt.txt" << std::endl;
      }  
    return u_opt;
  }

  HDSA::Ptr<HDSA::Vector<RealT> > Load_Optimal_z( ) const 
  {
    HDSA::Ptr<Std_Vector<RealT> > z_opt = HDSA::makePtr<Std_Vector<RealT> >(m_);

    RealT val = 0.0;
    // read in data
    std::ifstream in("z_opt.txt");           
    // read the elements in the file into a vector  
    // test file open   
    if (in) {   
      for(int i = 0; i < m_; i++)
	{
	  in >> val;
	  z_opt->Replace_Element(i,val);
	}   
    }
    else
      {
	std::cout << "Error loading the data from z_opt.txt" << std::endl;
      }  
    return z_opt;
  }

  HDSA::Ptr<HDSA::MultiVector<RealT> > Load_Z_Data( ) const 
  {
    HDSA::Ptr<Std_Vector<RealT> > z = HDSA::makePtr<Std_Vector<RealT> >(m_);
    HDSA::Ptr<HDSA::MultiVector<RealT> > Z = HDSA::makePtr<HDSA::MultiVector<RealT> >(2,*z);

    HDSA::Ptr<HDSA::Vector<RealT> > z0 = (*Z)[0];
    HDSA::Ptr<HDSA::Vector<RealT> > z1 = (*Z)[1];
    Std_Vector<RealT> z0_std = dynamic_cast<Std_Vector<RealT>&>(*z0);
    Std_Vector<RealT> z1_std = dynamic_cast<Std_Vector<RealT>&>(*z1);

    // read in data
    std::ifstream in_Z("Z.txt");           
    // read the elements in the file into a vector  
    // test file open   
    RealT val;
    if (in_Z) {   
      for(int i = 0; i < m_; i++)
	{
	  in_Z >> val;
	  z0_std.Replace_Element(i,val);

	  in_Z >> val;
	  z1_std.Replace_Element(i,val);
	}   
    }
    else
      {
	std::cout << "Error loading the data from Z.txt" << std::endl;
      } 

    return Z;
  }

  HDSA::Ptr<HDSA::MultiVector<RealT> > Load_Y_Data( ) const 
  {
    HDSA::Ptr<Std_Vector<RealT> > y = HDSA::makePtr<Std_Vector<RealT> >(m_);
    HDSA::Ptr<HDSA::MultiVector<RealT> > Y = HDSA::makePtr<HDSA::MultiVector<RealT> >(2,*y);

    HDSA::Ptr<HDSA::Vector<RealT> > y0 = (*Y)[0];
    HDSA::Ptr<HDSA::Vector<RealT> > y1 = (*Y)[1];
    Std_Vector<RealT> y0_std = dynamic_cast<Std_Vector<RealT>&>(*y0);
    Std_Vector<RealT> y1_std = dynamic_cast<Std_Vector<RealT>&>(*y1);

    // read in data
    std::ifstream in_Y("Y.txt");           
    // read the elements in the file into a vector  
    // test file open   
    RealT val;
    if (in_Y) {   
      for(int i = 0; i < m_; i++)
	{
	  in_Y >> val;
	  y0_std.Replace_Element(i,val);

	  in_Y >> val;
	  y1_std.Replace_Element(i,val);
	}   
    }
    else
      {
	std::cout << "Error loading the data from Y.txt" << std::endl;
      } 

    return Y;
  }

  HDSA::Ptr<HDSA::Vector<RealT> > Load_Matlab_z_Update( ) const 
  {
    HDSA::Ptr<Std_Vector<RealT> > z_opt = HDSA::makePtr<Std_Vector<RealT> >(m_);

    RealT val = 0.0;
    // read in data
     std::ifstream in("z_update_matlab_solution.txt");           
    // read the elements in the file into a vector  
    // test file open   
    if (in) {   
      for(int i = 0; i < m_; i++)
	{
	  in >> val;
	  z_opt->Replace_Element(i,val);
	}   
    }
    else
      {
	std::cout << "Error loading the data from z_update_matlab_solution.txt" << std::endl;
      }  
    return z_opt;
  }

};

#endif


