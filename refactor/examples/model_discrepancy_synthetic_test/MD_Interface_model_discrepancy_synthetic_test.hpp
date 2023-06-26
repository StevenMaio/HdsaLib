#ifndef HDSA_MD_INTERFACE_MODEL_DISCREPANCY_SYNTHETIC_TEST_HPP
#define HDSA_MD_INTERFACE_MODEL_DISCREPANCY_SYNTHETIC_TEST_HPP

template <class RealT>
class Model_Discrepancy_Interface_model_discrepancy_synthetic_test : public HDSA::Model_Discrepancy_Interface<RealT> {

private:
  int m_; // Mesh resolution
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > x_; // Mesh nodes on [0,1]
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > S_; // Stiffness matrix
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > M_; // Mass matrix
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Minv_; // Mass matrix inverse
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > L_; // State weighting matrix
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Gamma_; // Control weighting matrix

public:

  Model_Discrepancy_Interface_model_discrepancy_synthetic_test()
  {  
    m_ = 51;
    RealT h = 1.0/static_cast<RealT>(m_-1);
    x_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    for(int k = 0; k < m_; k++)
      {
	x_->Replace_Element(k,0,static_cast<RealT>(k)/static_cast<RealT>(m_-1));
      }

    S_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);
    M_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);

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

    L_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);
    Gamma_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > I = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);
    for(int k = 0; k < m_; k++)
      {
	I->Replace_Element(k,k,1.0);
      }
    Minv_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);
    HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*M_,*Minv_,*I); 

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > E_L = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > E_Gamma = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);
    for(int i = 0; i < m_; i++)
      {
	for(int j = 0; j < m_; j++)
	  {
	    RealT val = 2.0*( (5.e-2)*(*S_)(i,j) + (*M_)(i,j) );
	    E_L->Replace_Element(i,j,val);
	    val = (1.e2)*( (1.e-2)*(*S_)(i,j) + (*M_)(i,j) );
	    E_Gamma->Replace_Element(i,j,val);
	  }
      }

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > tmp = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);
    Minv_->Multiply(*tmp,*E_L);
    E_L->Multiply(*L_,*tmp);

    Minv_->Multiply(*tmp,*E_Gamma);
    E_Gamma->Multiply(*Gamma_,*tmp);
  }

  virtual ~Model_Discrepancy_Interface_model_discrepancy_synthetic_test()
  { }

  // Manipulate prior covariances using direct linear algebra since the dimension is small
  void Apply_Gamma_Mat_Inverse(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const 
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > b = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    const Std_Vector<RealT>& z_in_std = dynamic_cast<const Std_Vector<RealT>&>(z_in);
    Std_Vector<RealT>& z_out_std = dynamic_cast<Std_Vector<RealT>&>(z_out);
    for(int k = 0; k < m_; k++)
      {
	b->Replace_Element(k,0,z_in_std(k));
      }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > x = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*Gamma_,*x,*b);
    for(int k = 0; k < m_; k++)
      {
	z_out_std.Replace_Element(k,(*x)(k,0));
      }
  }

  void Apply_L_Plus_beta_Identity_Mat_Inverse(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in, const RealT beta) const 
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > A = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);
    for(int i = 0; i < m_; i++)
      {
	for(int j = 0; j < m_; j++)
	  {
	    A->Replace_Element(i,j,(*L_)(i,j));
	  }
	A->Replace_Element(i,i,(*L_)(i,i)+beta);
      }

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > b = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    const Std_Vector<RealT> u_in_std = dynamic_cast<const Std_Vector<RealT>&>(u_in);
    Std_Vector<RealT> u_out_std = dynamic_cast<Std_Vector<RealT>&>(u_out);
    for(int k = 0; k < m_; k++)
      {
	b->Replace_Element(k,0,u_in_std(k));
      }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > x = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*A,*x,*b);
    for(int k = 0; k < m_; k++)
      {
	u_out_std.Replace_Element(k,(*x)(k,0));
      }
  }

  void Apply_L_Mat_Inverse(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const 
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > b = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    const Std_Vector<RealT> u_in_std = dynamic_cast<const Std_Vector<RealT>&>(u_in);
    Std_Vector<RealT> u_out_std = dynamic_cast<Std_Vector<RealT>&>(u_out);
    for(int k = 0; k < m_; k++)
      {
	b->Replace_Element(k,0,u_in_std(k));
      }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > x = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*L_,*x,*b);
    for(int k = 0; k < m_; k++)
      {
	u_out_std.Replace_Element(k,(*x)(k,0));
      }
  }

  // Assume a constraint u = z^3 nodewise on the mesh defined by nodes in x_
  // Assume an objective (1/2)*(u-T)^t*M*(u-T) where T = (x_+1.0)^3 so that the optimal solution is u_opt=(x_+1.0)^3 and z_opt=x_+1.0
  // Assume a high-fidelity model u = z^3 + .2*z^2

  void Apply_Solution_Operator_z_Jacobian_Transpose(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & u_in, const HDSA::Vector<RealT> & z) const 
  {
    const Std_Vector<RealT> u_in_std = dynamic_cast<const Std_Vector<RealT>&>(u_in);
    const Std_Vector<RealT> z_std = dynamic_cast<const Std_Vector<RealT>&>(z);
    Std_Vector<RealT> z_out_std = dynamic_cast<const Std_Vector<RealT>&>(z_out);
    for(int k = 0; k < m_; k++)
      {
	z_out_std.Replace_Element(k,3.0*std::pow(z_std(k),2.0)*u_in_std(k));
      }
  }

  // This implementation assumes that it is evaluated at the optimal z so that the adjoint=0, a more general implementation would include a term multiplied by the adjoint variable
  void Apply_RS_Hessian(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in, const HDSA::Vector<RealT> & z) const 
  {
    const Std_Vector<RealT> z_std = dynamic_cast<const Std_Vector<RealT>&>(z);
    const Std_Vector<RealT> z_in_std = dynamic_cast<const Std_Vector<RealT>&>(z_in);
    Std_Vector<RealT> z_out_std = dynamic_cast<const Std_Vector<RealT>&>(z_out);
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > v = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    for(int k = 0; k < m_; k++)
      {
	v->Replace_Element(k,0,9.0*(z_in_std(k)*std::pow(z_std(k),2.0)));
      }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > M_v = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    M_->Multiply(*M_v,*v);
    for(int k = 0; k < m_; k++)
      {
	z_out_std.Replace_Element(k,(*M_v)(k,0)*std::pow(z_std(k),2.0));
      }						
  }

  void Misfit_Gradient(HDSA::Vector<RealT> & u_grad, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z) const 
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > v = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    const Std_Vector<RealT> u_std = dynamic_cast<const Std_Vector<RealT>&>(u);
    Std_Vector<RealT> u_grad_std = dynamic_cast<Std_Vector<RealT>&>(u_grad);
    for(int k = 0; k < m_; k++)
      {
	v->Replace_Element(k,0,u_std(k)-std::pow((*x_)(k,0)+1.0,3.0));
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
    for(int k = 0; k < m_; k++)
      {
	u_opt->Replace_Element(k,std::pow((*x_)(k,0)+1.0,3.0));
      }
    return u_opt;
  }

  HDSA::Ptr<HDSA::Vector<RealT> > Load_Optimal_z( ) const 
  {
    HDSA::Ptr<Std_Vector<RealT> > z_opt = HDSA::makePtr<Std_Vector<RealT> >(m_);
    for(int k = 0; k < m_; k++)
      {
	z_opt->Replace_Element(k,(*x_)(k,0)+1.0);
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

    for(int k = 0; k < m_; k++)
      {
	z0_std.Replace_Element(k,(*x_)(k,0)+1.0);
	z1_std.Replace_Element(k,(*x_)(k,0) + std::pow((*x_)(k,0),2.0));
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

    for(int k = 0; k < m_; k++)
      {
	y0_std.Replace_Element(k,0.2*std::pow((*x_)(k,0)+1.0,2.0));
	y1_std.Replace_Element(k,0.2*std::pow((*x_)(k,0) + std::pow((*x_)(k,0),2.0),2.0));
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


