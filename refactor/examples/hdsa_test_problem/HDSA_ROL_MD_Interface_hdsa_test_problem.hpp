#ifndef HDSA_ROL_MD_INTERFACE_HDSA_TEST_PROBLEM_HPP
#define HDSA_ROL_MD_INTERFACE_HDSA_TEST_PROBLEM_HPP

template <class RealT>
class ROL_Model_Discrepancy_Interface_hdsa_test_problem : public HDSA::ROL_Model_Discrepancy_Interface<RealT> {

private:
  int m_; // Mesh resolution
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > x_; // Mesh nodes on [0,1]
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > S_; // Stiffness matrix
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > M_; // Mass matrix
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Minv_; // Mass matrix inverse
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > E_L_; // State elliptic operator
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Gamma_; // Control weighting matrix
  
public:

  ROL_Model_Discrepancy_Interface_hdsa_test_problem(ROL::Ptr<ROL::Objective_SimOpt<RealT> > & obj_simopt, ROL::Ptr<ROL::Constraint_SimOpt<RealT> > & con_simopt, 
						    ROL::Ptr<ROL::Vector<RealT> > & u, ROL::Ptr<ROL::Vector<RealT> > & z, int m): 
    HDSA::ROL_Model_Discrepancy_Interface<RealT>(obj_simopt,con_simopt,u,z) 
  {  
    m_ = m;

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

    Gamma_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > I = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);
    for(int k = 0; k < m_; k++)
      {
	I->Replace_Element(k,k,1.0);
      }
    Minv_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);
    HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*M_,*Minv_,*I); 

    E_L_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > E_Gamma = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);
    for(int i = 0; i < m_; i++)
      {
	for(int j = 0; j < m_; j++)
	  {
	    RealT val = 2.0*( (5.e-2)*(*S_)(i,j) + (*M_)(i,j) );
	    E_L_->Replace_Element(i,j,val);
	    val = (1.e2)*( (1.e-2)*(*S_)(i,j) + (*M_)(i,j) );
	    E_Gamma->Replace_Element(i,j,val);
	  }
      }

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > tmp = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);
    Minv_->Multiply(*tmp,*E_Gamma);
    E_Gamma->Multiply(*Gamma_,*tmp);

    int num_sing_vals = 50;
    int oversampling = 1;
    int num_subspace_iters = 2;
    HDSA::Ptr<HDSA::Vector<RealT> > u_vec = HDSA::makePtr<HDSA::ROL_Vector<RealT> >(*u);
    HDSA::Model_Discrepancy_Interface_Elliptic_Prior<RealT>::Compute_Elliptic_GSVD(num_sing_vals,oversampling,num_subspace_iters,*u_vec);

  }

  virtual ~ROL_Model_Discrepancy_Interface_hdsa_test_problem()
  { }

  void Apply_u_Elliptic_Operator_Inverse(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const 
  {
    HDSA::ROL_Vector<RealT>& u_out_rol = dynamic_cast<HDSA::ROL_Vector<RealT>&>(u_out);
    ROL::Ptr<std::vector<RealT> > u_out_std = dynamic_cast<ROL::StdVector<RealT>&>(*u_out_rol.rol_vec).getVector();

    const HDSA::ROL_Vector<RealT>& u_in_rol = dynamic_cast<const HDSA::ROL_Vector<RealT>&>(u_in);
    ROL::Ptr<const std::vector<RealT> > u_in_std = dynamic_cast<const ROL::StdVector<RealT>&>(*u_in_rol.rol_vec).getVector();

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > b = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    for(int k = 0; k < m_; k++)
      {
	b->Replace_Element(k,0,(*u_in_std)[k]);
      }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > x = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*E_L_,*x,*b);
    for(int k = 0; k < m_; k++)
      {
	(*u_out_std)[k] = (*x)(k,0);
      }
  }

  void Apply_u_Elliptic_Operator_Inverse_Transpose(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const 
  {
    HDSA::ROL_Vector<RealT>& u_out_rol = dynamic_cast<HDSA::ROL_Vector<RealT>&>(u_out);
    ROL::Ptr<std::vector<RealT> > u_out_std = dynamic_cast<ROL::StdVector<RealT>&>(*u_out_rol.rol_vec).getVector();

    const HDSA::ROL_Vector<RealT>& u_in_rol = dynamic_cast<const HDSA::ROL_Vector<RealT>&>(u_in);
    ROL::Ptr<const std::vector<RealT> > u_in_std = dynamic_cast<const ROL::StdVector<RealT>&>(*u_in_rol.rol_vec).getVector();

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > b = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    for(int k = 0; k < m_; k++)
      {
    	b->Replace_Element(k,0,(*u_in_std)[k]);
      }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > x = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*E_L_,*x,*b);
    for(int k = 0; k < m_; k++)
      {
    	(*u_out_std)[k] = (*x)(k,0);
      }
  }

  void Apply_u_Mass_Mat(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const 
  {
    HDSA::ROL_Vector<RealT>& u_out_rol = dynamic_cast<HDSA::ROL_Vector<RealT>&>(u_out);
    ROL::Ptr<std::vector<RealT> > u_out_std = dynamic_cast<ROL::StdVector<RealT>&>(*u_out_rol.rol_vec).getVector();

    const HDSA::ROL_Vector<RealT>& u_in_rol = dynamic_cast<const HDSA::ROL_Vector<RealT>&>(u_in);
    ROL::Ptr<const std::vector<RealT> > u_in_std = dynamic_cast<const ROL::StdVector<RealT>&>(*u_in_rol.rol_vec).getVector();

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > b = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    for(int k = 0; k < m_; k++)
      {
    	b->Replace_Element(k,0,(*u_in_std)[k]);
      }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > x = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    M_->Multiply(*x,*b);
    for(int k = 0; k < m_; k++)
      {
    	(*u_out_std)[k] = (*x)(k,0);
      }
  }

  void Apply_u_Mass_Mat_Inverse(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const 
  {
    HDSA::ROL_Vector<RealT>& u_out_rol = dynamic_cast<HDSA::ROL_Vector<RealT>&>(u_out);
    ROL::Ptr<std::vector<RealT> > u_out_std = dynamic_cast<ROL::StdVector<RealT>&>(*u_out_rol.rol_vec).getVector();

    const HDSA::ROL_Vector<RealT>& u_in_rol = dynamic_cast<const HDSA::ROL_Vector<RealT>&>(u_in);
    ROL::Ptr<const std::vector<RealT> > u_in_std = dynamic_cast<const ROL::StdVector<RealT>&>(*u_in_rol.rol_vec).getVector();

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > b = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    for(int k = 0; k < m_; k++)
      {
    	b->Replace_Element(k,0,(*u_in_std)[k]);
      }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > x = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    Minv_->Multiply(*x,*b);
    for(int k = 0; k < m_; k++)
      {
    	(*u_out_std)[k] = (*x)(k,0);
      }
  }

  // Manipulate prior covariances using direct linear algebra since the dimension is small
  void Apply_Gamma_Mat_Inverse(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const 
  {
    HDSA::ROL_Vector<RealT>& z_out_rol = dynamic_cast<HDSA::ROL_Vector<RealT>&>(z_out);
    ROL::Ptr<std::vector<RealT> > z_out_std = dynamic_cast<ROL::StdVector<RealT>&>(*z_out_rol.rol_vec).getVector();

    const HDSA::ROL_Vector<RealT>& z_in_rol = dynamic_cast<const HDSA::ROL_Vector<RealT>&>(z_in);
    ROL::Ptr<const std::vector<RealT> > z_in_std = dynamic_cast<const ROL::StdVector<RealT>&>(*z_in_rol.rol_vec).getVector();

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > b = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    for(int k = 0; k < m_; k++)
      {
    	b->Replace_Element(k,0,(*z_in_std)[k]);
      }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > x = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*Gamma_,*x,*b);
    for(int k = 0; k < m_; k++)
      {
    	(*z_out_std)[k] = (*x)(k,0);
      }
  }

  HDSA::Ptr<HDSA::Vector<RealT> > Load_Optimal_u(void) const 
  {
    ROL::Ptr<std::vector<RealT> > u_ptr    = ROL::makePtr<std::vector<RealT> >(m_, 0.0);
    ROL::StdVector<RealT> u(u_ptr);
    ROL::Ptr<ROL::Vector<RealT> > up  = ROL::makePtrFromRef(u);
    for (int i=0;i<m_;i++){
      (*u_ptr)[i] = std::pow(1.0 + (*x_)(i,0),3.0);
    }
    HDSA::Ptr<HDSA::Vector<RealT> > u_opt = HDSA::makePtr<HDSA::ROL_Vector<RealT> >(*up);
    HDSA::ROL_Vector<RealT>& u_opt_rol = dynamic_cast<HDSA::ROL_Vector<RealT>&>(*u_opt);
    u_opt_rol.rol_vec->set(*up);
    return u_opt;
  }

  HDSA::Ptr<HDSA::Vector<RealT> > Load_Optimal_z(void) const 
  {
    ROL::Ptr<std::vector<RealT> > z_ptr    = ROL::makePtr<std::vector<RealT> >(m_, 0.0);
    ROL::StdVector<RealT> z(z_ptr);
    ROL::Ptr<ROL::Vector<RealT> > zp  = ROL::makePtrFromRef(z);
    for (int i=0;i<m_;i++){
      (*z_ptr)[i] = 1.0 + (*x_)(i,0);
    }
    HDSA::Ptr<HDSA::Vector<RealT> > z_opt = HDSA::makePtr<HDSA::ROL_Vector<RealT> >(*zp);
    HDSA::ROL_Vector<RealT>& z_opt_rol = dynamic_cast<HDSA::ROL_Vector<RealT>&>(*z_opt);
    z_opt_rol.rol_vec->set(*zp);
    return z_opt;
  }

  HDSA::Ptr<HDSA::MultiVector<RealT> > Load_Z_Data(void) const 
  {
    HDSA::Ptr<HDSA::Vector<RealT> > z_opt = Load_Optimal_z();
    HDSA::Ptr<HDSA::MultiVector<RealT> > Z = HDSA::makePtr<HDSA::MultiVector<RealT> >(2,*z_opt);
    (*Z)[0]->set(*z_opt);

    ROL::Ptr<std::vector<RealT> > z_ptr    = ROL::makePtr<std::vector<RealT> >(m_, 0.0);
    ROL::StdVector<RealT> z(z_ptr);
    ROL::Ptr<ROL::Vector<RealT> > zp  = ROL::makePtrFromRef(z);
    for (int i=0;i<m_;i++){
      (*z_ptr)[i] = (*x_)(i,0) + std::pow((*x_)(i,0),2.0);
    }
    HDSA::ROL_Vector<RealT>& z_opt_rol = dynamic_cast<HDSA::ROL_Vector<RealT>&>(*z_opt);
    z_opt_rol.rol_vec->set(*zp);
    (*Z)[1]->set(*z_opt);

    return Z;
  }

  HDSA::Ptr<HDSA::MultiVector<RealT> > Load_Y_Data(void) const 
  {
    ROL::Ptr<std::vector<RealT> > u_ptr    = ROL::makePtr<std::vector<RealT> >(m_, 0.0);
    ROL::StdVector<RealT> u(u_ptr);
    ROL::Ptr<ROL::Vector<RealT> > up  = ROL::makePtrFromRef(u);
    for (int i=0;i<m_;i++){
      (*u_ptr)[i] = 0.2*std::pow(1.0+(*x_)(i,0),2.0);
    }
    HDSA::Ptr<HDSA::Vector<RealT> > u_hdsa = HDSA::makePtr<HDSA::ROL_Vector<RealT> >(*up);
    HDSA::ROL_Vector<RealT>& u_hdsa_rol = dynamic_cast<HDSA::ROL_Vector<RealT>&>(*u_hdsa);
    u_hdsa_rol.rol_vec->set(*up);
    HDSA::Ptr<HDSA::MultiVector<RealT> > Y = HDSA::makePtr<HDSA::MultiVector<RealT> >(2,*u_hdsa);
    (*Y)[0]->set(*u_hdsa);

    for (int i=0;i<m_;i++){
      (*u_ptr)[i] = 0.2*std::pow((*x_)(i,0)+std::pow((*x_)(i,0),2.0),2.0);
    }
    u_hdsa_rol.rol_vec->set(*up);
    (*Y)[1]->set(*u_hdsa);

    return Y;
  }

};

#endif


