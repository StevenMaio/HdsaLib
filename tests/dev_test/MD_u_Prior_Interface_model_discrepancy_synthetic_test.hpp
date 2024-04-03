#ifndef HDSA_MD_U_PRIOR_INTERFACE_MODEL_DISCREPANCY_SYNTHETIC_TEST_HPP
#define HDSA_MD_U_PRIOR_INTERFACE_MODEL_DISCREPANCY_SYNTHETIC_TEST_HPP

template <class RealT>
class MD_u_Prior_Interface_model_discrepancy_synthetic_test : public HDSA::MD_u_Prior_Interface<RealT> {

private:
  int m_; // Mesh resolution                                                                                                                                                                                                    
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > x_; // Mesh nodes on [0,1]                                                                                                                                                              
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > S_; // Stiffness matrix                                                                                                                                                                  
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > M_; // Mass matrix                                                                                                                                                                       
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > W_u_; // State weighting matrix                                                                                                                                                         

public:
  MD_u_Prior_Interface_model_discrepancy_synthetic_test()
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

    W_u_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > I = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);
    for(int k = 0; k < m_; k++)
      {
        I->Replace_Element(k,k,1.0);
      }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Minv = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);
    HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*M_,*Minv,*I);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > E_u = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);
    for(int i = 0; i < m_; i++)
      {
        for(int j = 0; j < m_; j++)
          {
            RealT val = 2.0*( (5.e-2)*(*S_)(i,j) + (*M_)(i,j) );
            E_u->Replace_Element(i,j,val);
	  }
      }

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > tmp = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);
    Minv->Multiply(*tmp,*E_u);
    E_u->Multiply(*W_u_,*tmp);
  }

  virtual ~MD_u_Prior_Interface_model_discrepancy_synthetic_test()
  { }

  void Apply_M_u(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const 
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

  void Apply_W_u_Plus_scalar_M_u_Inverse(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in, const RealT & beta) const 
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > b = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    const Std_Vector<RealT>& u_in_std = dynamic_cast<const Std_Vector<RealT>&>(u_in);
    Std_Vector<RealT>& u_out_std = dynamic_cast<Std_Vector<RealT>&>(u_out);
    for(int k = 0; k < m_; k++)
      {
	b->Replace_Element(k,0,u_in_std(k));
      }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > x = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Wu_beta_Mu = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);
    for(int i = 0; i < m_; i++)
      {
	for(int j = 0; j < m_; j++)
	  {
	    Wu_beta_Mu->Replace_Element(i,j,(*W_u_)(i,j) + beta*(*M_)(i,j));
	  }
      }
    HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*Wu_beta_Mu,*x_,*b);
    for(int k = 0; k < m_; k++)
      {
	u_out_std.Replace_Element(k,(*x)(k,0));
      }
  }

  void Apply_W_u_Inverse(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const 
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > b = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    const Std_Vector<RealT>& u_in_std = dynamic_cast<const Std_Vector<RealT>&>(u_in);
    Std_Vector<RealT>& u_out_std = dynamic_cast<Std_Vector<RealT>&>(u_out);
    for(int k = 0; k < m_; k++)
      {
        b->Replace_Element(k,0,u_in_std(k));
      }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > x = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*W_u_,*x_,*b);
    for(int k = 0; k < m_; k++)
      {
        u_out_std.Replace_Element(k,(*x)(k,0));
      }
  }


  };

#endif

