#ifndef ELLIPTIC_U_PRIOR_INTERFACE_SIMOPT_TEST_PROBLEM_HPP
#define ELLIPTIC_U_PRIOR_INTERFACE_SIMOPT_TEST_PROBLEM_HPP


template <class RealT>
class Elliptic_u_Prior_Interface_SimOptTestProb : public HDSA::MD_Elliptic_u_Prior_Interface<RealT> {
  
private:
  int m_; // Mesh resolution
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > S_; // Stiffness matrix
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > M_; // Mass matrix
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > E_u_; // State elliptic operator

  
public:
  Elliptic_u_Prior_Interface_SimOptTestProb(RealT & alpha_u, RealT & beta_u, int & m): HDSA::MD_Elliptic_u_Prior_Interface<RealT>(alpha_u)
  {
    m_ = m;

    RealT h = 1.0/static_cast<RealT>(m_-1);

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

    E_u_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);
    for(int i = 0; i < m_; i++)
      {
	for(int j = 0; j < m_; j++)
	  {
	    RealT val = beta_u*(*S_)(i,j) + (*M_)(i,j);
	    E_u_->Replace_Element(i,j,val);
	  }
      }
  }
  
  virtual ~Elliptic_u_Prior_Interface_SimOptTestProb()
  { }                                    
  
  void Apply_E_u_Inverse(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const
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
    HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*E_u_,*x,*b);
    for(int k = 0; k < m_; k++)
      {
    	(*u_out_std)[k] = (*x)(k,0);
      }
 
  }
  
  void Apply_E_u_Inverse_Transpose(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const
  {
    Apply_E_u_Inverse(u_out,u_in);
  }

  void Apply_M_u(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const
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

 

  };


#endif

