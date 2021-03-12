#ifndef HDSA_CDR_MODEL_ERROR_TEST_RS_OBJECTIVE_HPP
#define HDSA_CDR_MODEL_ERROR_TEST_RS_OBJECTIVE_HPP

template <class RealT>
class CDR_model_error_test_RS_Objective: public HDSA::RS_Objective<RealT> {

private:
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > A_;
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > B_;
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > D_; 
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > R_;
  HDSA::Ptr<HDSA::Vector<RealT> > data_;
  RealT beta_;

public:

  CDR_model_error_test_RS_Objective(HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & A, HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & B, HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & D, 
				    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & R, HDSA::Ptr<HDSA::Vector<RealT> > & data, RealT beta): A_(A), B_(B), D_(D), R_(R), data_(data), beta_(beta)
  { 

  }

  ~CDR_model_error_test_RS_Objective()
  { }

  // evaluate objective function
  RealT value(const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    RealT val = 0.0;

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Bz_matvec = B_->Multiply(z,false);
    HDSA::Ptr<HDSA::Vector<RealT> > Bz = data_->Clone();
    Bz_matvec->Write_Column_to_Vector(0,Bz);

    HDSA::Ptr<HDSA::Vector<RealT> > r_tmp = data_->Clone();
    r_tmp->basis(r_tmp->dimension()-1);
    Bz->plus(*r_tmp);

    HDSA::Ptr<HDSA::Linear_Operator<RealT> > A_Op = HDSA::makePtr<Jac_Operator<RealT> >(A_);
    HDSA::Ptr<HDSA::Vector<RealT> > Ainv_Bz = data_->Clone();
    HDSA::Linear_Algebra::Iterative_Linear_Solve<RealT>(*Ainv_Bz,*Bz,A_Op,1.e-12,"GMRES");

    HDSA::Ptr<HDSA::Vector<RealT> > m = data_->Clone();
    m->set(*Ainv_Bz);
    m->axpy(-1.0,*data_);
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > D_matvec = D_->Multiply(*m,false);
    HDSA::Ptr<HDSA::Vector<RealT> > m_tmp = data_->Clone();
    D_matvec->Write_Column_to_Vector(0,*m_tmp);

    val += 0.5*m_tmp->dot(*m);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > R_matvec = R_->Multiply(z,false);
    HDSA::Ptr<HDSA::Vector<RealT> > Rz = z.Clone();
    R_matvec->Write_Column_to_Vector(0,*Rz);
    val += 0.5*beta_*Rz->dot(*Rz);

    return val;
  }

  // evaluate the gradient with respect to z
  void gradient_z(HDSA::Vector<RealT> & grad, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true) 
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Bz_matvec = B_->Multiply(z,false);
    HDSA::Ptr<HDSA::Vector<RealT> > Bz = data_->Clone();
    Bz_matvec->Write_Column_to_Vector(0,Bz);

    HDSA::Ptr<HDSA::Vector<RealT> > r_tmp = data_->Clone();
    r_tmp->basis(r_tmp->dimension()-1);
    Bz->plus(*r_tmp);

    HDSA::Ptr<HDSA::Linear_Operator<RealT> > A_Op = HDSA::makePtr<Jac_Operator<RealT> >(A_);
    HDSA::Ptr<HDSA::Vector<RealT> > Ainv_Bz = data_->Clone();
    HDSA::Linear_Algebra::Iterative_Linear_Solve<RealT>(*Ainv_Bz,*Bz,A_Op,1.e-12,"GMRES");

    HDSA::Ptr<HDSA::Vector<RealT> > m = data_->Clone();
    m->set(*Ainv_Bz);
    m->axpy(-1.0,*data_);
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > D_matvec = D_->Multiply(*m,false);
    HDSA::Ptr<HDSA::Vector<RealT> > m_tmp = data_->Clone();
    D_matvec->Write_Column_to_Vector(0,*m_tmp);

    HDSA::Ptr<HDSA::Linear_Operator<RealT> > At_Op = HDSA::makePtr<Jac_Trans_Operator<RealT> >(A_);
    HDSA::Ptr<HDSA::Vector<RealT> > Atinv_v = data_->Clone();
    HDSA::Linear_Algebra::Iterative_Linear_Solve<RealT>(*Atinv_v,*m_tmp,At_Op,1.e-12,"GMRES");

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Btv_matvec = B_->Multiply(*Atinv_v,true);
    Btv_matvec->Write_Column_to_Vector(0,grad);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > R_matvec = R_->Multiply(z,false);
    HDSA::Ptr<HDSA::Vector<RealT> > Rz = z.Clone();
    R_matvec->Write_Column_to_Vector(0,*Rz);
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > RtR_matvec = R_->Multiply(*Rz,true);
    HDSA::Ptr<HDSA::Vector<RealT> > grad_reg = grad.Clone();
    RtR_matvec->Write_Column_to_Vector(0,*grad_reg);

    grad.axpy(beta_,*grad_reg);
  } 
 

  // Overload HDSA::Linear_Operator to take matrix vector products for the state jacobian solve
  template <class ScalarType>
  class Jac_Operator : public HDSA::Linear_Operator<ScalarType>
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > A_;
    
    public:
    
    Jac_Operator(HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & A)
    { 
      A_ = A;
    }
    
    //! Dtor
    ~Jac_Operator()
    {}
      
    void matvec(HDSA::Ptr<HDSA::Vector<RealT> > & y, const HDSA::Ptr<HDSA::Vector<RealT> > & x) const 
    {
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Ax_matvec = A_->Multiply(x,false);
      Ax_matvec->Write_Column_to_Vector(0,*y);
    }
      
  };

  // Overload HDSA::Linear_Operator to take matrix vector products for the state jacobian transpose solve
  template <class ScalarType>
  class Jac_Trans_Operator : public HDSA::Linear_Operator<ScalarType>
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > A_;
    
    public:
      
    Jac_Trans_Operator(HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & A)
    { 
      A_ = A;
    }
    
    //! Dtor
    ~Jac_Trans_Operator()
    {}
      
    void matvec(HDSA::Ptr<HDSA::Vector<RealT> > & y, const HDSA::Ptr<HDSA::Vector<RealT> > & x) const 
    {
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Atx_matvec = A_->Multiply(x,true);
      Atx_matvec->Write_Column_to_Vector(0,*y);
    }
      
  };

};


#endif
