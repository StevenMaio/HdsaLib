#ifndef HDSA_CDR_MODEL_ERROR_TEST_CONSTRAINT_HPP
#define HDSA_CDR_MODEL_ERROR_TEST_CONSTRAINT_HPP

template <class RealT>
class CDR_model_error_test_Constraint: public HDSA::Constraint<RealT> {

private:
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > A_;
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > B_;

public:

  CDR_model_error_test_Constraint(HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & A, HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & B): A_(A), B_(B)
  { }

  ~CDR_model_error_test_Constraint()
  { }

  // evaluate constraint residual
  void value(HDSA::Vector<RealT> & r, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Au_matvec = A_->Multiply(u,false);
    HDSA::Ptr<HDSA::Vector<RealT> > Au = u.Clone();
    Au_matvec->Write_Column_to_Vector(0,*Au);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Bz_matvec = B_->Multiply(z,false);
    HDSA::Ptr<HDSA::Vector<RealT> > Bz = u.Clone();
    Bz_matvec->Write_Column_to_Vector(0,*Bz);

    HDSA::Ptr<HDSA::Vector<RealT> > r_tmp = r.Clone();
    r_tmp->basis(r_tmp->dimension()-1);
    Bz->plus(*r_tmp);

    r.set(*Au);
    r.axpy(-1.0,*Bz);
  }

  // evaluate the jacobian with respect to u vector product
  void jacobian_u(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, 
		  const HDSA::Vector<RealT> & theta, const bool & update = true) 
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Av_matvec = A_->Multiply(v,false);
    Av_matvec->Write_Column_to_Vector(0,Jv);
  } 

  // evaluate the inverse jacobian with respect to u vector product
  void jacobian_u_inverse(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, 
			  const HDSA::Vector<RealT> & theta, const bool & update = true) 
  {
    HDSA::Ptr<HDSA::Linear_Operator<RealT> > A_Op = HDSA::makePtr<Jac_Operator<RealT> >(A_);
    HDSA::Linear_Algebra::Iterative_Linear_Solve<RealT>(Jv,v,A_Op,1.e-12,"GMRES");
  } 

  // evaluate the jacobian with respect to z vector product
  void jacobian_z(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, 
		  const HDSA::Vector<RealT> & theta, const bool & update = true) 
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Bv_matvec = B_->Multiply(v,false);
    Bv_matvec->Write_Column_to_Vector(0,Jv);
    Jv.scale(-1.0);
  } 

  // evaluate the adjoint jacobian with respect to u vector product
  void jacobian_u_adjoint(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, 
			  const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Atv_matvec = A_->Multiply(v,true);
    Atv_matvec->Write_Column_to_Vector(0,Jv);
  } 

  // evaluate the inverse adjoint jacobian with respect to u vector product
  void jacobian_u_adjoint_inverse(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, 
				  const HDSA::Vector<RealT> & theta, const bool & update = true) 
  {
    HDSA::Ptr<HDSA::Linear_Operator<RealT> > At_Op = HDSA::makePtr<Jac_Trans_Operator<RealT> >(A_);
    HDSA::Linear_Algebra::Iterative_Linear_Solve<RealT>(Jv,v,At_Op,1.e-12,"GMRES");
  }

  // evaluate the adjoint jacobian with respect to z vector product
  void jacobian_z_adjoint(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, 
			  const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Btv_matvec = B_->Multiply(v,true);
    Btv_matvec->Write_Column_to_Vector(0,Jv);
    Jv.scale(-1.0);
  }

  // evaluate the adjoint hessian with respect to u,u vector product, i.e. to compute second derivatives of lambda^T*c(u,z,theta) with respect to u,u
  void hessian_u_u_adjoint(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & lambda, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, 
			   const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true) 
  {
    Jv.zero();
  }

  // evaluate the adjoint hessian with respect to u,z vector product, i.e. to compute second derivatives of lambda^T*c(u,z,theta) with respect to u,z
  void hessian_u_z_adjoint(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & lambda, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, 
			   const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    Jv.zero();
  } 

  // evaluate the adjoint hessian with respect to z,z vector product, i.e. to compute second derivatives of lambda^T*c(u,z,theta) with respect to z,z
  void hessian_z_z_adjoint(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & lambda, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, 
			   const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    Jv.zero();
  }     

  // evaluate the adjoint hessian with respect to z,u vector product, i.e. to compute second derivatives of lambda^T*c(u,z,theta) with respect to z,u
  void hessian_z_u_adjoint(HDSA::Vector<RealT> & Jv, const HDSA::Vector<RealT> & lambda, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & u, 
			   const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true) 
  {
    Jv.zero();
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
