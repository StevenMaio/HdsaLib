#ifndef HDSA_CDR_MODEL_ERROR_TEST_FS_OBJECTIVE_HPP
#define HDSA_CDR_MODEL_ERROR_TEST_FS_OBJECTIVE_HPP

template <class RealT>
class CDR_model_error_test_FS_Objective: public HDSA::FS_Objective<RealT> {

private:
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > D_; 
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > R_;
  HDSA::Ptr<HDSA::Vector<RealT> > data_;
  RealT beta_;

public:

  CDR_model_error_test_FS_Objective( HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & D, HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & R, HDSA::Ptr<HDSA::Vector<RealT> > & data, RealT beta):
    D_(D), R_(R), data_(data), beta_(beta)
  { }

  ~CDR_model_error_test_FS_Objective()
  { }

  // evaluate objective function
  RealT value(const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    RealT val = 0.0;

    HDSA::Ptr<HDSA::Vector<RealT> > m = u.Clone();
    m->set(u);
    m->axpy(-1.0,*data_);
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > D_matvec = D_->Multiply(*m,false);
    HDSA::Ptr<HDSA::Vector<RealT> > Dm = u.Clone();
    D_matvec->Write_Column_to_Vector(0,*Dm);
    val += 0.5*Dm->dot(*m);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > R_matvec = R_->Multiply(z,false);
    HDSA::Ptr<HDSA::Vector<RealT> > Rz = z.Clone();
    R_matvec->Write_Column_to_Vector(0,*Rz);
    val += 0.5*beta_*Rz->dot(*Rz);

    return val;    
  }

  // evaluate the gradient with respect to u
  void gradient_u(HDSA::Vector<RealT> & grad, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    HDSA::Ptr<HDSA::Vector<RealT> > m = u.Clone();
    m->set(u);
    m->axpy(-1.0,*data_);
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > D_matvec = D_->Multiply(*m,false);
    D_matvec->Write_Column_to_Vector(0,grad);
  }

  // evaluate the gradient with respect to z
  void gradient_z(HDSA::Vector<RealT> & grad, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > R_matvec = R_->Multiply(z,false);
    HDSA::Ptr<HDSA::Vector<RealT> > Rz = z.Clone();
    R_matvec->Write_Column_to_Vector(0,*Rz);
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > RtR_matvec = R_->Multiply(*Rz,true);
    RtR_matvec->Write_Column_to_Vector(0,grad);
  }
  
};


#endif
