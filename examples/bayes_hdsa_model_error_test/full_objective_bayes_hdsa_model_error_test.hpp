#ifndef HDSA_BAYES_HDSA_MODEL_ERROR_TEST_FS_OBJECTIVE_HPP
#define HDSA_BAYES_HDSA_MODEL_ERROR_TEST_FS_OBJECTIVE_HPP

template <class RealT>
class bayes_hdsa_model_error_test_FS_Objective: public HDSA::FS_Objective<RealT> {

private:
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > D_;
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > W_;
  HDSA::Ptr<HDSA::Vector<RealT> > data_;
  RealT beta_;

public:
  
  bayes_hdsa_model_error_test_FS_Objective(HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & D, HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & W, HDSA::Ptr<HDSA::Vector<RealT> > & data, RealT & beta): 
    D_(D), W_(W), data_(data), beta_(beta)
  { }

  ~bayes_hdsa_model_error_test_FS_Objective()
  { }

  // evaluate objective function
  RealT value(const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    RealT val = 0.0;

    HDSA::Ptr<HDSA::Vector<RealT> > m = u.Clone();
    m->set(u);
    m->axpy(-1.0,*data_);
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > W_matvec = W_->Multiply(*m,false);
    HDSA::Ptr<HDSA::Vector<RealT> > Wm = u.Clone();
    W_matvec->Write_Column_to_Vector(0,*Wm);
    val += 0.5*Wm->dot(*m);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > D_matvec = D_->Multiply(z,false);
    HDSA::Ptr<HDSA::Vector<RealT> > Dz = z.Clone();
    D_matvec->Write_Column_to_Vector(0,*Dz);
    val += 0.5*beta_*z.dot(*Dz);

    return val;    
  }

  // evaluate the gradient with respect to u
  void gradient_u(HDSA::Vector<RealT> & grad, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    HDSA::Ptr<HDSA::Vector<RealT> > m = u.Clone();
    m->set(u);
    m->axpy(-1.0,*data_);
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > W_matvec = W_->Multiply(*m,false);
    W_matvec->Write_Column_to_Vector(0,grad);
  }

  // evaluate the gradient with respect to z
  void gradient_z(HDSA::Vector<RealT> & grad, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > D_matvec = D_->Multiply(z,false);
    HDSA::Ptr<HDSA::Vector<RealT> > Dz = z.Clone();
    D_matvec->Write_Column_to_Vector(0,*Dz);
    Dz->scale(beta_);
  }
  
};


#endif
