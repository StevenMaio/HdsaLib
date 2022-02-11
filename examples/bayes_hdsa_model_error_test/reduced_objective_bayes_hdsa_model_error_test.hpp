#ifndef HDSA_BAYES_HDSA_MODEL_ERROR_TEST_RS_OBJECTIVE_HPP
#define HDSA_BAYES_HDSA_MODEL_ERROR_TEST_RS_OBJECTIVE_HPP

template <class RealT>
class bayes_hdsa_model_error_test_RS_Objective: public HDSA::RS_Objective<RealT> {

private:
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > D_; 
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > W_;
  HDSA::Ptr<HDSA::Vector<RealT> > data_;
  RealT beta_;

public:

  bayes_hdsa_model_error_test_RS_Objective(HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & D, HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & W, HDSA::Ptr<HDSA::Vector<RealT> > & data, RealT beta):
    D_(D), W_(W), data_(data), beta_(beta)
  { 

  }

  ~bayes_hdsa_model_error_test_RS_Objective()
  { }

  // evaluate objective function
  RealT value(const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    RealT val = 0.0;

    HDSA::Ptr<HDSA::Vector<RealT> > m = z.Clone();
    for(int k = 0; k < z.dimension(); k++)
      {
	m->Replace_Element(k,std::pow(z(k),3.0));
      }
    m->axpy(-1.0,*data_);
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > W_matvec = W_->Multiply(*m,false);
    HDSA::Ptr<HDSA::Vector<RealT> > Wm = z.Clone();
    W_matvec->Write_Column_to_Vector(0,*Wm);
    val += 0.5*Wm->dot(*m);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > D_matvec = D_->Multiply(z,false);
    HDSA::Ptr<HDSA::Vector<RealT> > Dz = z.Clone();
    D_matvec->Write_Column_to_Vector(0,*Dz);
    val += 0.5*beta_*z.dot(*Dz);

    return val;
  }

  // evaluate the gradient with respect to z
  void gradient_z(HDSA::Vector<RealT> & grad, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true) 
  {
    HDSA::Ptr<HDSA::Vector<RealT> > m = z.Clone();
    for(int k = 0; k < z.dimension(); k++)
      {
	m->Replace_Element(k,std::pow(z(k),3.0));
      }
    m->axpy(-1.0,*data_);
    HDSA::Ptr<HDSA::Vector<RealT> > dm = z.Clone();
    for(int k = 0; k < z.dimension(); k++)
      {
	dm->Replace_Element(k,(*m)(k)*3.0*std::pow(z(k),2.0));
      }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Wd_matvec = W_->Multiply(*dm,false);
    HDSA::Ptr<HDSA::Vector<RealT> > Wdm = z.Clone();
    Wd_matvec->Write_Column_to_Vector(0,grad);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > D_matvec = D_->Multiply(z,false);
    HDSA::Ptr<HDSA::Vector<RealT> > Dz = z.Clone();
    D_matvec->Write_Column_to_Vector(0,*Dz);
    grad.axpy(beta_,*Dz);
  } 

  // evaluate the z,z hessian vector product
  void hessVec_z_z(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta,
		   const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr)
  {
  HDSA::Ptr<HDSA::Vector<RealT> > m = z.Clone();
    for(int k = 0; k < z.dimension(); k++)
      {
	m->Replace_Element(k,std::pow(z(k),3.0));
      }
    m->axpy(-1.0,*data_);
    HDSA::Ptr<HDSA::Vector<RealT> > dm = z.Clone();
    for(int k = 0; k < z.dimension(); k++)
      {
	dm->Replace_Element(k,((*m)(k)*6.0*z(k)+9.0*std::pow(z(k),4.0))*v(k));
      }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Wd_matvec = W_->Multiply(*dm,false);
    HDSA::Ptr<HDSA::Vector<RealT> > Wdm = z.Clone();
    Wd_matvec->Write_Column_to_Vector(0,hv);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > D_matvec = D_->Multiply(v,false);
    HDSA::Ptr<HDSA::Vector<RealT> > Dv = z.Clone();
    D_matvec->Write_Column_to_Vector(0,*Dv);
    hv.axpy(beta_,*Dv);
  }

};


#endif
