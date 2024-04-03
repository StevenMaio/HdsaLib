#ifndef HDSA_MD_DATA_INTERFACE_MODEL_DISCREPANCY_SYNTHETIC_TEST_HPP
#define HDSA_MD_DATA_INTERFACE_MODEL_DISCREPANCY_SYNTHETIC_TEST_HPP

template <class RealT>
class MD_Data_Interface_model_discrepancy_synthetic_test : public HDSA::MD_Data_Interface<RealT> {

private:
  int m_; // Mesh resolution
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > x_; // Mesh nodes on [0,1]

public:
  MD_Data_Interface_model_discrepancy_synthetic_test()
  {  
    m_ = 51;
    RealT h = 1.0/static_cast<RealT>(m_-1);
    x_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    for(int k = 0; k < m_; k++)
      {
	x_->Replace_Element(k,0,static_cast<RealT>(k)/static_cast<RealT>(m_-1));
      }
  }

  virtual ~MD_Data_Interface_model_discrepancy_synthetic_test()
  { }

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

  HDSA::Ptr<HDSA::MultiVector<RealT> > Load_D_Data( ) const 
  {
    HDSA::Ptr<Std_Vector<RealT> > d = HDSA::makePtr<Std_Vector<RealT> >(m_);
    HDSA::Ptr<HDSA::MultiVector<RealT> > D = HDSA::makePtr<HDSA::MultiVector<RealT> >(2,*d);

    HDSA::Ptr<HDSA::Vector<RealT> > d0 = (*D)[0];
    HDSA::Ptr<HDSA::Vector<RealT> > d1 = (*D)[1];
    Std_Vector<RealT> d0_std = dynamic_cast<Std_Vector<RealT>&>(*d0);
    Std_Vector<RealT> d1_std = dynamic_cast<Std_Vector<RealT>&>(*d1);

    for(int k = 0; k < m_; k++)
      {
	d0_std.Replace_Element(k,0.2*std::pow((*x_)(k,0)+1.0,2.0));
	d1_std.Replace_Element(k,0.2*std::pow((*x_)(k,0) + std::pow((*x_)(k,0),2.0),2.0));
      }
    return D;
  }

};

#endif


