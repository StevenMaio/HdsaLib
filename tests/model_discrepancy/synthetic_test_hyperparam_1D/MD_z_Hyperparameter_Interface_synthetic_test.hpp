#ifndef HDSA_MD_Z_HYPERPARAMETER_INTERFACE_SYNTHETIC_TEST_HPP
#define HDSA_MD_Z_HYPERPARAMETER_INTERFACE_SYNTHETIC_TEST_HPP


  template <class RealT>
  class MD_z_Hyperparameter_Interface_synthetic_test : public HDSA::MD_z_Hyperparameter_Interface<RealT>
  {

  private:
    int m_;

  public:
    HDSA::Ptr<HDSA::Vector<RealT>> Load_Spatial_Node_Data(void) const
    {
      HDSA::Ptr<Std_Vector<RealT>> x = HDSA::makePtr<Std_Vector<RealT>>(m_);
      for (int k = 0; k < m_; k++)
      {
        x->Replace_Element(k, static_cast<RealT>(k)/static_cast<RealT>(m_-1));
      }
      return x;
    }

    MD_z_Hyperparameter_Interface_synthetic_test() : HDSA::MD_z_Hyperparameter_Interface<RealT>("spatial field")
    {
      m_ = 51;
    }

    virtual ~MD_z_Hyperparameter_Interface_synthetic_test()
    {
    }
  };


#endif