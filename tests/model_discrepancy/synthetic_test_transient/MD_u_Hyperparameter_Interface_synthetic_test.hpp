#ifndef HDSA_MD_U_HYPERPARAMETER_INTERFACE_SYNTHETIC_TEST_HPP
#define HDSA_MD_U_HYPERPARAMETER_INTERFACE_SYNTHETIC_TEST_HPP

template <class RealT>
class MD_u_Hyperparameter_Interface_synthetic_test : public HDSA::MD_u_Hyperparameter_Interface<RealT>
{

private:
  int n_y_;
  int n_t_;

public:
  HDSA::Ptr<HDSA::Vector<RealT>> Load_Spatial_Node_Data(void) const
  {
    HDSA::Ptr<Std_Vector<RealT>> x = HDSA::makePtr<Std_Vector<RealT>>(n_y_);
    for (int k = 0; k < n_y_; k++)
    {
      x->Replace_Element(k, static_cast<RealT>(k) / static_cast<RealT>(n_y_ - 1));
    }
    return x;
  }

  std::vector<RealT> Load_Time_Node_Data(void) const
  {
    std::vector<RealT> vec = std::vector<RealT>(n_t_);
    for (int k = 0; k < n_t_; k++)
    {
      vec[k] = static_cast<RealT>(k) / static_cast<RealT>(n_t_ - 1);
    }
    return vec;
  }

  MD_u_Hyperparameter_Interface_synthetic_test(int n_y, int n_t) : HDSA::MD_u_Hyperparameter_Interface<RealT>(true)
  {
    n_y_ = n_y;
    n_t_ = n_t;
  }

  virtual ~MD_u_Hyperparameter_Interface_synthetic_test()
  {
  }
};

#endif
