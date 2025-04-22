#ifndef HDSA_MD_NUMERIC_LAPLACIAN_U_PRIOR_INTERFACE_HPP
#define HDSA_MD_NUMERIC_LAPLACIAN_U_PRIOR_INTERFACE_HPP

namespace HDSA
{

  template <class RealT>
  class MD_Numeric_Laplacian_u_Prior_Interface
  {

  private:
    const HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> S_;
    const HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> M_;
    const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> data_interface_;
    const HDSA::Ptr<HDSA::MD_u_Hyperparameter_Interface<RealT>> u_hyperparam_interface_;
    HDSA::Ptr<HDSA::MD_Determine_u_Hyperparameters<RealT> > & determine_u_hyperparams_;
    RealT beta_u_;
    HDSA::Ptr<HDSA::Sparse_Matrix<RealT> > E_u_;

  public:
    MD_Numeric_Laplacian_u_Prior_Interface(const HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> &S, const HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> &M, const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> &data_interface, const HDSA::Ptr<HDSA::MD_u_Hyperparameter_Interface<RealT>> &u_hyperparam_interface) : S_(S), M_(M), data_interface_(data_interface), u_hyperparam_interface_(u_hyperparam_interface)
    {
      E_u_ = M->clone();
    }

    virtual ~MD_Numeric_Laplacian_u_Prior_Interface()
    {
    }
  };

}

#endif
