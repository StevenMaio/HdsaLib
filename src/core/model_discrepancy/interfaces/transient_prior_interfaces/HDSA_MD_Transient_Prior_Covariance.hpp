#ifndef HDSA_MD_TRANSIENT_PRIOR_COVARIANCE_HPP
#define HDSA_MD_TRANSIENT_PRIOR_COVARIANCE_HPP

#include "HDSA_MD_Data_Interface.hpp"
#include "HDSA_MD_u_Hyperparameter_Interface.hpp"
#include "HDSA_MD_Determine_u_Hyperparameters_Def.hpp"

namespace HDSA
{

  template <class RealT>
  class MD_Transient_Prior_Covariance
  {

  private:
    const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> data_interface_;
    const HDSA::Ptr<HDSA::MD_u_Hyperparameter_Interface<RealT>> u_hyperparam_interface_;
    HDSA::Ptr<HDSA::MD_Determine_u_Hyperparameters<RealT>> determine_u_hyperparams_;
    RealT T_;
    int n_t_;
    int n_y_;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> S_t_;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> M_t_;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> E_t_;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> W_t_inv_;
    std::vector<RealT> alpha_t_;
    RealT beta_t_;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> evecs_;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> evals_;

  public:
    MD_Transient_Prior_Covariance(const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> &data_interface, const HDSA::Ptr<HDSA::MD_u_Hyperparameter_Interface<RealT>> &u_hyperparam_interface, RealT T, int n_t, int n_y) : data_interface_(data_interface), u_hyperparam_interface_(u_hyperparam_interface), T_(T), n_t_(n_t), n_y_(n_y)
    {
      RealT h = 1.0 / static_cast<RealT>(n_t_ - 1);
      S_t_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(n_t_, n_t_);
      M_t_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(n_t_, n_t_);
      E_t_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(n_t_, n_t_);

      S_t_->Replace_Element(0, 0, 1.0 / h);
      S_t_->Replace_Element(0, 1, -1.0 / h);
      for (int i = 1; i < n_t_ - 1; i++)
      {
        S_t_->Replace_Element(i, i, 2.0 / h);
        S_t_->Replace_Element(i, i - 1, -1.0 / h);
        S_t_->Replace_Element(i, i + 1, -1.0 / h);
      }
      S_t_->Replace_Element(n_t_ - 1, n_t_ - 2, -1.0 / h);
      S_t_->Replace_Element(n_t_ - 1, n_t_ - 1, 1.0 / h);

      M_t_->Replace_Element(0, 0, (1.0 / 3.0) * h);
      M_t_->Replace_Element(0, 1, (1.0 / 6.0) * h);
      for (int i = 1; i < n_t_ - 1; i++)
      {
        M_t_->Replace_Element(i, i, (2.0 / 3.0) * h);
        M_t_->Replace_Element(i, i - 1, (1.0 / 6.0) * h);
        M_t_->Replace_Element(i, i + 1, (1.0 / 6.0) * h);
      }
      M_t_->Replace_Element(n_t_ - 1, n_t_ - 2, (1.0 / 6.0) * h);
      M_t_->Replace_Element(n_t_ - 1, n_t_ - 1, (1.0 / 3.0) * h);

      determine_u_hyperparams_ = HDSA::makePtr<HDSA::MD_Determine_u_Hyperparameters<RealT>>(data_interface_, u_hyperparam_interface_);
      if (u_hyperparam_interface_->Get_beta_t() == 0.0)
      {
        std::cout << "Error: the value of beta_t must be specificed" << std::endl;
      }
      Set_beta_t(u_hyperparam_interface_->Get_beta_t());
    }

    virtual ~MD_Transient_Prior_Covariance()
    {
    }

    HDSA::Ptr<HDSA::MD_u_Hyperparameter_Interface<RealT>> Get_Hyperparameter_Interface(void) const
    {
      return u_hyperparam_interface_;
    }

    HDSA::Ptr<HDSA::MD_Determine_u_Hyperparameters<RealT>> Get_Determine_Hyperparameters(void) const
    {
      return determine_u_hyperparams_;
    }

    int Get_Num_Time_Nodes(void) const
    {
      return n_t_;
    }

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> Get_M_t(void) const
    {
      return M_t_;
    }


    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> Get_W_t_Inverse(void) const
    {
      return W_t_inv_;
    }

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> Get_Evecs(void) const
    {
      return evecs_;
    }

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> Get_Evals(void) const
    {
      return evals_;
    }

    void Set_alpha_t(std::vector<RealT> alpha_t_new)
    {
      alpha_t_ = alpha_t_new;
      Compute_Time_Covariance_GEVP();
    }

    void Set_beta_t(RealT beta_t_new)
    {
      beta_t_ = beta_t_new;
      for (int i = 0; i < n_t_; i++)
      {
        for (int j = 0; j < n_t_; j++)
        {
          RealT val = beta_t_ * (*S_t_)(i, j) + (*M_t_)(i, j);
          E_t_->Replace_Element(i, j, val);
        }
      }
    }

    void Compute_Time_Covariance_GEVP(void)
    {
      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> A;
      if (alpha_t_.size() == 1)
      {
        A = E_t_;
      }
      else
      {
        A = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(n_t_, n_t_);
        for (int i = 0; i < n_t_; i++)
        {
          for (int j = 0; j < n_t_; j++)
          {
            RealT val = (*E_t_)(i, j) / std::sqrt(alpha_t_[i] * alpha_t_[j]);
            A->Replace_Element(i, j, val);
          }
        }
      }
      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> V = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(n_t_, n_t_);
      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> S = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(n_t_, 1);
      HDSA::Linear_Algebra::Symmetric_Gen_Eig_Decomposition<RealT>(*A, *M_t_, *V, *S);
      evecs_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(n_t_, n_t_);
      evals_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(n_t_, 1);
      for (int j = 0; j < n_t_; j++)
      {
        evals_->Replace_Element(j, 0, 1.0 / (*S)(n_t_ - 1 - j, 0));
        for (int i = 0; i < n_t_; i++)
        {
          evecs_->Replace_Element(i, j, (*V)(i, n_t_ - 1 - j));
        }
      }
      W_t_inv_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(n_t_, n_t_);
      for (int i = 0; i < n_t_; i++)
      {
        for (int j = 0; j < n_t_; j++)
        {
          RealT val = 0.0;
          for (int k = 0; k < n_t_; k++)
          {
            val += (*evecs_)(i,k) * (*evecs_)(j,k) * (*evals_)(k,0);
          }
          W_t_inv_->Replace_Element(i,j,val);
        }
      }
    }
  };

}

#endif
