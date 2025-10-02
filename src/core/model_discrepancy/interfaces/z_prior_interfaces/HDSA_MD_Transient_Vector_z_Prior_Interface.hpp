#ifndef HDSA_MD_TRANSIENT_VECTOR_Z_PRIOR_INTERFACE_HPP
#define HDSA_MD_TRANSIENT_VECTOR_Z_PRIOR_INTERFACE_HPP

#include "HDSA_MD_Elliptic_z_Prior_Interface.hpp"
#include "HDSA_MD_z_Hyperparameter_Interface.hpp"
#include "HDSA_MD_Determine_z_Hyperparameters_Decl.hpp"
#include "HDSA_Std_Vector.hpp"

namespace HDSA
{

  template <class RealT>
  class MD_Transient_Vector_z_Prior_Interface : public HDSA::MD_Elliptic_z_Prior_Interface<RealT>
  {

  private:
    const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> data_interface_;
    const HDSA::Ptr<HDSA::MD_z_Hyperparameter_Interface<RealT>> z_hyperparam_interface_;
    const HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT>> u_prior_interface_;
    int n_t_;
    RealT T_;
    int num_controls_;
    HDSA::Ptr<HDSA::MD_Determine_z_Hyperparameters<RealT>> determine_z_hyperparams_;
    RealT beta_t_;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> S_t_;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> M_t_;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> E_t_;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> V_;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> Lambda_;

  public:
    void Apply_E_z_Inverse(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in) const override
    {
      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> b = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(n_t_, num_controls_);
      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> x = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(n_t_, num_controls_);

      const HDSA::Transient_Vector<RealT> z_in_trans = dynamic_cast<const HDSA::Transient_Vector<RealT> &>(z_in);
      HDSA::Transient_Vector<RealT> z_out_trans = dynamic_cast<HDSA::Transient_Vector<RealT> &>(z_out);

      for (int i = 0; i < n_t_; i++)
      {
        const HDSA::Std_Vector<RealT> z_in_i_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(*z_in_trans[i]);
        for (int j = 0; j < num_controls_; j++)
        {
          b->Set_Entry(i, j, z_in_i_std(j));
        }
      }

      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> tmp = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(n_t_, num_controls_);
      V_->Multiply(*tmp, *b, true);
      for (int i = 0; i < n_t_; i++)
      {
        for (int j = 0; j < num_controls_; j++)
        {
          RealT val = (*tmp)(i, j) / std::sqrt((*Lambda_)(i, 0));
          tmp->Set_Entry(i, j, val);
        }
      }

      V_->Multiply(*x, *tmp);

      for (int i = 0; i < n_t_; i++)
      {
        HDSA::Std_Vector<RealT> z_out_i_std = dynamic_cast<HDSA::Std_Vector<RealT> &>(*z_out_trans[i]);
        for (int j = 0; j < num_controls_; j++)
        {
          z_out_i_std.Set_Entry(j, (*x)(i, j));
        }
      }
    }

    void Apply_E_z_Inverse_Transpose(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in) const override
    {
      Apply_E_z_Inverse(z_out, z_in);
    }

    void Apply_M_z(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in) const override
    {
      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> b = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(n_t_, num_controls_);
      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> x = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(n_t_, num_controls_);

      const HDSA::Transient_Vector<RealT> z_in_trans = dynamic_cast<const HDSA::Transient_Vector<RealT> &>(z_in);
      HDSA::Transient_Vector<RealT> z_out_trans = dynamic_cast<HDSA::Transient_Vector<RealT> &>(z_out);

      for (int i = 0; i < n_t_; i++)
      {
        const HDSA::Std_Vector<RealT> z_in_i_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(*z_in_trans[i]);
        for (int j = 0; j < num_controls_; j++)
        {
          b->Set_Entry(i, j, z_in_i_std(j));
        }
      }

      M_t_->Multiply(*x, *b, true);

      for (int i = 0; i < n_t_; i++)
      {
        HDSA::Std_Vector<RealT> z_out_i_std = dynamic_cast<HDSA::Std_Vector<RealT> &>(*z_out_trans[i]);
        for (int j = 0; j < num_controls_; j++)
        {
          z_out_i_std.Set_Entry(j, (*x)(i, j));
        }
      }
    }

    void Sample_with_Covariance_W_z_Acute_Inverse(HDSA::MultiVector<RealT> &samples) const override
    {
      int num_samps = samples.Number_of_Vectors();
      HDSA::Ptr<HDSA::Vector<RealT>> omega = samples[0]->Clone();
      for (int k = 0; k < num_samps; k++)
      {
        HDSA::Ptr<HDSA::Vector<RealT>> vec = samples[k];
        omega->Randomize_Standard_Normal();

        HDSA::Ptr<HDSA::Dense_Matrix<RealT>> b = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(n_t_, num_controls_);
        HDSA::Ptr<HDSA::Dense_Matrix<RealT>> x = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(n_t_, num_controls_);

        HDSA::Transient_Vector<RealT> omega_trans = dynamic_cast<HDSA::Transient_Vector<RealT> &>(*omega);
        HDSA::Transient_Vector<RealT> vec_trans = dynamic_cast<HDSA::Transient_Vector<RealT> &>(*vec);

        for (int i = 0; i < n_t_; i++)
        {
          HDSA::Std_Vector<RealT> omega_i_std = dynamic_cast<HDSA::Std_Vector<RealT> &>(*omega_trans[i]);
          for (int j = 0; j < num_controls_; j++)
          {
            b->Set_Entry(i, j, std::sqrt((*Lambda_)(i, 0)) * omega_i_std(j));
          }
        }

        V_->Multiply(*x, *b);

        for (int i = 0; i < n_t_; i++)
        {
          HDSA::Std_Vector<RealT> vec_i_std = dynamic_cast<HDSA::Std_Vector<RealT> &>(*vec_trans[i]);
          for (int j = 0; j < num_controls_; j++)
          {
            vec_i_std.Set_Entry(j, (*x)(i, j));
          }
        }
      }
    }

    void Apply_E_z(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in) const override
    {
      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> b = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(n_t_, num_controls_);
      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> x = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(n_t_, num_controls_);

      const HDSA::Transient_Vector<RealT> z_in_trans = dynamic_cast<const HDSA::Transient_Vector<RealT> &>(z_in);
      HDSA::Transient_Vector<RealT> z_out_trans = dynamic_cast<HDSA::Transient_Vector<RealT> &>(z_out);

      for (int i = 0; i < n_t_; i++)
      {
        const HDSA::Std_Vector<RealT> z_in_i_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(*z_in_trans[i]);
        for (int j = 0; j < num_controls_; j++)
        {
          b->Set_Entry(i, j, z_in_i_std(j));
        }
      }

      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> tmp = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(n_t_, num_controls_);
      V_->Multiply(*tmp, *b, true);
      for (int i = 0; i < n_t_; i++)
      {
        for (int j = 0; j < num_controls_; j++)
        {
          RealT val = (*tmp)(i, j) * std::sqrt((*Lambda_)(i, 0));
          tmp->Set_Entry(i, j, val);
        }
      }

      V_->Multiply(*x, *tmp);

      for (int i = 0; i < n_t_; i++)
      {
        HDSA::Std_Vector<RealT> z_out_i_std = dynamic_cast<HDSA::Std_Vector<RealT> &>(*z_out_trans[i]);
        for (int j = 0; j < num_controls_; j++)
        {
          z_out_i_std.Set_Entry(j, (*x)(i, j));
        }
      }
    }

    void Apply_E_z_Transpose(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in) const override
    {
      Apply_E_z(z_out, z_in);
    }

    void Apply_M_z_Inverse(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in) const override
    {
      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> b = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(n_t_, num_controls_);
      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> x = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(n_t_, num_controls_);

      const HDSA::Transient_Vector<RealT> z_in_trans = dynamic_cast<const HDSA::Transient_Vector<RealT> &>(z_in);
      HDSA::Transient_Vector<RealT> z_out_trans = dynamic_cast<HDSA::Transient_Vector<RealT> &>(z_out);

      for (int i = 0; i < n_t_; i++)
      {
        const HDSA::Std_Vector<RealT> z_in_i_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(*z_in_trans[i]);
        for (int j = 0; j < num_controls_; j++)
        {
          b->Set_Entry(i, j, z_in_i_std(j));
        }
      }

      M_t_->Multiply(*x, *b, true);
      HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*M_t_, *x, *b);

      for (int i = 0; i < n_t_; i++)
      {
        HDSA::Std_Vector<RealT> z_out_i_std = dynamic_cast<HDSA::Std_Vector<RealT> &>(*z_out_trans[i]);
        for (int j = 0; j < num_controls_; j++)
        {
          z_out_i_std.Set_Entry(j, (*x)(i, j));
        }
      }
    }

    MD_Transient_Vector_z_Prior_Interface(const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> &data_interface, const HDSA::Ptr<HDSA::MD_z_Hyperparameter_Interface<RealT>> &z_hyperparam_interface, const HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT>> &u_prior_interface, int n_t, RealT T, int num_controls) : HDSA::MD_Elliptic_z_Prior_Interface<RealT>(z_hyperparam_interface->Get_alpha_z()), data_interface_(data_interface), z_hyperparam_interface_(z_hyperparam_interface), u_prior_interface_(u_prior_interface), n_t_(n_t), T_(T), num_controls_(num_controls)
    {
      RealT h = 1.0 / static_cast<RealT>(n_t_ - 1);
      S_t_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(n_t_, n_t_);
      M_t_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(n_t_, n_t_);
      E_t_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(n_t_, n_t_);

      S_t_->Set_Entry(0, 0, 1.0 / h);
      S_t_->Set_Entry(0, 1, -1.0 / h);
      for (int i = 1; i < n_t_ - 1; i++)
      {
        S_t_->Set_Entry(i, i, 2.0 / h);
        S_t_->Set_Entry(i, i - 1, -1.0 / h);
        S_t_->Set_Entry(i, i + 1, -1.0 / h);
      }
      S_t_->Set_Entry(n_t_ - 1, n_t_ - 2, -1.0 / h);
      S_t_->Set_Entry(n_t_ - 1, n_t_ - 1, 1.0 / h);

      M_t_->Set_Entry(0, 0, (1.0 / 3.0) * h);
      M_t_->Set_Entry(0, 1, (1.0 / 6.0) * h);
      for (int i = 1; i < n_t_ - 1; i++)
      {
        M_t_->Set_Entry(i, i, (2.0 / 3.0) * h);
        M_t_->Set_Entry(i, i - 1, (1.0 / 6.0) * h);
        M_t_->Set_Entry(i, i + 1, (1.0 / 6.0) * h);
      }
      M_t_->Set_Entry(n_t_ - 1, n_t_ - 2, (1.0 / 6.0) * h);
      M_t_->Set_Entry(n_t_ - 1, n_t_ - 1, (1.0 / 3.0) * h);

      determine_z_hyperparams_ = HDSA::makePtr<HDSA::MD_Determine_z_Hyperparameters<RealT>>(data_interface_, z_hyperparam_interface_, u_prior_interface_);

      if (z_hyperparam_interface_->Get_beta_t() == 0.0)
      {
        HDSA_TEST_FOR_EXCEPTION(true, std::logic_error,
                                "Error in HDSA::MD_Transient_Vector_z_Prior_Interface: The value of beta_t must be specificed" << std::endl);
      }
      Set_beta_t(z_hyperparam_interface_->Get_beta_t());

      if (z_hyperparam_interface_->Get_alpha_z() == 0.0)
      {
        determine_z_hyperparams_->Determine_alpha_z(this);
      }
      this->Set_alpha_z(z_hyperparam_interface_->Get_alpha_z());
    }

    virtual ~MD_Transient_Vector_z_Prior_Interface()
    {
    }

    void Set_beta_t(RealT beta_t_new)
    {
      beta_t_ = beta_t_new;
      for (int i = 0; i < n_t_; i++)
      {
        for (int j = 0; j < n_t_; j++)
        {
          RealT val = beta_t_ * (*S_t_)(i, j) + (*M_t_)(i, j);
          E_t_->Set_Entry(i, j, val);
        }
      }

      V_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(n_t_, n_t_);
      Lambda_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(n_t_, 1);
      HDSA::Linear_Algebra::Symmetric_Gen_Eig_Decomposition<RealT>(*E_t_, *M_t_, *V_, *Lambda_);
    }

    int Get_n_t(void) const
    {
      return n_t_;
    }

    RealT Get_T(void) const
    {
      return T_;
    }
  };

}

#include "HDSA_MD_Determine_z_Hyperparameters_Def.hpp"

#endif
