#ifndef HDSA_MD_NUMERIC_LAPLACIAN_Z_PRIOR_INTERFACE_HPP
#define HDSA_MD_NUMERIC_LAPLACIAN_Z_PRIOR_INTERFACE_HPP

#include "HDSA_MD_Elliptic_z_Prior_Interface.hpp"
#include "HDSA_MD_Data_Interface.hpp"
#include "HDSA_Sparse_Matrix.hpp"
#include "HDSA_Sparse_Matrix_Solver.hpp"
#include "HDSA_MD_u_Prior_Interface.hpp"
#include "HDSA_MD_z_Hyperparameter_Interface.hpp"
#include "HDSA_MD_Determine_z_Hyperparameters_Decl.hpp"
#include "HDSA_Operator_Sqrt.hpp"

namespace HDSA
{

  template <class RealT>
  class MD_Numeric_Laplacian_z_Prior_Interface : public HDSA::MD_Elliptic_z_Prior_Interface<RealT>
  {

  private:
    const HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> S_;
    const HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> M_;
    const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> data_interface_;
    const HDSA::Ptr<HDSA::MD_z_Hyperparameter_Interface<RealT>> z_hyperparam_interface_;
    const HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT>> u_prior_interface_;
    HDSA::Ptr<HDSA::MD_Determine_z_Hyperparameters<RealT>> determine_z_hyperparams_;
    RealT beta_z_;
    HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> E_z_;
    HDSA::Ptr<HDSA::Sparse_Matrix_Solver<RealT>> E_z_solver_;
    HDSA::Ptr<HDSA::Sparse_Matrix_Solver<RealT>> M_z_solver_;

  public:
    void Apply_E_z_Inverse(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in) const
    {
      E_z_solver_->Apply_A_Inverse(z_out, z_in);
    }

    void Apply_E_z_Inverse_Transpose(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in) const
    {
      E_z_solver_->Apply_A_Inverse(z_out, z_in);
    }

    void Apply_M_z(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in) const
    {
      M_->Apply(z_out, z_in);
    }

    virtual void Sample_with_Covariance_W_z_Acute_Inverse(HDSA::MultiVector<RealT> &samples) const
    {
      HDSA::Ptr<M_z_Sqrt<RealT>> M_sqrt = HDSA::makePtr<M_z_Sqrt<RealT>>(M_);
      samples.zeros();
      for (int k = 0; k < samples.Number_of_Vectors(); k++)
      {
        HDSA::Ptr<HDSA::Vector<RealT>> omega = samples[k]->clone();
        omega->randomize_standard_normal();
        HDSA::Ptr<HDSA::Vector<RealT>> vec = samples[k]->clone();
        M_sqrt->Apply_Operator_Sqrt(*vec, *omega);
        Apply_E_z_Inverse(*samples[k], *vec);
      }
    }

    virtual void Apply_E_z(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in) const
    {
      E_z_->Apply(z_out, z_in);
    }

    virtual void Apply_E_z_Transpose(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in) const
    {
      E_z_->Apply(z_out, z_in);
    }

    virtual void Apply_M_z_Inverse(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in) const
    {
      M_z_solver_->Apply_A_Inverse(z_out, z_in);
    }

    MD_Numeric_Laplacian_z_Prior_Interface(const HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> &S, const HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> &M, const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> &data_interface, const HDSA::Ptr<HDSA::MD_z_Hyperparameter_Interface<RealT>> &z_hyperparam_interface, const HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT>> &u_prior_interface) : HDSA::MD_Elliptic_z_Prior_Interface<RealT>(z_hyperparam_interface->Get_alpha_z()), S_(S), M_(M), data_interface_(data_interface), z_hyperparam_interface_(z_hyperparam_interface), u_prior_interface_(u_prior_interface)
    {
      E_z_ = M_->clone();
      determine_z_hyperparams_ = HDSA::makePtr<HDSA::MD_Determine_z_Hyperparameters<RealT>>(data_interface_, z_hyperparam_interface_, u_prior_interface_);

      if (z_hyperparam_interface_->Get_beta_z() == 0.0)
      {
        HDSA_TEST_FOR_EXCEPTION(true, std::logic_error,
                                "Error in HDSA::MD_Numeric_Laplacian_z_Prior_Interface: The value of beta_z must be specificed" << std::endl);
      }
      Set_beta_z(z_hyperparam_interface_->Get_beta_z());

      if (z_hyperparam_interface_->Get_alpha_z() == 0.0)
      {
        determine_z_hyperparams_->Determine_alpha_z(this);
      }
      this->Set_alpha_z(z_hyperparam_interface_->Get_alpha_z());

      M_z_solver_ = HDSA::makePtr<HDSA::Sparse_Matrix_Solver<RealT>>(M_);
    }

    virtual ~MD_Numeric_Laplacian_z_Prior_Interface()
    {
    }

    void Set_beta_z(RealT beta_z_new)
    {
      E_z_->set(*M_);
      E_z_->axpy(beta_z_new, *S_);
      beta_z_ = beta_z_new;
      E_z_solver_ = HDSA::makePtr<HDSA::Sparse_Matrix_Solver<RealT>>(E_z_);
    }

    template <class ScalarType>
    class M_z_Sqrt : public HDSA::Operator_Sqrt<ScalarType>
    {
    private:
      const HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> M_;

    public:
      M_z_Sqrt(const HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> &M) : M_(M)
      {
      }

      ~M_z_Sqrt()
      {
      }

      void Apply_Operator(HDSA::Vector<ScalarType> &vec_out, const HDSA::Vector<ScalarType> &vec_in) const
      {
        M_->Apply(vec_out, vec_in);
      }
    };
  };

}

#include "HDSA_MD_Determine_z_Hyperparameters_Def.hpp"

#endif
