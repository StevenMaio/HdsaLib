#ifndef HDSA_MD_LUMPED_MASS_U_PRIOR_INTERFACE_HPP
#define HDSA_MD_LUMPED_MASS_U_PRIOR_INTERFACE_HPP

#include "HDSA_MD_Scaled_u_Prior_Interface.hpp"
#include "HDSA_MD_Data_Interface.hpp"
#include "HDSA_Linear_Operator.hpp"
#include "HDSA_Linear_Algebra.hpp"
#include "HDSA_Sparse_Matrix.hpp"
#include "HDSA_Sparse_Matrix_Solver.hpp"
#include "HDSA_Sparse_Matrix_Sqrt.hpp"
#include "HDSA_MD_u_Hyperparameter_Interface.hpp"
#include "HDSA_MD_Determine_u_Hyperparameters_Decl.hpp"

namespace HDSA
{

  template <class RealT>
  class MD_Lumped_Mass_u_Prior_Interface : public HDSA::MD_Scaled_u_Prior_Interface<RealT>
  {

  private:
    const HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> S_;
    const HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> M_;
    const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> data_interface_;
    const HDSA::Ptr<HDSA::MD_u_Hyperparameter_Interface<RealT>> u_hyperparam_interface_;
    const HDSA::Ptr<const HDSA::Comm<int>> comm_;
    const HDSA::Ptr<HDSA::Random_Number_Generator<RealT>> random_number_generator_;
    int verbosity_;
    HDSA::Ptr<HDSA::MD_Determine_u_Hyperparameters<RealT>> determine_u_hyperparams_;
    RealT beta_u_;
    HDSA::Ptr<HDSA::Vector<RealT>> M_lumped_;
    HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> E_u_;
    HDSA::Ptr<HDSA::Sparse_Matrix_Solver<RealT>> E_u_solver_;
    HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> W_u_acute_;
    std::vector<RealT> scalars_;
    std::vector<HDSA::Ptr<HDSA::Sparse_Matrix<RealT>>> W_u_acute_plus_scalar_M_u_;
    std::vector<HDSA::Ptr<HDSA::Sparse_Matrix_Solver<RealT>>> W_u_acute_plus_scalar_M_u_solver_;
    std::vector<HDSA::Ptr<HDSA::Sparse_Matrix_Sqrt<RealT>>> W_u_acute_plus_scalar_M_u_sqrt_;

  public:
    void Apply_M_u(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &u_in) const override
    {
      M_->Apply(u_out, u_in);
    }

    void Apply_W_u_Acute_Inverse(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &u_in) const override
    {
      HDSA::Ptr<HDSA::Vector<RealT>> u_tmp = u_in.Clone();
      Apply_E_u_Inverse(*u_tmp, u_in);
      for (int k = 0; k < u_in.Dimension(); k++)
      {
        RealT val = u_tmp->Get_Entry(k) * M_lumped_->Get_Entry(k);
        u_tmp->Set_Entry(k, val);
      }
      Apply_E_u_Inverse(u_out, *u_tmp);
    }

    void Apply_W_u_Acute_Plus_scalar_M_u_Inverse(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &u_in, const RealT &scalar) const override
    {
      int i = -1;
      for (int k = 0; k < scalars_.size(); k++)
      {
        if (std::abs(scalars_[k] - scalar) < 1.e-12)
        {
          i = k;
          break;
        }
      }

      if (i < 0)
      {
        HDSA_TEST_FOR_EXCEPTION(true, std::logic_error,
                                "Error in HDSA::MD_Lumped_Mass_u_Prior_Interface::Apply_W_u_Acute_Plus_scalar_M_u_Inverse: scalar has not been set by Precompute_W_u_Plus_scalar_M_u_Data" << std::endl);
      }

      W_u_acute_plus_scalar_M_u_solver_[i]->Apply_A_Inverse(u_out, u_in);
    }

    void Sample_with_Covariance_W_u_Acute_Inverse(HDSA::MultiVector<RealT> &samples) const override
    {
      int num_samples = samples.Number_of_Vectors();
      HDSA::Ptr<HDSA::Vector<RealT>> u_tmp = samples[0]->Clone();
      for (int k = 0; k < num_samples; k++)
      {
        u_tmp->Randomize_Standard_Normal();
        for (int k = 0; k < u_tmp->Dimension(); k++)
        {
          RealT val = u_tmp->Get_Entry(k) * std::sqrt(M_lumped_->Get_Entry(k));
          u_tmp->Set_Entry(k, val);
        }
        Apply_E_u_Inverse(*samples[k], *u_tmp);
      }
    }

    void Sample_with_Covariance_W_u_Acute_Plus_scalar_M_u_Inverse(HDSA::MultiVector<RealT> &samples, const RealT &scalar) const override
    {
      int i = -1;
      for (int k = 0; k < scalars_.size(); k++)
      {
        if (std::abs(scalars_[k] - scalar) < 1.e-12)
        {
          i = k;
          break;
        }
      }

      if (i < 0)
      {
        HDSA_TEST_FOR_EXCEPTION(true, std::logic_error,
                                "Error in HDSA::MD_Lumped_Mass_u_Prior_Interface::Sample_with_Covariance_W_u_Acute_Plus_scalar_M_u_Inverse: scalar has not been set by Precompute_W_u_Plus_scalar_M_u_Data" << std::endl);
      }

      HDSA::Ptr<HDSA::Vector<RealT>> vec_rand = samples[0]->Clone();
      HDSA::Ptr<HDSA::Vector<RealT>> vec = samples[0]->Clone();
      for (int s = 0; s < samples.Number_of_Vectors(); s++)
      {
        vec_rand->Randomize_Standard_Normal();
        W_u_acute_plus_scalar_M_u_sqrt_[i]->Matrix_Sqrt_Apply(*vec, *vec_rand);
        Apply_W_u_Acute_Plus_scalar_M_u_Inverse(*samples[s], *vec, scalar);
      }
    }

    void Precompute_W_u_Plus_scalar_M_u_Data(RealT &scalar) override
    {
      RealT scalar_shift = scalar * u_hyperparam_interface_->Get_alpha_u();
      scalars_.push_back(scalar_shift);

      HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> A = W_u_acute_->Clone();
      A->Set(*W_u_acute_);
      A->Scaled_Plus(scalar_shift, *M_);
      A->Set_Symmetric();
      W_u_acute_plus_scalar_M_u_.push_back(A);

      HDSA::Ptr<HDSA::Incomplete_Chol_Factor<RealT>> L = HDSA::makePtr<HDSA::Incomplete_Chol_Factor<RealT>>(A);

      // HDSA::Ptr<HDSA::Sparse_Matrix_Solver<RealT>> A_solver = HDSA::makePtr<HDSA::Sparse_Matrix_Solver<RealT>>(A);
      HDSA::Ptr<HDSA::Sparse_Matrix_Solver<RealT>> A_solver = HDSA::makePtr<HDSA::Sparse_Matrix_Solver<RealT>>(A, false);
      A_solver->Set_Incomplete_Factor(L);
      W_u_acute_plus_scalar_M_u_solver_.push_back(A_solver);

      HDSA::Ptr<HDSA::Sparse_Matrix_Sqrt<RealT>> A_sqrt = HDSA::makePtr<HDSA::Sparse_Matrix_Sqrt<RealT>>(A);
      A_sqrt->Set_Incomplete_Factor(L);
      W_u_acute_plus_scalar_M_u_sqrt_.push_back(A_sqrt);
    }

    MD_Lumped_Mass_u_Prior_Interface(const HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> &S, const HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> &M, const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> &data_interface, const HDSA::Ptr<HDSA::MD_u_Hyperparameter_Interface<RealT>> &u_hyperparam_interface, const HDSA::Ptr<const HDSA::Comm<int>> &comm) : HDSA::MD_Scaled_u_Prior_Interface<RealT>(u_hyperparam_interface->Get_alpha_u()), S_(S), M_(M), data_interface_(data_interface), u_hyperparam_interface_(u_hyperparam_interface), comm_(comm), random_number_generator_(HDSA::makePtr<HDSA::Random_Number_Generator<RealT>>())
    {
      verbosity_ = 0;
      Auxillary_Constructor();
    }

    MD_Lumped_Mass_u_Prior_Interface(const HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> &S, const HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> &M, const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> &data_interface, const HDSA::Ptr<HDSA::MD_u_Hyperparameter_Interface<RealT>> &u_hyperparam_interface, const HDSA::Ptr<const HDSA::Comm<int>> &comm, int seed) : HDSA::MD_Scaled_u_Prior_Interface<RealT>(u_hyperparam_interface->Get_alpha_u()), S_(S), M_(M), data_interface_(data_interface), u_hyperparam_interface_(u_hyperparam_interface), comm_(comm), random_number_generator_(HDSA::makePtr<HDSA::Random_Number_Generator<RealT>>(seed))
    {
      verbosity_ = 0;
      Auxillary_Constructor();
    }

    MD_Lumped_Mass_u_Prior_Interface(const HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> &S, const HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> &M, const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> &data_interface, const HDSA::Ptr<HDSA::MD_u_Hyperparameter_Interface<RealT>> &u_hyperparam_interface, const HDSA::Ptr<const HDSA::Comm<int>> &comm, const HDSA::Ptr<HDSA::Random_Number_Generator<RealT>> &random_number_generator) : HDSA::MD_Scaled_u_Prior_Interface<RealT>(u_hyperparam_interface->Get_alpha_u()), S_(S), M_(M), data_interface_(data_interface), u_hyperparam_interface_(u_hyperparam_interface), comm_(comm), random_number_generator_(random_number_generator)
    {
      verbosity_ = 0;
      Auxillary_Constructor();
    }

    void Set_Verbosity(int verbosity)
    {
      verbosity_ = verbosity;
    }

    void Disable_Sampling_Preconditioner(void)
    {
      for (int k = 0; k < W_u_acute_plus_scalar_M_u_sqrt_.size(); k++)
      {
        W_u_acute_plus_scalar_M_u_sqrt_[k]->Disable_Incomplete_Factorization();
      }
    }

    void Auxillary_Constructor()
    {
      if (HDSA::Ptr<const HDSA::Transient_Vector<RealT>> u_opt_trans = HDSA::dynamicPtrCast<const HDSA::Transient_Vector<RealT>>(data_interface_->Get_u_opt()))
      {
        M_lumped_ = (*u_opt_trans)[0]->Clone();
      }
      else
      {
        M_lumped_ = data_interface_->Get_u_opt()->Clone();
      }
      HDSA::Ptr<HDSA::Vector<RealT>> u_tmp = M_lumped_->Clone();
      u_tmp->Set_Scalar(1.0);
      M_->Apply(*M_lumped_, *u_tmp);

      E_u_ = M_->Clone();
      determine_u_hyperparams_ = HDSA::makePtr<HDSA::MD_Determine_u_Hyperparameters<RealT>>(data_interface_, u_hyperparam_interface_);

      if (u_hyperparam_interface_->Get_beta_u() == 0.0)
      {
        HDSA_TEST_FOR_EXCEPTION(true, std::logic_error,
                                "Error in HDSA::MD_Numeric_Lumped_Mass_u_Prior_Interface: The value of beta_u must be specificed" << std::endl);
      }
      Set_beta_u(u_hyperparam_interface_->Get_beta_u());

      if (!u_hyperparam_interface_->Is_Transient())
      {
        if (u_hyperparam_interface_->Get_alpha_u() == 0.0)
        {
          determine_u_hyperparams_->Determine_alpha_u(this);
        }
        this->Set_alpha_u(u_hyperparam_interface_->Get_alpha_u());
      }
    }

    virtual ~MD_Lumped_Mass_u_Prior_Interface()
    {
    }

    void Set_beta_u(RealT beta_u_new)
    {
      beta_u_ = beta_u_new;
      E_u_->Set(*M_);
      E_u_->Scaled_Plus(beta_u_new, *S_);

      // E_u_solver_ = HDSA::makePtr<HDSA::Sparse_Matrix_Solver<RealT>>(E_u_);
      E_u_->Set_Symmetric();
      HDSA::Ptr<HDSA::Incomplete_Chol_Factor<RealT>> L = HDSA::makePtr<HDSA::Incomplete_Chol_Factor<RealT>>(E_u_);
      E_u_solver_ = HDSA::makePtr<HDSA::Sparse_Matrix_Solver<RealT>>(E_u_, false, false);
      E_u_solver_->Set_Incomplete_Factor(L);
      Assemble_W_u_Acute();
    }

    void Apply_E_u_Inverse(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &u_in) const
    {
      E_u_solver_->Apply_A_Inverse(u_out, u_in);
    }

    void Assemble_W_u_Acute(void)
    {
      int max_nonzeros_per_row = E_u_->Get_Max_Nonzeros_Per_Row();
      W_u_acute_ = M_->Clone(3 * max_nonzeros_per_row); // This may fail in spatial dimensions > 2. Need to test this and possibly modify.
      int n = W_u_acute_->Get_Number_of_Rows();

      W_u_acute_->Begin_Fill();

      std::vector<std::vector<int>> column_indices;
      column_indices.resize(n);
      std::vector<std::vector<RealT>> vals;
      vals.resize(n);
      for (int i = 0; i < n; i++)
      {
        E_u_->Get_Global_Row(i, column_indices[i], vals[i]);
      }

      for (int i = 0; i < n; i++)
      {
        std::vector<int> column_indices_i = column_indices[i];
        std::vector<RealT> vals_i = vals[i];
        int i_owner_rank = Get_Row_Owner_Rank(i);

        int col_indices_i_dim = column_indices_i.size();
        char *i_dim_buff = (char *)(&col_indices_i_dim);
        comm_->broadcast(i_owner_rank, 4, i_dim_buff);

        std::vector<int> col_indices_i;
        col_indices_i.resize(col_indices_i_dim);
        if (comm_->getRank() == i_owner_rank)
        {
          for (int ii = 0; ii < col_indices_i_dim; ii++)
          {
            col_indices_i[ii] = column_indices_i[ii];
          }
        }
        char *i_vals_buff = (char *)(&col_indices_i[0]);
        comm_->broadcast(i_owner_rank, 4 * col_indices_i_dim, i_vals_buff);

        std::vector<RealT> weights = std::vector<RealT>(col_indices_i_dim);
        for (int ii = 0; ii < col_indices_i_dim; ii++)
        {
          RealT weight_i = M_lumped_->Get_Entry(col_indices_i[ii]);
          char *send_buff = (char *)(&weight_i);
          std::vector<RealT> recv_vec = std::vector<RealT>(comm_->getSize(), 0.0);
          char *recv_buff = (char *)(&recv_vec[0]);
          comm_->gatherAll(8, send_buff, comm_->getSize() * 8, recv_buff);
          weights[ii] = std::accumulate(recv_vec.begin(), recv_vec.end(), 0.0);
        }

        for (int j = 0; j < n; j++)
        {
          std::vector<int> col_indices_j = column_indices[j];
          std::vector<RealT> vals_j = vals[j];
          int j_owner_rank = Get_Row_Owner_Rank(j);

          if (comm_->getRank() == j_owner_rank)
          {
            int dim = vals_j.size();
            char *dim_send_buff = (char *)(&dim);
            comm_->send(4, dim_send_buff, i_owner_rank);
            char *indices_send_buff = (char *)(&col_indices_j[0]);
            comm_->send(4 * dim, indices_send_buff, i_owner_rank);
            char *vals_send_buff = (char *)(&vals_j[0]);
            comm_->send(8 * dim, vals_send_buff, i_owner_rank);
          }

          if (comm_->getRank() == i_owner_rank)
          {
            int dim = 0;
            char *dim_recv_buff = (char *)(&dim);
            comm_->receive(j_owner_rank, 4, dim_recv_buff);
            col_indices_j.resize(dim);
            vals_j.resize(dim);
            char *indices_recv_buff = (char *)(&col_indices_j[0]);
            comm_->receive(j_owner_rank, 4 * dim, indices_recv_buff);
            char *val_recv_buff = (char *)(&vals_j[0]);
            comm_->receive(j_owner_rank, 8 * dim, val_recv_buff);
          }

          RealT val = 0.0;

          if (comm_->getRank() == i_owner_rank)
          {
            for (int ii = 0; ii < col_indices_i.size(); ii++)
            {
              for (int jj = 0; jj < col_indices_j.size(); jj++)
              {
                if (col_indices_i[ii] == col_indices_j[jj])
                {
                  val += vals_i[ii] * vals_j[jj] / weights[ii];
                  break;
                }
              }
            }
          }

          if (val != 0.0)
          {
            W_u_acute_->Set_Entry(i, j, val);
          }
        }
      }
      W_u_acute_->End_Fill();
    }

    int Get_Row_Owner_Rank(int i)
    {
      bool is_i_owned = E_u_->Is_Row_Owned(i);
      int i_owner_rank = 0;
      int send_int = 0;
      if (is_i_owned)
      {
        i_owner_rank = comm_->getRank();
        send_int = i_owner_rank;
      }
      char *send_buff = (char *)(&send_int);
      std::vector<int> recv_vec = std::vector<int>(comm_->getSize(), 0.0);
      char *recv_buff = (char *)(&recv_vec[0]);
      comm_->gatherAll(4, send_buff, comm_->getSize() * 4, recv_buff);
      if (!is_i_owned)
      {
        for (int j = 0; j < comm_->getSize(); j++)
        {
          i_owner_rank += recv_vec[j];
        }
      }
      return i_owner_rank;
    }
  };
}

#include "HDSA_MD_Determine_u_Hyperparameters_Def.hpp"

#endif
