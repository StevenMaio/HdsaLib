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

  public:
    void Apply_M_u(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &u_in) const
    {
      M_->Apply(u_out, u_in);
    }

    void Apply_W_u_Acute_Inverse(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &u_in) const
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

    void Apply_W_u_Acute_Plus_scalar_M_u_Inverse(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &u_in, const RealT &scalar) const
    {
      HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> A = W_u_acute_->Clone();
      A->Set(*W_u_acute_);
      A->Scaled_Plus(scalar, *M_);
      HDSA::Ptr<HDSA::Sparse_Matrix_Solver<RealT>> A_solver = HDSA::makePtr<HDSA::Sparse_Matrix_Solver<RealT>>(A);
      A_solver->Apply_A_Inverse(u_out, u_in);
    }

    void Sample_with_Covariance_W_u_Acute_Inverse(HDSA::MultiVector<RealT> &samples) const
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

    void Sample_with_Covariance_W_u_Acute_Plus_scalar_M_u_Inverse(HDSA::MultiVector<RealT> &samples, const RealT &scalar) const
    {
      HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> A = W_u_acute_->Clone();
      A->Set(*W_u_acute_);
      A->Scaled_Plus(scalar, *M_);
      HDSA::Ptr<HDSA::Sparse_Matrix_Sqrt<RealT>> A_inv_sqrt = HDSA::makePtr<HDSA::Sparse_Matrix_Sqrt<RealT>>(A, true);
      int M = samples.Number_of_Vectors();
      HDSA::Ptr<HDSA::Vector<RealT>> vec = samples[0]->Clone();
      for (int s = 0; s < M; s++)
      {
        vec->Randomize_Standard_Normal();
        A_inv_sqrt->Apply_Sqrt(*samples[s], *vec);
      }
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

    void Auxillary_Constructor()
    {
      HDSA::Ptr<HDSA::Vector<RealT>> u_tmp = data_interface_->Get_u_opt()->Clone();
      u_tmp->Set_Scalar(1.0);
      M_lumped_ = u_tmp->Clone();
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
      E_u_solver_ = HDSA::makePtr<HDSA::Sparse_Matrix_Solver<RealT>>(E_u_);
      Assemble_W_u_Acute();
    }

    void Apply_E_u_Inverse(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &u_in) const
    {
      E_u_solver_->Apply_A_Inverse(u_out, u_in);
    }

    void Assemble_W_u_Acute(void)
    {
      int max_nonzeros_per_row = E_u_->Get_Max_Nonzeros_Per_Row();
      W_u_acute_ = M_->Clone(2 * max_nonzeros_per_row); // This may fail in spatial dimensions > 2. Need to test this and possibly modify.
      int n = W_u_acute_->Get_Number_of_Rows();

      W_u_acute_->Begin_Fill();

      for (int i = 0; i < n; i++)
      {
        std::vector<int> column_indices_i;
        std::vector<RealT> vals_i;
        E_u_->Get_Global_Row(i, column_indices_i, vals_i);
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
          std::vector<int> col_indices_j;
          std::vector<RealT> vals_j;
          E_u_->Get_Global_Row(j, col_indices_j, vals_j);
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
