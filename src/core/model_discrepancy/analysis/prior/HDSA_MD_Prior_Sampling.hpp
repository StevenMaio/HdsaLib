#ifndef HDSA_MD_PRIOR_SAMPLING_HPP
#define HDSA_MD_PRIOR_SAMPLING_HPP

#include "HDSA_MD_Data_Interface.hpp"
#include "HDSA_MD_u_Prior_Interface.hpp"
#include "HDSA_MD_z_Prior_Interface.hpp"
#include "HDSA_MD_Elliptic_z_Prior_Interface.hpp"
#include "HDSA_MD_u_Hyperparameter_Interface.hpp"
#include "HDSA_MD_z_Hyperparameter_Interface.hpp"
#include "HDSA_MD_Multi_State_u_Hyperparameter_Interface.hpp"
#include "HDSA_Ensemble_Vector.hpp"
#include "HDSA_Linear_Algebra.hpp"

namespace HDSA
{

  template <class RealT>
  class MD_Prior_Sampling
  {

  private:
    const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> data_interface_;
    const HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT>> u_prior_interface_;
    const HDSA::Ptr<HDSA::MD_z_Prior_Interface<RealT>> z_prior_interface_;
    HDSA::Ptr<HDSA::MultiVector<RealT>> prior_delta_z_opt_;
    std::vector<std::vector<std::vector<RealT>>> prior_delta_z_opt_time_evol_; // [num_samps][num_time_steps][num_states]
    std::vector<std::vector<RealT>> prior_discrep_data_time_evol_;             // [num_time_steps][num_states]
    std::vector<HDSA::Ptr<HDSA::Vector<RealT>>> prior_z_pert_;
    std::vector<RealT> prior_z_pert_evals_;
    std::vector<HDSA::Ptr<HDSA::MultiVector<RealT>>> prior_delta_z_pert_;

  public:
    MD_Prior_Sampling(const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> &data_interface, const HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT>> &u_prior_interface, const HDSA::Ptr<HDSA::MD_z_Prior_Interface<RealT>> &z_prior_interface) : data_interface_(data_interface), u_prior_interface_(u_prior_interface), z_prior_interface_(z_prior_interface)
    {
    }

    virtual ~MD_Prior_Sampling()
    {
    }

    void Generate_Prior_Discrepancy_Sample_Data(int &num_samps, const HDSA::Ptr<HDSA::MD_u_Hyperparameter_Interface<RealT>> &u_hyperparam_interface, const HDSA::Ptr<HDSA::MD_z_Hyperparameter_Interface<RealT>> &z_hyperparam_interface, const HDSA::Ptr<HDSA::MultiVector<RealT>> &spatial_coords)
    {
      Generate_Prior_Discrepancy_z_opt_Sample_Data(num_samps, u_hyperparam_interface);
      Generate_Prior_Discrepancy_z_pert_Sample_Data(num_samps, z_hyperparam_interface, spatial_coords);
    }

    HDSA::Ptr<HDSA::MultiVector<RealT>> Get_prior_delta_z_opt(void) const
    {
      return prior_delta_z_opt_;
    }

    std::vector<std::vector<std::vector<RealT>>> Get_prior_delta_z_opt_time_evol(void) const
    {
      return prior_delta_z_opt_time_evol_;
    }

    std::vector<std::vector<RealT>> Get_prior_discrep_data_time_evol(void) const
    {
      return prior_discrep_data_time_evol_;
    }

    std::vector<HDSA::Ptr<HDSA::Vector<RealT>>> Get_prior_z_pert(void) const
    {
      return prior_z_pert_;
    }

    std::vector<HDSA::Ptr<HDSA::MultiVector<RealT>>> Get_prior_delta_z_pert(void) const
    {
      return prior_delta_z_pert_;
    }

    void Generate_Prior_Discrepancy_z_opt_Sample_Data(int &num_samps, const HDSA::Ptr<HDSA::MD_u_Hyperparameter_Interface<RealT>> &u_hyperparam_interface)
    {
      prior_delta_z_opt_ = Prior_Discrepancy_Samples_at_z_opt(num_samps);
      bool is_multi_state = u_hyperparam_interface->Is_Multi_State_Interface();
      int num_states = 1;
      if (is_multi_state)
      {
        HDSA::MD_Multi_State_u_Hyperparameter_Interface<RealT> *multi_state_u_hyperparam_interface = dynamic_cast<HDSA::MD_Multi_State_u_Hyperparameter_Interface<RealT> *>(&(*u_hyperparam_interface));
        num_states = multi_state_u_hyperparam_interface->Get_Number_of_States();
      }
      if (u_hyperparam_interface->Is_Transient())
      {
        prior_delta_z_opt_time_evol_.resize(num_samps);
        for (int k = 0; k < num_samps; k++)
        {
          HDSA::Ptr<HDSA::Vector<RealT>> dk = (*prior_delta_z_opt_)[k];

          HDSA::Ptr<HDSA::Vector<RealT>> dk_tmp;
          if (HDSA::Ensemble_Vector<RealT> *dk_ens = dynamic_cast<HDSA::Ensemble_Vector<RealT> *>(&(*dk)))
          {
            dk_tmp = dk_ens->Ensemble_Average();
          }
          else
          {
            dk_tmp = dk;
          }
          HDSA::Transient_Vector<RealT> *dk_trans = dynamic_cast<HDSA::Transient_Vector<RealT> *>(&(*dk_tmp));
          int n_t = dk_trans->Get_n_t();
          prior_delta_z_opt_time_evol_[k].resize(n_t);
          for (int i = 0; i < n_t; i++)
          {
            prior_delta_z_opt_time_evol_[k][i].resize(num_states);
            HDSA::Ptr<HDSA::Vector<RealT>> tmp = (*dk_trans)[i]->clone();
            u_prior_interface_->Apply_M_u(*tmp, *(*dk_trans)[i]);
            for (int j = 0; j < num_states; j++)
            {
              HDSA::Ptr<const HDSA::Vector<RealT>> dk_transi_j = data_interface_->Extract_State_Component(*(*dk_trans)[i], j);
              HDSA::Ptr<const HDSA::Vector<RealT>> tmp_j = data_interface_->Extract_State_Component(*tmp, j);
              prior_delta_z_opt_time_evol_[k][i][j] = std::sqrt(dk_transi_j->dot(*tmp_j));
            }
          }
        }
        int n_t = prior_delta_z_opt_time_evol_[0].size();
        prior_discrep_data_time_evol_.resize(n_t);

        HDSA::Ptr<HDSA::Vector<RealT>> D0 = (*data_interface_->get_D())[0];
        HDSA::Ptr<HDSA::Vector<RealT>> D0_tmp;
        if (HDSA::Ensemble_Vector<RealT> *D0_ens = dynamic_cast<HDSA::Ensemble_Vector<RealT> *>(&(*D0)))
        {
          D0_tmp = D0_ens->Ensemble_Average();
        }
        else
        {
          D0_tmp = D0;
        }

        HDSA::Transient_Vector<RealT> *d_trans = dynamic_cast<HDSA::Transient_Vector<RealT> *>(&(*D0_tmp));
        for (int i = 0; i < n_t; i++)
        {
          HDSA::Ptr<HDSA::Vector<RealT>> tmp = (*d_trans)[i]->clone();
          u_prior_interface_->Apply_M_u(*tmp, *(*d_trans)[i]);
          prior_discrep_data_time_evol_[i].resize(num_states);
          for (int j = 0; j < num_states; j++)
          {
            HDSA::Ptr<const HDSA::Vector<RealT>> d_transi_j = data_interface_->Extract_State_Component(*(*d_trans)[i], j);
            HDSA::Ptr<const HDSA::Vector<RealT>> tmp_j = data_interface_->Extract_State_Component(*tmp, j);
            prior_discrep_data_time_evol_[i][j] = std::sqrt(d_transi_j->dot(*tmp_j));
          }
        }
      }
    }

    void Generate_Prior_Discrepancy_z_pert_Sample_Data(int &num_samps, const HDSA::Ptr<HDSA::MD_z_Hyperparameter_Interface<RealT>> &z_hyperparam_interface, const HDSA::Ptr<HDSA::MultiVector<RealT>> &spatial_coords)
    {
      if (z_hyperparam_interface->Get_z_type() == "spatial field")
      {
        prior_z_pert_.resize(2);
        prior_z_pert_evals_.resize(2);
        prior_delta_z_pert_.resize(2);

        HDSA::Ptr<HDSA::Vector<RealT>> z_tmp = data_interface_->get_z_opt()->clone();
        z_prior_interface_->Apply_M_z(*z_tmp, *data_interface_->get_z_opt());
        RealT scaling = 0.3 * std::sqrt(data_interface_->get_z_opt()->dot(*z_tmp));

        prior_z_pert_[0] = data_interface_->get_z_opt()->clone();
        HDSA::Ptr<HDSA::Vector<RealT>> z1 = prior_z_pert_[0];
        z1->setScalar(1.0);
        z_tmp->zeros();
        z_prior_interface_->Apply_M_z(*z_tmp, *z1);
        RealT tmp1 = std::sqrt(z1->dot(*z_tmp));
        z1->scale(scaling / tmp1);
        prior_z_pert_evals_[0] = 1.0;

        prior_z_pert_[1] = data_interface_->get_z_opt()->clone();
        HDSA::Ptr<HDSA::Vector<RealT>> z2 = prior_z_pert_[1];
        int spatial_dim = spatial_coords->Number_of_Vectors();
        RealT omega = 0.0;
        std::vector<std::vector<RealT>> coord_bounds = z_hyperparam_interface->Spatial_Domain_Bounds();
        for (int j = 0; j < spatial_dim; j++)
        {
          omega += std::pow(coord_bounds[j][1] - coord_bounds[j][0], 2.0);
        }
        omega *= z_hyperparam_interface->Get_beta_z();
        omega = 1.5 * std::sqrt(1.0 / omega) / M_PI;
        omega = std::round(omega);
        int num_spatial_nodes = (*spatial_coords)[0]->dimension();
        for (int i = 0; i < num_spatial_nodes; i++)
        {
          RealT val = 1.0;
          for (int j = 0; j < spatial_dim; j++)
          {
            RealT x = ((*spatial_coords)[j]->Get_Entry(i) - coord_bounds[j][0]) / (coord_bounds[j][1] - coord_bounds[j][0]);
            val *= std::cos(2.0 * M_PI * omega * x);
          }
          z2->Set_Entry(i, val);
        }
        z_tmp->zeros();
        z_prior_interface_->Apply_M_z(*z_tmp, *z2);
        RealT tmp2 = std::sqrt(z2->dot(*z_tmp));
        z2->scale(1.0 / tmp2);

        MD_Elliptic_z_Prior_Interface<RealT> *elliptic_z_prior_interface = dynamic_cast<MD_Elliptic_z_Prior_Interface<RealT> *>(&(*z_prior_interface_));
        elliptic_z_prior_interface->Apply_E_z(*z_tmp, *z2);
        prior_z_pert_evals_[1] = 1.0 / (z2->dot(*z_tmp));
        z2->scale(scaling);

        for (int k = 0; k < 2; k++)
        {
          prior_delta_z_pert_[k] = HDSA::makePtr<HDSA::MultiVector<RealT>>(num_samps, *data_interface_->get_u_opt());
          u_prior_interface_->Sample_with_Covariance_W_u_Inverse(*prior_delta_z_pert_[k]);
          prior_delta_z_pert_[k]->scale(scaling * std::sqrt(z_hyperparam_interface->Get_alpha_z()) * prior_z_pert_evals_[k]);
        }
      }
      else
      {
        prior_z_pert_.resize(1);
        prior_delta_z_pert_.resize(1);

        HDSA::Ptr<HDSA::Vector<RealT>> z_tmp = data_interface_->get_z_opt()->clone();
        z_prior_interface_->Apply_M_z(*z_tmp, *data_interface_->get_z_opt());
        RealT scaling = 0.3 * std::sqrt(data_interface_->get_z_opt()->dot(*z_tmp));

        prior_z_pert_[0] = data_interface_->get_z_opt()->clone();
        HDSA::Ptr<HDSA::Vector<RealT>> z1 = prior_z_pert_[0];
        z1->setScalar(1.0);
        z_tmp->zeros();
        z_prior_interface_->Apply_M_z(*z_tmp, *z1);
        RealT tmp1 = std::sqrt(z1->dot(*z_tmp));
        z1->scale(scaling / tmp1);

        prior_delta_z_pert_[0] = HDSA::makePtr<HDSA::MultiVector<RealT>>(num_samps, *data_interface_->get_u_opt());
        u_prior_interface_->Sample_with_Covariance_W_u_Inverse(*prior_delta_z_pert_[0]);
        prior_delta_z_pert_[0]->scale(scaling * std::sqrt(z_hyperparam_interface->Get_alpha_z()));
      }
    }

    HDSA::Ptr<HDSA::MultiVector<RealT>> Prior_Discrepancy_Samples_at_z_opt(const int &num_samples)
    {
      HDSA::Ptr<HDSA::MultiVector<RealT>> delta_samples = HDSA::makePtr<HDSA::MultiVector<RealT>>(num_samples, *data_interface_->get_u_opt());
      u_prior_interface_->Sample_with_Covariance_W_u_Inverse(*delta_samples);
      for (int k = 0; k < delta_samples->Number_of_Vectors(); k++)
      {
        (*delta_samples)[k]->plus(*data_interface_->get_data_shift());
      }
      return delta_samples;
    }

    std::vector<HDSA::Ptr<HDSA::MultiVector<RealT>>> Prior_Discrepancy_Samples(const HDSA::MultiVector<RealT> &z, const int &num_samples)
    {
      int N = z.Number_of_Vectors();

      std::vector<HDSA::Ptr<HDSA::MultiVector<RealT>>> delta_samples;
      delta_samples.resize(num_samples);

      // Compute Sigma = Z^T*M_z*W_z_inv*M_z*Z
      HDSA::Ptr<HDSA::MultiVector<RealT>> Z = HDSA::makePtr<HDSA::MultiVector<RealT>>(N, *data_interface_->get_z_opt());
      HDSA::Ptr<HDSA::MultiVector<RealT>> M_z_Z = HDSA::makePtr<HDSA::MultiVector<RealT>>(N, *data_interface_->get_z_opt());
      HDSA::Ptr<HDSA::MultiVector<RealT>> W_z_inv_M_z_Z = HDSA::makePtr<HDSA::MultiVector<RealT>>(N, *data_interface_->get_z_opt());
      for (int k = 0; k < N; k++)
      {
        HDSA::Ptr<HDSA::Vector<RealT>> z_tmp1 = (*Z)[k];
        z_tmp1->set(*z[k]);
        z_tmp1->axpy(-1.0, *data_interface_->get_z_opt());

        HDSA::Ptr<HDSA::Vector<RealT>> z_tmp2 = (*M_z_Z)[k];
        z_prior_interface_->Apply_M_z(*z_tmp2, *z_tmp1);

        HDSA::Ptr<HDSA::Vector<RealT>> z_tmp3 = (*W_z_inv_M_z_Z)[k];
        z_prior_interface_->Apply_W_z_Inverse(*z_tmp3, *z_tmp2);
      }
      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> Sigma = M_z_Z->MatMat(*W_z_inv_M_z_Z);

      // Factorize Sigma = R*R^T
      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> R = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(N, N);
      HDSA::Linear_Algebra::Cholesky_Factorization(*Sigma, *R);

      // Loop over samples
      for (int i = 0; i < num_samples; i++)
      {
        delta_samples[i] = HDSA::makePtr<HDSA::MultiVector<RealT>>(N, *data_interface_->get_u_opt());

        // Generate random sample i
        HDSA::Ptr<HDSA::MultiVector<RealT>> u_rand = HDSA::makePtr<HDSA::MultiVector<RealT>>(N + 1, *data_interface_->get_u_opt());
        u_prior_interface_->Sample_with_Covariance_W_u_Inverse(*u_rand);

        for (int k = 0; k < N; k++)
        {
          HDSA::Ptr<HDSA::Vector<RealT>> u_out = (*delta_samples[i])[k];
          u_out->set(*(*u_rand)[N]);
          for (int j = 0; j < k + 1; j++)
          {
            u_out->axpy((*R)(j, k), *(*u_rand)[j]);
          }
          u_out->plus(*data_interface_->get_data_shift());
        }
      }

      return delta_samples;
    }
  };

}

#endif
