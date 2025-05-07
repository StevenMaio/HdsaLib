#ifndef HDSA_MD_DETERMINE_Z_HYPERPARAMETERS_DEF_HPP
#define HDSA_MD_DETERMINE_Z_HYPERPARAMETERS_DEF_HPP

namespace HDSA
{

  template <class RealT>
  MD_Determine_z_Hyperparameters<RealT>::MD_Determine_z_Hyperparameters(const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> &data_interface, const HDSA::Ptr<HDSA::MD_z_Hyperparameter_Interface<RealT>> &z_hyperparam_interface, const HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT>> &u_prior_interface) : data_interface_(data_interface), z_hyperparam_interface_(z_hyperparam_interface), u_prior_interface_(u_prior_interface)
  {
  }

  template <class RealT>
  MD_Determine_z_Hyperparameters<RealT>::~MD_Determine_z_Hyperparameters()
  {
  }

  template <class RealT>
  void MD_Determine_z_Hyperparameters<RealT>::Determine_alpha_z(HDSA::MD_z_Prior_Interface<RealT> *z_prior_interface) const
  {
    std::cout << "Need to test Determine_alpha_z" << std::endl;
    HDSA::Ptr<const HDSA::Vector<RealT>> z_opt = data_interface_->get_z_opt();
    HDSA::Ptr<HDSA::Vector<RealT>> tmp = z_opt->clone();
    z_prior_interface->Apply_M_z(*tmp, *z_opt);
    RealT z_opt_norm = std::sqrt(tmp->dot(*z_opt));

    std::vector<RealT> evals = Compute_Eigenvalues(z_prior_interface);
    std::vector<RealT> truncated_evals;
    RealT threshold = 1.e-3;
    std::copy_if(evals.begin(), evals.end(), std::back_inserter(truncated_evals), [threshold](RealT value)
                 { return value >= threshold; });

    int samples = 1000;
    int e_dim = truncated_evals.size();
    RealT MC_est = 0.0;
    for (int k = 0; k < samples; k++)
    {
      RealT num = 0.0;
      RealT denom = 0.0;
      for (int i = 0; i < e_dim; i++)
      {
        RealT nu = std::pow(z_hyperparam_interface_->Get_Random_Number_Generator()->Generate_Standard_Normal_Sample(), 2.0);
        num += nu * std::pow(truncated_evals[i], 4.0);
        denom += nu * std::pow(truncated_evals[i], 2.0);
      }
      MC_est += num / denom;
    }
    MC_est = MC_est / static_cast<RealT>(samples);

    if (z_hyperparam_interface_->Get_Discrepancy_Percent_z_Variation() == 1.0)
    {
      int num_state_solves = z_hyperparam_interface_->Get_Num_State_Solves();
      if (num_state_solves > 0)
      {
        HDSA::Ptr<const HDSA::Vector<RealT>> u_nom = data_interface_->get_u_opt();
        HDSA::Ptr<HDSA::MultiVector<RealT>> z_samples = HDSA::makePtr<HDSA::MultiVector<RealT>>(num_state_solves, *data_interface_->get_z_opt());
        MD_Numeric_Laplacian_z_Prior_Interface<RealT> z_elliptic = dynamic_cast<MD_Numeric_Laplacian_z_Prior_Interface<RealT> &>(*z_prior_interface);
        z_elliptic.Sample_with_Covariance_W_z_Acute_Inverse(*z_samples);
        HDSA::Ptr<HDSA::Vector<RealT>> z_tmp = data_interface_->get_z_opt()->clone();
        for (int k = 0; k < num_state_solves; k++)
        {
          z_elliptic.Apply_M_z(*z_tmp, *(*z_samples)[k]);
          RealT val = std::sqrt(z_tmp->dot(*(*z_samples)[k]));
          (*z_samples)[k]->scale(z_opt_norm / val);
          (*z_samples)[k]->plus(*data_interface_->get_z_opt());
        }
        HDSA::Ptr<HDSA::MultiVector<RealT>> u_samples = HDSA::makePtr<HDSA::MultiVector<RealT>>(num_state_solves, *data_interface_->get_u_opt());
        for (int k = 0; k < num_state_solves; k++)
        {
          z_hyperparam_interface_->State_Solve(*(*u_samples)[k], *(*z_samples)[k]);
          (*u_samples)[k]->axpy(-1.0, *data_interface_->get_u_opt());
        }
        std::vector<RealT> e_norm_sq = std::vector<RealT>(num_state_solves, 0.0);
        RealT val = 0.0;
        HDSA::Ptr<HDSA::Vector<RealT>> u_tmp = data_interface_->get_u_opt()->clone();
        HDSA::Ptr<const HDSA::Vector<RealT>> d1 = (*data_interface_->get_D())[0];
        u_prior_interface_->Apply_M_u(*u_tmp, *d1);
        RealT d1_norm_sq = u_tmp->dot(*d1);
        for (int k = 0; k < num_state_solves; k++)
        {
          u_tmp->zeros();
          u_prior_interface_->Apply_M_u(*u_tmp, *(*u_samples)[k]);
          val += (u_tmp->dot(*(*u_samples)[k])) / d1_norm_sq;
        }
        val = val / static_cast<RealT>(num_state_solves);
        val = std::sqrt(val);
        z_hyperparam_interface_->Set_Discrepancy_Percent_z_Variation(val);
      }
      else
      {
        int N = data_interface_->get_D()->Number_of_Vectors();
        if (N > 1)
        {
          HDSA::Ptr<const HDSA::Vector<RealT>> z1 = (*data_interface_->get_Z())[0];
          HDSA::Ptr<HDSA::Vector<RealT>> z_tmp = z1->clone();
          z_prior_interface->Apply_M_z(*z_tmp, *z1);
          RealT z1_norm_sq = z_tmp->dot(*z1);
          std::vector<RealT> z_pert_norm_sq = std::vector<RealT>(N - 1, 0.0);
          for (int k = 0; k < N - 1; k++)
          {
            HDSA::Ptr<HDSA::Vector<RealT>> zk = z1->clone();
            zk->set(*(*data_interface_->get_Z())[k + 1]);
            zk->axpy(-1.0, *z1);
            z_prior_interface->Apply_M_z(*z_tmp, *zk);
            z_pert_norm_sq[k] = z_tmp->dot(*zk);
          }

          HDSA::Ptr<const HDSA::Vector<RealT>> d1 = (*data_interface_->get_D())[0];
          HDSA::Ptr<HDSA::Vector<RealT>> d_tmp = d1->clone();
          u_prior_interface_->Apply_M_u(*d_tmp, *d1);
          RealT d1_norm_sq = d_tmp->dot(*d1);
          std::vector<RealT> d_pert_norm_sq = std::vector<RealT>(N - 1, 0.0);
          for (int k = 0; k < N - 1; k++)
          {
            HDSA::Ptr<HDSA::Vector<RealT>> dk = d1->clone();
            dk->set(*(*data_interface_->get_D())[k + 1]);
            dk->axpy(-1.0, *d1);
            u_prior_interface_->Apply_M_u(*d_tmp, *dk);
            d_pert_norm_sq[k] = d_tmp->dot(*dk);
          }
          RealT val = 0.0;
          for (int k = 0; k < N - 1; k++)
          {
            RealT valk = (d_pert_norm_sq[k] / d1_norm_sq) / (z_pert_norm_sq[k] / z1_norm_sq);
            val += std::sqrt(valk);
          }
          val = val / static_cast<RealT>(N - 1);
          z_hyperparam_interface_->Set_Discrepancy_Percent_z_Variation(val);
        }
      }

      RealT alpha_z_new = std::pow(z_hyperparam_interface_->Get_Discrepancy_Percent_z_Variation(), 2.0) / (MC_est * std::pow(z_opt_norm, 2.0));
      z_hyperparam_interface_->Set_alpha_z(alpha_z_new);
    }
  }

  template <class RealT>
  std::vector<RealT> MD_Determine_z_Hyperparameters<RealT>::Compute_Eigenvalues(HDSA::MD_z_Prior_Interface<RealT> *z_prior_interface) const
  {
    std::vector<RealT> evals;
    if (z_hyperparam_interface_->Get_z_type() == "spatial field")
    {
      std::vector<std::vector<RealT>> spatial_bounds = z_hyperparam_interface_->Spatial_Domain_Bounds();
      int d = spatial_bounds.size();
      int n_z = data_interface_->get_z_opt()->dimension();
      if (d == 1)
      {
        RealT Lx = spatial_bounds[0][1] - spatial_bounds[0][0];
        int n = n_z - 1;
        evals.resize(n + 1);
        for (int i = 0; i < n + 1; i++)
        {
          evals[i] = 1.0 / (1.0 + z_hyperparam_interface_->Get_beta_z() * std::pow(M_PI / Lx, 2.0) * std::pow(static_cast<RealT>(i), 2.0));
        }
      }
      else if (d == 2)
      {
        RealT Lx = spatial_bounds[0][1] - spatial_bounds[0][0];
        RealT Ly = spatial_bounds[1][1] - spatial_bounds[1][0];
        int n = std::round(std::sqrt(n_z)) - 1;
        evals.resize((n + 1) * (n + 1));
        int count = 0;
        for (int i = 0; i < n + 1; i++)
        {
          for (int j = 0; j < n + 1; j++)
          {
            evals[count] = 1.0 / (1.0 + z_hyperparam_interface_->Get_beta_z() * (std::pow(M_PI / Lx, 2.0) * std::pow(static_cast<RealT>(i), 2.0) + std::pow(M_PI / Ly, 2.0) * std::pow(static_cast<RealT>(j), 2.0)));
            count += 1;
          }
        }
      }
      else if (d == 3)
      {
        RealT Lx = spatial_bounds[0][1] - spatial_bounds[0][0];
        RealT Ly = spatial_bounds[1][1] - spatial_bounds[1][0];
        RealT Lz = spatial_bounds[2][1] - spatial_bounds[2][0];
        int n = std::round(std::pow(n_z, 1.0 / 3.0)) - 1;
        evals.resize((n + 1) * (n + 1) * (n + 1));
        int count = 0;
        for (int i = 0; i < n + 1; i++)
        {
          for (int j = 0; j < n + 1; j++)
          {
            for (int k = 0; k < n + 1; k++)
            {
              evals[count] = 1.0 / (1.0 + z_hyperparam_interface_->Get_beta_z() * (std::pow(M_PI / Lx, 2.0) * std::pow(static_cast<RealT>(i), 2.0) + std::pow(M_PI / Ly, 2.0) * std::pow(static_cast<RealT>(j), 2.0)) + std::pow(M_PI / Lz, 2.0) * std::pow(static_cast<RealT>(k), 2.0));
              count += 1;
            }
          }
        }
      }
      else
      {
        std::cout << "Error: the code only supports spatial dimensions s=1,2,3" << std::endl;
      }
      std::sort(evals.begin(), evals.end(), std::greater<RealT>());
    }
    else
    {
      std::cout << "Error: the code only supports z_type=spatial field at this time." << std::endl;
    }
    return evals;
  }

  template <class RealT>
  void MD_Determine_z_Hyperparameters<RealT>::Determine_beta_z(void) const
  {
    RealT beta_z_new = 0.008;
    z_hyperparam_interface_->Set_beta_z(beta_z_new);
  }

  template <class RealT>
  void MD_Determine_z_Hyperparameters<RealT>::Determine_beta_t(void) const
  {
    RealT beta_t_new = 0.027;
    z_hyperparam_interface_->Set_beta_t(beta_t_new);
  }

}

#endif
