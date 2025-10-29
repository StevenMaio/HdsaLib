#ifndef HDSA_MD_DETERMINE_U_HYPERPARAMETERS_DEF_HPP
#define HDSA_MD_DETERMINE_U_HYPERPARAMETERS_DEF_HPP

#include "HDSA_MD_Determine_u_Hyperparameters_Decl.hpp"
#include "HDSA_MD_Data_Interface.hpp"
#include "HDSA_MD_Numeric_Laplacian_u_Prior_Interface.hpp"
#include "HDSA_MD_Lumped_Mass_u_Prior_Interface.hpp"
#include "HDSA_MD_Transient_Elliptic_u_Prior_Interface.hpp"
#include "HDSA_MD_Laplacian_Like_Operator_Properties.hpp"

namespace HDSA
{

  template <class RealT>
  MD_Determine_u_Hyperparameters<RealT>::MD_Determine_u_Hyperparameters(const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> &data_interface, const HDSA::Ptr<HDSA::MD_u_Hyperparameter_Interface<RealT>> &u_hyperparam_interface) : data_interface_(data_interface), u_hyperparam_interface_(u_hyperparam_interface)
  {
    data_interface_->Load_Data();
    if (u_hyperparam_interface_->Center_Data())
    {
      data_interface_->Center_Data();
    }
    component_id_ = u_hyperparam_interface_->Get_Component_ID();
    is_transient_ = u_hyperparam_interface_->Is_Transient();
  }

  template <class RealT>
  MD_Determine_u_Hyperparameters<RealT>::~MD_Determine_u_Hyperparameters()
  {
  }

  template <class RealT>
  void MD_Determine_u_Hyperparameters<RealT>::Determine_alpha_u(HDSA::MD_u_Prior_Interface<RealT> *u_prior_interface) const
  {
    HDSA::Ptr<const HDSA::MultiVector<RealT>> D = data_interface_->Get_D();
    HDSA::Ptr<HDSA::Vector<RealT>> delta = (*D)[0];
    HDSA::Ptr<const HDSA::Vector<RealT>> delta_k = data_interface_->Extract_State_Component(*delta, component_id_, true);
    HDSA::Ptr<HDSA::Vector<RealT>> tmp1 = delta_k->Clone();
    u_prior_interface->Apply_M_u(*tmp1, *delta_k);
    RealT d1_Norm_sq = delta_k->Dot(*tmp1);

    RealT u_op_trace = 0.0;
    if (is_transient_)
    {

      MD_Transient_Elliptic_u_Prior_Interface<RealT> u_transient = dynamic_cast<MD_Transient_Elliptic_u_Prior_Interface<RealT> &>(*u_prior_interface);
      HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT>> u_spatial = u_transient.Get_Spatial_Cov();

      if (HDSA::Ptr<HDSA::MD_Numeric_Laplacian_u_Prior_Interface<RealT>> u_elliptic = HDSA::dynamicPtrCast<HDSA::MD_Numeric_Laplacian_u_Prior_Interface<RealT>>(u_spatial))
      {
        HDSA::Ptr<HDSA::Dense_Matrix<RealT>> sing_vals = u_elliptic->Get_Sing_Vals();
        int m = sing_vals->Number_of_Rows();
        for (int k = 0; k < m; k++)
        {
          u_op_trace += std::pow((*sing_vals)(k, 0), 2.0);
        }
      }
      else if (HDSA::MD_Lumped_Mass_u_Prior_Interface<RealT> *u_lumped = dynamic_cast<HDSA::MD_Lumped_Mass_u_Prior_Interface<RealT> *>(&(*u_prior_interface)))
      {
        HDSA::Ptr<HDSA::MD_Laplacian_Like_Operator_Properties<RealT>> laplacian_op_prop = HDSA::makePtr<HDSA::MD_Laplacian_Like_Operator_Properties<RealT>>();
        int trace_estimator_sample_size = u_hyperparam_interface_->Get_trace_estimator_sample_size();
        if (trace_estimator_sample_size > 0)
        {
          HDSA::Ptr<HDSA::Vector<RealT>> u_vec = data_interface_->Get_u_opt()->Clone();
          u_op_trace = laplacian_op_prop->Randomized_Inv_Operator_Trace_Estimation(*u_lumped, u_vec, trace_estimator_sample_size);
        }
        else
        {
          int n_u = data_interface_->Extract_State_Component(*data_interface_->Get_u_opt(), component_id_)->Dimension();
          RealT beta_u = u_hyperparam_interface_->Get_beta_u();
          std::vector<std::vector<RealT>> spatial_domain_bounds = u_hyperparam_interface_->Spatial_Domain_Bounds();
          u_op_trace = laplacian_op_prop->Get_Rectangular_Domain_Squared_Inv_Operator_Trace(beta_u, spatial_domain_bounds, n_u);
        }
      }

      HDSA::Ptr<HDSA::MD_Transient_Prior_Covariance<RealT>> u_time = u_transient.Get_Time_Cov();
      RealT tmp = 0.0;
      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> evals = u_time->Get_Evals();
      int n = evals->Number_of_Rows();
      for (int k = 0; k < n; k++)
      {
        tmp += (*evals)(k, 0);
      }
      u_op_trace *= tmp;
    }
    else
    {

      if (HDSA::MD_Numeric_Laplacian_u_Prior_Interface<RealT> *u_elliptic = dynamic_cast<HDSA::MD_Numeric_Laplacian_u_Prior_Interface<RealT> *>(&(*u_prior_interface)))
      {
        HDSA::Ptr<HDSA::Dense_Matrix<RealT>> sing_vals = u_elliptic->Get_Sing_Vals();
        int m = sing_vals->Number_of_Rows();
        for (int k = 0; k < m; k++)
        {
          u_op_trace += std::pow((*sing_vals)(k, 0), 2.0);
        }
      }
      else if (HDSA::MD_Lumped_Mass_u_Prior_Interface<RealT> *u_lumped = dynamic_cast<HDSA::MD_Lumped_Mass_u_Prior_Interface<RealT> *>(&(*u_prior_interface)))
      {
        HDSA::Ptr<HDSA::MD_Laplacian_Like_Operator_Properties<RealT>> laplacian_op_prop = HDSA::makePtr<HDSA::MD_Laplacian_Like_Operator_Properties<RealT>>();
        int trace_estimator_sample_size = u_hyperparam_interface_->Get_trace_estimator_sample_size();
        if (trace_estimator_sample_size > 0)
        {
          HDSA::Ptr<HDSA::Vector<RealT>> u_vec = data_interface_->Get_u_opt()->Clone();
          u_op_trace = laplacian_op_prop->Randomized_Inv_Operator_Trace_Estimation(*u_lumped, u_vec, trace_estimator_sample_size);
        }
        else
        {
          int n_u = data_interface_->Extract_State_Component(*data_interface_->Get_u_opt(), component_id_)->Dimension();
          RealT beta_u = u_hyperparam_interface_->Get_beta_u();
          std::vector<std::vector<RealT>> spatial_domain_bounds = u_hyperparam_interface_->Spatial_Domain_Bounds();
          u_op_trace = laplacian_op_prop->Get_Rectangular_Domain_Squared_Inv_Operator_Trace(beta_u, spatial_domain_bounds, n_u);
        }
      }
    }

    RealT alpha_u_new = d1_Norm_sq / u_op_trace;
    u_hyperparam_interface_->Set_alpha_u(alpha_u_new);
  }

  template <class RealT>
  void MD_Determine_u_Hyperparameters<RealT>::Determine_alpha_t(HDSA::MD_u_Prior_Interface<RealT> *u_prior_interface) const
  {
    HDSA::Ptr<const HDSA::MultiVector<RealT>> D = data_interface_->Get_D();
    HDSA::Transient_Vector<RealT> tmp = dynamic_cast<HDSA::Transient_Vector<RealT> &>(*(*D)[0]);
    int n_t = tmp.Get_n_t();
    int N = D->Number_of_Vectors();
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> weights = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(n_t, N);
    for (int j = 0; j < N; j++)
    {
      HDSA::Transient_Vector<RealT> d_trans = dynamic_cast<HDSA::Transient_Vector<RealT> &>(*(*D)[j]);
      RealT max_j = 0.0;
      for (int k = 0; k < n_t; k++)
      {
        HDSA::Ptr<const HDSA::Vector<RealT>> d = data_interface_->Extract_State_Component(*d_trans[k], component_id_);
        HDSA::Ptr<HDSA::Vector<RealT>> vec = d->Clone();
        u_prior_interface->Apply_M_u(*vec, *d);
        RealT val = vec->Dot(*d);
        max_j = std::max(max_j, val);
        weights->Set_Entry(k, j, val);
      }
      for (int k = 0; k < n_t; k++)
      {
        RealT val = (*weights)(k, j) / max_j;
        val = (val + u_hyperparam_interface_->Get_Time_Variance_Inflation()) / (1.0 + u_hyperparam_interface_->Get_Time_Variance_Inflation());
        weights->Set_Entry(k, j, val);
      }
    }
    std::vector<RealT> alpha_t_new = std::vector<RealT>(n_t, 0.0);
    for (int k = 0; k < n_t; k++)
    {
      RealT val = 0.0;
      for (int j = 0; j < N; j++)
      {
        val += (*weights)(k, j);
      }
      alpha_t_new[k] = val / static_cast<RealT>(N);
    }
    u_hyperparam_interface_->Set_alpha_t(alpha_t_new);
  }

  template <class RealT>
  void MD_Determine_u_Hyperparameters<RealT>::Determine_alpha_d(HDSA::MD_u_Prior_Interface<RealT> *u_prior_interface) const
  {
    HDSA::Ptr<const HDSA::MultiVector<RealT>> D = data_interface_->Get_D();
    int N = D->Number_of_Vectors();
    RealT mags = 0.0;
    for (int j = 0; j < N; j++)
    {
      HDSA::Ptr<HDSA::Vector<RealT>> delta_tmp = (*D)[j];
      HDSA::Ptr<HDSA::Vector<RealT>> delta;
      if (HDSA::Ptr<HDSA::Transient_Vector<RealT>> delta_tmp_trans = HDSA::dynamicPtrCast<HDSA::Transient_Vector<RealT>>(delta_tmp))
      {
        delta = (*delta_tmp_trans)[0]->Clone();
        for (int k = 0; k < delta_tmp_trans->Get_n_t(); k++)
        {
          delta->Plus(*(*delta_tmp_trans)[k]);
        }
        delta->Scale(1.0 / static_cast<RealT>(delta_tmp_trans->Get_n_t()));
      }
      else
      {
        delta = delta_tmp;
      }
      HDSA::Ptr<const HDSA::Vector<RealT>> delta_k = data_interface_->Extract_State_Component(*delta, component_id_, true);
      HDSA::Ptr<HDSA::Vector<RealT>> tmp1 = delta_k->Clone();
      u_prior_interface->Apply_M_u(*tmp1, *delta_k);
      mags += std::sqrt(tmp1->Dot(*delta_k));
    }
    mags = mags / static_cast<RealT>(N);
    RealT alpha_d_new = std::pow(mags * u_hyperparam_interface_->Get_Data_Noise_Percent(), 2.0);
    u_hyperparam_interface_->Set_alpha_d(alpha_d_new);
  }

  template <class RealT>
  void MD_Determine_u_Hyperparameters<RealT>::Determine_GSVD_Hyperparameters(void) const
  {
    std::vector<std::vector<RealT>> spatial_domain_bounds = u_hyperparam_interface_->Spatial_Domain_Bounds();
    int spatial_dim = spatial_domain_bounds.size();
    std::vector<int> modes_per_dim = std::vector<int>(spatial_dim);
    int num_sing_vals = 1;
    int n_u = data_interface_->Extract_State_Component(*data_interface_->Get_u_opt(), 0)->Dimension();
    for (int k = 0; k < spatial_dim; k++)
    {
      RealT L = spatial_domain_bounds[k][1] - spatial_domain_bounds[k][0];
      RealT tmp = std::pow(L / M_PI, 2.0) * (1.0 / u_hyperparam_interface_->Get_beta_u()) * (1.0 / u_hyperparam_interface_->Get_W_u_Inv_Spectal_Gap() - 1.0);
      modes_per_dim[k] = std::round(std::sqrt(tmp));
      num_sing_vals *= modes_per_dim[k];
    }
    num_sing_vals = std::min(num_sing_vals, n_u);
    int oversampling = std::min(10, n_u - num_sing_vals);

    u_hyperparam_interface_->Set_GSVD_Hyperparameters(num_sing_vals, oversampling, 1);
  }

}

#endif
