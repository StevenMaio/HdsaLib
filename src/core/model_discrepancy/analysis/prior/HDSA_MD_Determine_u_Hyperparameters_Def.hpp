#ifndef HDSA_MD_DETERMINE_U_HYPERPARAMETERS_DEF_HPP
#define HDSA_MD_DETERMINE_U_HYPERPARAMETERS_DEF_HPP

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
    HDSA::Ptr<const HDSA::MultiVector<RealT>> D = data_interface_->get_D();
    HDSA::Ptr<HDSA::Vector<RealT>> delta = (*D)[0];
    HDSA::Ptr<const HDSA::Vector<RealT>> delta_k = data_interface_->Extract_State_Component(*delta, component_id_, true);
    HDSA::Ptr<HDSA::Vector<RealT>> tmp1 = delta_k->clone();
    u_prior_interface->Apply_M_u(*tmp1, *delta_k);
    RealT d1_norm_sq = delta_k->dot(*tmp1);

    RealT u_op_trace = 0.0;
    if (is_transient_)
    {
      MD_Transient_Elliptic_u_Prior_Interface<RealT> u_transient = dynamic_cast<MD_Transient_Elliptic_u_Prior_Interface<RealT> &>(*u_prior_interface);
      HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT>> u_spatial = u_transient.Get_Spatial_Cov();
      MD_Numeric_Laplacian_u_Prior_Interface<RealT> u_elliptic = dynamic_cast<MD_Numeric_Laplacian_u_Prior_Interface<RealT> &>(*u_spatial);

      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> sing_vals = u_elliptic.Get_Sing_Vals();
      int m = sing_vals->numRows();
      for (int k = 0; k < m; k++)
      {
        u_op_trace += std::pow((*sing_vals)(k, 0), 2.0);
      }

      HDSA::Ptr<HDSA::MD_Transient_Prior_Covariance<RealT>> u_time = u_transient.Get_Time_Cov();
      RealT tmp = 0.0;
      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> evals = u_time->Get_Evals();
      int n = evals->numRows();
      for (int k = 0; k < n; k++)
      {
        tmp += (*evals)(k, 0);
      }
      u_op_trace *= tmp;
    }
    else
    {
      MD_Numeric_Laplacian_u_Prior_Interface<RealT> u_elliptic = dynamic_cast<MD_Numeric_Laplacian_u_Prior_Interface<RealT> &>(*u_prior_interface);
      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> sing_vals = u_elliptic.Get_Sing_Vals();
      int m = sing_vals->numRows();
      for (int k = 0; k < m; k++)
      {
        u_op_trace += std::pow((*sing_vals)(k, 0), 2.0);
      }
    }

    RealT alpha_u_new = d1_norm_sq / u_op_trace;
    u_hyperparam_interface_->Set_alpha_u(alpha_u_new);
  }

  template <class RealT>
  void MD_Determine_u_Hyperparameters<RealT>::Determine_alpha_t(HDSA::MD_u_Prior_Interface<RealT> *u_prior_interface) const
  {
    HDSA::Ptr<const HDSA::MultiVector<RealT>> D = data_interface_->get_D();
    Transient_Vector<RealT> tmp = dynamic_cast<Transient_Vector<RealT> &>(*(*D)[0]);
    int n_t = tmp.Get_n_t();
    int N = D->Number_of_Vectors();
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> weights = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(n_t, N);
    for (int j = 0; j < N; j++)
    {
      Transient_Vector<RealT> d_trans = dynamic_cast<Transient_Vector<RealT> &>(*(*D)[j]);
      RealT max_j = 0.0;
      for (int k = 0; k < n_t; k++)
      {
        HDSA::Ptr<const HDSA::Vector<RealT>> d = data_interface_->Extract_State_Component(*d_trans[k], component_id_);
        HDSA::Ptr<HDSA::Vector<RealT>> vec = d->clone();
        u_prior_interface->Apply_M_u(*vec, *d);
        RealT val = vec->dot(*d);
        max_j = std::max(max_j, val);
        weights->Replace_Element(k, j, val);
      }
      for (int k = 0; k < n_t; k++)
      {
        RealT val = (*weights)(k, j) / max_j;
        val = (val + u_hyperparam_interface_->Get_Time_Variance_Inflation()) / (1.0 + u_hyperparam_interface_->Get_Time_Variance_Inflation());
        weights->Replace_Element(k, j, val);
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
    HDSA::Ptr<const HDSA::MultiVector<RealT>> D = data_interface_->get_D();
    int N = D->Number_of_Vectors();
    RealT mags = 0.0;
    for (int j = 0; j < N; j++)
    {
      HDSA::Ptr<HDSA::Vector<RealT>> delta = (*D)[j];
      HDSA::Ptr<const HDSA::Vector<RealT>> delta_k = data_interface_->Extract_State_Component(*delta, component_id_, true);
      HDSA::Ptr<HDSA::Vector<RealT>> tmp1 = delta_k->clone();
      u_prior_interface->Apply_M_u(*tmp1, *delta_k);
      mags += std::sqrt(tmp1->dot(*delta_k));
    }
    mags = mags / static_cast<RealT>(N);
    RealT alpha_d_new = std::pow(mags * u_hyperparam_interface_->Get_Data_Noise_Percent(), 2.0);
    u_hyperparam_interface_->Set_alpha_d(alpha_d_new);
  }

  template <class RealT>
  void MD_Determine_u_Hyperparameters<RealT>::Determine_GSVD_Hyperparameters(void) const
  {
    std::cout << "Need to implement Determine_GSVD_Hyperparameters" << std::endl;
  }

}

#endif
