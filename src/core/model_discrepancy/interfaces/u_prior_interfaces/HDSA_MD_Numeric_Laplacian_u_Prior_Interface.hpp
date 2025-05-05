#ifndef HDSA_MD_NUMERIC_LAPLACIAN_U_PRIOR_INTERFACE_HPP
#define HDSA_MD_NUMERIC_LAPLACIAN_U_PRIOR_INTERFACE_HPP

namespace HDSA
{

  template <class RealT>
  class MD_Numeric_Laplacian_u_Prior_Interface : public HDSA::MD_Elliptic_u_Prior_Interface<RealT>
  {

  private:
    const HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> S_;
    const HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> M_;
    const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> data_interface_;
    const HDSA::Ptr<HDSA::MD_u_Hyperparameter_Interface<RealT>> u_hyperparam_interface_;
    HDSA::Ptr<HDSA::MD_Determine_u_Hyperparameters<RealT>> determine_u_hyperparams_;
    RealT beta_u_;
    HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> E_u_;
    HDSA::Ptr<HDSA::Sparse_Matrix_Solver<RealT>> E_u_solver_;

  public:
    void Apply_E_u_Inverse(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &u_in) const
    {
      E_u_solver_->Apply_A_Inverse(u_out, u_in);
    }

    void Apply_E_u_Inverse_Transpose(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &u_in) const
    {
      E_u_solver_->Apply_A_Inverse(u_out, u_in);
    }

    void Apply_M_u(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &u_in) const
    {
      M_->Apply(u_out, u_in);
    }

    MD_Numeric_Laplacian_u_Prior_Interface(const HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> &S, const HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> &M, const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> &data_interface, const HDSA::Ptr<HDSA::MD_u_Hyperparameter_Interface<RealT>> &u_hyperparam_interface) : HDSA::MD_Elliptic_u_Prior_Interface<RealT>(u_hyperparam_interface->Get_alpha_u()), S_(S), M_(M), data_interface_(data_interface), u_hyperparam_interface_(u_hyperparam_interface)
    {
      Auxillary_Constructor();
    }

    MD_Numeric_Laplacian_u_Prior_Interface(const HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> &S, const HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> &M, const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> &data_interface, const HDSA::Ptr<HDSA::MD_u_Hyperparameter_Interface<RealT>> &u_hyperparam_interface, int seed) : HDSA::MD_Elliptic_u_Prior_Interface<RealT>(u_hyperparam_interface->Get_alpha_u(), seed), S_(S), M_(M), data_interface_(data_interface), u_hyperparam_interface_(u_hyperparam_interface)
    {
      Auxillary_Constructor();
    }

    MD_Numeric_Laplacian_u_Prior_Interface(const HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> &S, const HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> &M, const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> &data_interface, const HDSA::Ptr<HDSA::MD_u_Hyperparameter_Interface<RealT>> &u_hyperparam_interface, const HDSA::Ptr<HDSA::Random_Number_Generator<RealT>> &random_number_generator) : HDSA::MD_Elliptic_u_Prior_Interface<RealT>(u_hyperparam_interface->Get_alpha_u(), random_number_generator), S_(S), M_(M), data_interface_(data_interface), u_hyperparam_interface_(u_hyperparam_interface)
    {
      Auxillary_Constructor();
    }

    void Auxillary_Constructor()
    {
      E_u_ = M_->clone();
      determine_u_hyperparams_ = HDSA::makePtr<HDSA::MD_Determine_u_Hyperparameters<RealT>>(data_interface_, u_hyperparam_interface_);

      if (u_hyperparam_interface_->Get_beta_u() == 0.0)
      {
        determine_u_hyperparams_->Determine_beta_u();
      }
      Set_beta_u(u_hyperparam_interface_->Get_beta_u());

      if (u_hyperparam_interface_->Get_gsvd_num_sing_vals() == 0)
      {
        determine_u_hyperparams_->Determine_GSVD_Hyperparameters();
      }
      HDSA::Ptr<const HDSA::Vector<RealT>> u_opt = data_interface_->get_u_opt();
      HDSA::Ptr<const HDSA::Vector<RealT>> u_vec;
      if (const Transient_Vector<RealT>* u_trans = dynamic_cast<const Transient_Vector<RealT>*>(&(*u_opt)))
      {
        u_vec = data_interface_->Extract_State_Component(*(*u_trans)[0],u_hyperparam_interface_-> Get_Componenet_ID());
      }
      else
      {
        u_vec = data_interface_->Extract_State_Component(*u_opt,u_hyperparam_interface_-> Get_Componenet_ID());
      }
      this->Compute_E_u_Inverse_GSVD(u_hyperparam_interface_->Get_gsvd_num_sing_vals(), u_hyperparam_interface_->Get_gsvd_oversampling(), u_hyperparam_interface_->Get_gsvd_num_subspace_iter(), *u_vec);

      if (!u_hyperparam_interface_->Is_Transient())
      {
        if (u_hyperparam_interface_->Get_alpha_u() == 0.0)
        {
          determine_u_hyperparams_->Determine_alpha_u(this);
        }
        this->Set_alpha_u(u_hyperparam_interface_->Get_alpha_u());
      }
      determine_u_hyperparams_->Determine_alpha_d();
    }

    virtual ~MD_Numeric_Laplacian_u_Prior_Interface()
    {
    }

    void Set_beta_u(RealT beta_u_new)
    {
      E_u_->set(*M_);
      E_u_->axpy(beta_u_new, *S_);
      beta_u_ = beta_u_new;
      E_u_solver_ = HDSA::makePtr<HDSA::Sparse_Matrix_Solver<RealT>>(E_u_);
    }
  };
}

#endif
