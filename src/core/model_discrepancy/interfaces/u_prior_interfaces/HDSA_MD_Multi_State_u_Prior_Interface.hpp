/***********************************************************************
 HdsaLib - A library for Hyper-differential Sensitivity Analysis
 
 Questions? Contact Joseph Hart (joshart@sandia.gov)
************************************************************************/

#ifndef HDSA_MD_MULTI_STATE_U_PRIOR_INTERFACE_HPP
#define HDSA_MD_MULTI_STATE_U_PRIOR_INTERFACE_HPP

namespace HDSA
{

  template <class RealT>
  class MD_Multi_State_u_Prior_Interface : public HDSA::MD_u_Prior_Interface<RealT>
  {

  private:
    const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> data_interface_;
    const std::vector<HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT>>> u_prior_interface_std_;
    int num_components_;

  public:
    MD_Multi_State_u_Prior_Interface(const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> &data_interface, const std::vector<HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT>>> &u_prior_interface_std) : data_interface_(data_interface), u_prior_interface_std_(u_prior_interface_std)
    {
      num_components_ = u_prior_interface_std_.size();
    }

    virtual ~MD_Multi_State_u_Prior_Interface()
    {
    }

    void Apply_M_u(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &u_in) const override
    {
      for (int k = 0; k < num_components_; k++)
      {
        HDSA::Ptr<const HDSA::Vector<RealT>> u_in_k = data_interface_->Extract_State_Component(u_in, k, true);
        HDSA::Ptr<HDSA::Vector<RealT>> u_out_k = u_in_k->Clone();
        u_prior_interface_std_[k]->Apply_M_u(*u_out_k, *u_in_k);
        data_interface_->Set_State_Component(u_out, *u_out_k, k, true);
      }
    }

    void Apply_W_u_Inverse(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &u_in) const override
    {
      for (int k = 0; k < num_components_; k++)
      {
        HDSA::Ptr<const HDSA::Vector<RealT>> u_in_k = data_interface_->Extract_State_Component(u_in, k, true);
        HDSA::Ptr<HDSA::Vector<RealT>> u_out_k = u_in_k->Clone();
        u_prior_interface_std_[k]->Apply_W_u_Inverse(*u_out_k, *u_in_k);
        data_interface_->Set_State_Component(u_out, *u_out_k, k, true);
      }
    }

    void Apply_W_u_Plus_scalar_M_u_Inverse(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &u_in, const RealT &scalar) const override
    {
      for (int k = 0; k < num_components_; k++)
      {
        HDSA::Ptr<const HDSA::Vector<RealT>> u_in_k = data_interface_->Extract_State_Component(u_in, k, true);
        HDSA::Ptr<HDSA::Vector<RealT>> u_out_k = u_in_k->Clone();
        u_prior_interface_std_[k]->Apply_W_u_Plus_scalar_M_u_Inverse(*u_out_k, *u_in_k, scalar);
        data_interface_->Set_State_Component(u_out, *u_out_k, k, true);
      }
    }

    void Sample_with_Covariance_W_u_Inverse(HDSA::MultiVector<RealT> &samples) const override
    {
      int num_vecs = samples.Number_of_Vectors();
      for (int k = 0; k < num_components_; k++)
      {
        HDSA::Ptr<const HDSA::Vector<RealT>> vec = data_interface_->Extract_State_Component(*(samples[0]), k, true);
        HDSA::Ptr<HDSA::MultiVector<RealT>> samples_k = HDSA::makePtr<HDSA::MultiVector<RealT>>(num_vecs, *vec);
        u_prior_interface_std_[k]->Sample_with_Covariance_W_u_Inverse(*samples_k);
        for (int j = 0; j < num_vecs; j++)
        {
          data_interface_->Set_State_Component(*(samples[j]), *((*samples_k)[j]), k, true);
        }
      }
    }

    void Sample_with_Covariance_W_u_Plus_scalar_M_u_Inverse(HDSA::MultiVector<RealT> &samples, const RealT &scalar) const override
    {
      int num_vecs = samples.Number_of_Vectors();
      for (int k = 0; k < num_components_; k++)
      {
        HDSA::Ptr<const HDSA::Vector<RealT>> vec = data_interface_->Extract_State_Component(*(samples[0]), k, true);
        HDSA::Ptr<HDSA::MultiVector<RealT>> samples_k = HDSA::makePtr<HDSA::MultiVector<RealT>>(num_vecs, *vec);
        u_prior_interface_std_[k]->Sample_with_Covariance_W_u_Plus_scalar_M_u_Inverse(*samples_k, scalar);
        for (int j = 0; j < num_vecs; j++)
        {
          data_interface_->Set_State_Component(*(samples[j]), *((*samples_k)[j]), k, true);
        }
      }
    }

    void Precompute_W_u_Plus_scalar_M_u_Data(RealT &scalar) override
    {
      for (int k = 0; k < num_components_; k++)
      {
        u_prior_interface_std_[k]->Precompute_W_u_Plus_scalar_M_u_Data(scalar);
      }
    }

  };

}

#endif
