/***********************************************************************
 HdsaLib - A library for Hyper-differential Sensitivity Analysis
 
 Questions? Contact Joseph Hart (joshart@sandia.gov)
************************************************************************/

#ifndef HDSA_MD_OUU_HYPERPARAMETER_DATA_INTERFACE_HPP
#define HDSA_MD_OUU_HYPERPARAMETER_DATA_INTERFACE_HPP

#include "HDSA_MD_Data_Interface.hpp"

namespace HDSA
{

  template <class RealT>
  class MD_OUU_Hyperparameter_Data_Interface : public HDSA::MD_Data_Interface<RealT>
  {

  private:
    int ens_size_;
    HDSA::Ptr<MD_OUU_Data_Interface<RealT>> md_ouu_data_interface_;

  public:
    MD_OUU_Hyperparameter_Data_Interface(HDSA::Ptr<MD_OUU_Data_Interface<RealT>> &md_ouu_data_interface) : ens_size_(md_ouu_data_interface->Get_Ensemble_Size()), md_ouu_data_interface_(md_ouu_data_interface)
    {
    }

    virtual ~MD_OUU_Hyperparameter_Data_Interface()
    {
    }

    HDSA::Ptr<HDSA::Vector<RealT>> Load_Optimal_u(void) const
    {
      HDSA::Ptr<HDSA::Vector<RealT>> u_opt = md_ouu_data_interface_->Load_Optimal_us(0)->Clone();
      u_opt->Zeros();
      for (int s = 0; s < ens_size_; s++)
      {
        HDSA::Ptr<HDSA::Vector<RealT>> us_opt = md_ouu_data_interface_->Load_Optimal_us(s);
        u_opt->Plus(*us_opt);
      }
      u_opt->Scale(1.0 / static_cast<RealT>(ens_size_));
      return u_opt;
    }

    HDSA::Ptr<HDSA::Vector<RealT>> Load_Optimal_z(void) const
    {
      HDSA::Ptr<HDSA::Vector<RealT>> z_opt = md_ouu_data_interface_->Load_Optimal_z();
      return z_opt;
    }

    HDSA::Ptr<HDSA::MultiVector<RealT>> Load_Z_Data(void) const
    {
      HDSA::Ptr<HDSA::MultiVector<RealT>> Z = md_ouu_data_interface_->Load_Z_Data();
      return Z;
    }

    HDSA::Ptr<HDSA::MultiVector<RealT>> Load_D_Data(void) const
    {
      HDSA::Ptr<HDSA::MultiVector<RealT>> D = md_ouu_data_interface_->Load_Ds_Data(0)->Clone();
      for (int s = 0; s < ens_size_; s++)
      {
        D->Plus(*md_ouu_data_interface_->Load_Ds_Data(s));
      }
      D->Scale(1.0 / static_cast<RealT>(ens_size_));
      return D;
    }
  };

}

#endif
