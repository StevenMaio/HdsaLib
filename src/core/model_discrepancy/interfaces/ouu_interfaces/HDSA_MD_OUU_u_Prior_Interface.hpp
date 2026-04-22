/***********************************************************************
 HdsaLib - A library for Hyper-differential Sensitivity Analysis
 
 Questions? Contact Joseph Hart (joshart@sandia.gov)
************************************************************************/

#ifndef HDSA_MD_OUU_U_PRIOR_INTERFACE_HPP
#define HDSA_MD_OUU_U_PRIOR_INTERFACE_HPP

#include "HDSA_MD_u_Prior_Interface.hpp"
#include "HDSA_Dense_Matrix.hpp"
#include "HDSA_Linear_Algebra.hpp"
#include "HDSA_Ensemble_Vector.hpp"
#include "HDSA_MD_OUU_Ensemble_Weighting_Matrix.hpp"

namespace HDSA
{

  template <class RealT>
  class MD_OUU_u_Prior_Interface : public HDSA::MD_u_Prior_Interface<RealT>
  {

  private:
    const HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT>> us_prior_interface_;
    const HDSA::Ptr<HDSA::MD_OUU_Ensemble_Weighting_Matrix<RealT>> ensemble_weighting_;
    int ens_size_;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> W_s_;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> W_s_inv_;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> R_inv_;

  public:
    MD_OUU_u_Prior_Interface(const HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT>> &us_prior_interface, const HDSA::Ptr<HDSA::MD_OUU_Ensemble_Weighting_Matrix<RealT>> &ensemble_weighting) : us_prior_interface_(us_prior_interface), ensemble_weighting_(ensemble_weighting)
    {
      ens_size_ = ensemble_weighting_->Get_ens_size();
      W_s_ = ensemble_weighting_->Get_W_s();
      W_s_inv_ = ensemble_weighting_->Get_W_s_inv();
      R_inv_ = ensemble_weighting_->Get_R_inv();
    }

    virtual ~MD_OUU_u_Prior_Interface()
    {
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Base class pure virtual function implementations
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void Apply_M_u(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &u_in) const
    {
      u_out.Zeros();

      if (const HDSA::Ensemble_Vector<RealT> *u_in_ens = dynamic_cast<const HDSA::Ensemble_Vector<RealT> *>(&u_in))
      {
        HDSA::Ptr<HDSA::Vector<RealT>> u_tmp = u_out.Clone();

        HDSA::Ensemble_Vector<RealT> u_out_ens = dynamic_cast<HDSA::Ensemble_Vector<RealT> &>(u_out);
        HDSA::Ensemble_Vector<RealT> u_tmp_ens = dynamic_cast<HDSA::Ensemble_Vector<RealT> &>(*u_tmp);

        for (int s = 0; s < ens_size_; s++)
        {
          us_prior_interface_->Apply_M_u(*u_tmp_ens[s], *(*u_in_ens)[s]);
        }

        for (int s = 0; s < ens_size_; s++)
        {
          for (int i = 0; i < ens_size_; i++)
          {
            u_out_ens[s]->Scaled_Plus((*W_s_)(s, i), *u_tmp_ens[i]);
          }
        }
      }
      else
      {
        us_prior_interface_->Apply_M_u(u_out, u_in);
      }
    }

    void Apply_W_u_Inverse(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &u_in) const
    {
      u_out.Zeros();
      HDSA::Ptr<HDSA::Vector<RealT>> u_tmp = u_out.Clone();

      const HDSA::Ensemble_Vector<RealT> u_in_ens = dynamic_cast<const HDSA::Ensemble_Vector<RealT> &>(u_in);
      HDSA::Ensemble_Vector<RealT> u_out_ens = dynamic_cast<HDSA::Ensemble_Vector<RealT> &>(u_out);
      HDSA::Ensemble_Vector<RealT> u_tmp_ens = dynamic_cast<HDSA::Ensemble_Vector<RealT> &>(*u_tmp);

      for (int s = 0; s < ens_size_; s++)
      {
        us_prior_interface_->Apply_W_u_Inverse(*u_tmp_ens[s], *u_in_ens[s]);
      }

      for (int s = 0; s < ens_size_; s++)
      {
        for (int i = 0; i < ens_size_; i++)
        {
          u_out_ens[s]->Scaled_Plus((*W_s_inv_)(s, i), *u_tmp_ens[i]);
        }
      }
    }

    void Apply_W_u_Plus_scalar_M_u_Inverse(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &u_in, const RealT &scalar) const
    {
      u_out.Zeros();
      HDSA::Ptr<HDSA::Vector<RealT>> u_tmp = u_out.Clone();

      const HDSA::Ensemble_Vector<RealT> u_in_ens = dynamic_cast<const HDSA::Ensemble_Vector<RealT> &>(u_in);
      HDSA::Ensemble_Vector<RealT> u_out_ens = dynamic_cast<HDSA::Ensemble_Vector<RealT> &>(u_out);
      HDSA::Ensemble_Vector<RealT> u_tmp_ens = dynamic_cast<HDSA::Ensemble_Vector<RealT> &>(*u_tmp);

      for (int s = 0; s < ens_size_; s++)
      {
        us_prior_interface_->Apply_W_u_Plus_scalar_M_u_Inverse(*u_tmp_ens[s], *u_in_ens[s], scalar);
      }

      for (int s = 0; s < ens_size_; s++)
      {
        for (int i = 0; i < ens_size_; i++)
        {
          u_out_ens[s]->Scaled_Plus((*W_s_inv_)(s, i), *u_tmp_ens[i]);
        }
      }
    }

    void Sample_with_Covariance_W_u_Inverse(HDSA::MultiVector<RealT> &samples) const
    {
      samples.Zeros();
      int M = samples.Number_of_Vectors();
      for (int i = 0; i < M; i++)
      {
        HDSA::Ensemble_Vector<RealT> vec_ens = dynamic_cast<HDSA::Ensemble_Vector<RealT> &>(*samples[i]);
        HDSA::Ptr<HDSA::MultiVector<RealT>> ind_samples = HDSA::makePtr<HDSA::MultiVector<RealT>>(ens_size_, *vec_ens[0]);
        us_prior_interface_->Sample_with_Covariance_W_u_Inverse(*ind_samples);
        for (int s = 0; s < ens_size_; s++)
        {
          for (int j = s; j < ens_size_; j++)
          {
            vec_ens[s]->Scaled_Plus((*R_inv_)(s, j), *(*ind_samples)[j]);
          }
        }
      }
    }

    void Sample_with_Covariance_W_u_Plus_scalar_M_u_Inverse(HDSA::MultiVector<RealT> &samples, const RealT &scalar) const
    {
      samples.Zeros();
      int M = samples.Number_of_Vectors();
      for (int i = 0; i < M; i++)
      {
        HDSA::Ensemble_Vector<RealT> vec_ens = dynamic_cast<HDSA::Ensemble_Vector<RealT> &>(*samples[i]);
        HDSA::Ptr<HDSA::MultiVector<RealT>> ind_samples = HDSA::makePtr<HDSA::MultiVector<RealT>>(ens_size_, *vec_ens[0]);
        us_prior_interface_->Sample_with_Covariance_W_u_Plus_scalar_M_u_Inverse(*ind_samples, scalar);
        for (int s = 0; s < ens_size_; s++)
        {
          for (int j = s; j < ens_size_; j++)
          {
            vec_ens[s]->Scaled_Plus((*R_inv_)(s, j), *(*ind_samples)[j]);
          }
        }
      }
    }
  };

}

#endif
