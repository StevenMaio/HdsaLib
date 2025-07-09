#ifndef HDSA_MD_OUU_U_PRIOR_INTERFACE_HPP
#define HDSA_MD_OUU_U_PRIOR_INTERFACE_HPP

namespace HDSA
{

  template <class RealT>
  class MD_OUU_u_Prior_Interface : public HDSA::MD_u_Prior_Interface<RealT>
  {

  private:
    const HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT>> us_prior_interface_;
    const HDSA::Ptr<HDSA::Dense_Matrix<RealT>> K_;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> C_;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> Cinv_;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> R_;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> Rinv_;
    RealT scaling_;
    int ens_size_;

  public:
    MD_OUU_u_Prior_Interface(const HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT>> &us_prior_interface, const HDSA::Ptr<HDSA::Dense_Matrix<RealT>> &K) : us_prior_interface_(us_prior_interface), K_(K)
    {
      ens_size_ = K->numRows();

      C_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(ens_size_, ens_size_);
      for (int i = 0; i < ens_size_; i++)
      {
        for (int j = 0; j < ens_size_; j++)
        {
          RealT val = 0.0;
          if (i == j)
          {
            val = -(*K)(i, i);
            for (int s = 0; s < ens_size_; s++)
            {
              val += 2.0 * (*K)(i, s);
            }
          }
          else
          {
            val = -2.0 * (*K)(i, j);
          }
          C_->Replace_Element(i, j, val);
        }
      }

      R_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(ens_size_, ens_size_);
      HDSA::Linear_Algebra::Cholesky_Factorization<RealT>(*C_, *R_);

      Rinv_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(ens_size_, ens_size_);
      for (int j = 0; j < ens_size_; j++)
      {
        HDSA::Ptr<HDSA::Dense_Matrix<RealT>> x = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(ens_size_, 1);
        HDSA::Ptr<HDSA::Dense_Matrix<RealT>> b = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(ens_size_, 1);
        b->Replace_Element(j, 0, 1.0);
        HDSA::Linear_Algebra::Upper_Tri_Solve<RealT>(*x, *b, *R_);
        for (int i = 0; i < ens_size_; i++)
        {
          Rinv_->Replace_Element(i, j, (*x)(i, 0));
        }
      }

      Cinv_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(ens_size_, ens_size_);
      Rinv_->Multiply(*Cinv_, *Rinv_, false, true);

      RealT tmp = 0.0;
      for (int s = 0; s < ens_size_; s++)
      {
        tmp += (*Cinv_)(s, s);
      }
      scaling_ = static_cast<RealT>(ens_size_) / tmp;

      for (int i = 0; i < ens_size_; i++)
      {
        for (int j = 0; j < ens_size_; j++)
        {
          RealT val = (1.0 / scaling_) * (*C_)(i, j);
          C_->Replace_Element(i, j, val);
          val = scaling_ * (*Cinv_)(i, j);
          Cinv_->Replace_Element(i, j, val);
          val = (1.0 / std::sqrt(scaling_)) * (*R_)(i, j);
          R_->Replace_Element(i, j, val);
          val = std::sqrt(scaling_) * (*Rinv_)(i, j);
          Rinv_->Replace_Element(i, j, val);
        }
      }
    }

    virtual ~MD_OUU_u_Prior_Interface()
    {
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Base class pure virtual function implementations
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void Apply_M_u(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &u_in) const
    {
      u_out.zeros();
      HDSA::Ptr<HDSA::Vector<RealT>> u_tmp = u_out.clone();

      const HDSA::Ensemble_Vector<RealT> u_in_ens = dynamic_cast<const HDSA::Ensemble_Vector<RealT> &>(u_in);
      HDSA::Ensemble_Vector<RealT> u_out_ens = dynamic_cast<HDSA::Ensemble_Vector<RealT> &>(u_out);
      HDSA::Ensemble_Vector<RealT> u_tmp_ens = dynamic_cast<HDSA::Ensemble_Vector<RealT> &>(*u_tmp);

      for (int s = 0; s < ens_size_; s++)
      {
        us_prior_interface_->Apply_M_u(*u_tmp_ens[s], *u_in_ens[s]);
      }

      for (int s = 0; s < ens_size_; s++)
      {
        for (int i = 0; i < ens_size_; i++)
        {
          u_out_ens[s]->axpy((*C_)(s, i), *u_tmp_ens[i]);
        }
      }
    }

    void Apply_W_u_Inverse(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &u_in) const
    {
      u_out.zeros();
      HDSA::Ptr<HDSA::Vector<RealT>> u_tmp = u_out.clone();

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
          u_out_ens[s]->axpy((*Cinv_)(s, i), *u_tmp_ens[i]);
        }
      }
    }

    void Apply_W_u_Plus_scalar_M_u_Inverse(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &u_in, const RealT &scalar) const
    {
      u_out.zeros();
      HDSA::Ptr<HDSA::Vector<RealT>> u_tmp = u_out.clone();

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
          u_out_ens[s]->axpy((*Cinv_)(s, i), *u_tmp_ens[i]);
        }
      }
    }

    void Sample_with_Covariance_W_u_Inverse(HDSA::MultiVector<RealT> &samples) const
    {
      samples.zeros();
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
            vec_ens[s]->axpy((*Rinv_)(s, j), *(*ind_samples)[j]);
          }
        }
      }
    }

    void Sample_with_Covariance_W_u_Plus_scalar_M_u_Inverse(HDSA::MultiVector<RealT> &samples, const RealT &scalar) const
    {
      samples.zeros();
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
            vec_ens[s]->axpy((*Rinv_)(s, j), *(*ind_samples)[j]);
          }
        }
      }
    }
  };

}

#endif
