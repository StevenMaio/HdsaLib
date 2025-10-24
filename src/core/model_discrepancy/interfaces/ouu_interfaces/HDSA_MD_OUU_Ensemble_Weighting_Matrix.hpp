#ifndef HDSA_MD_OUU_ENSEMBLE_WEIGHTING_HPP
#define HDSA_MD_OUU_ENSEMBLE_WEIGHTING_HPP

#include "HDSA_MD_Data_Interface.hpp"
#include "HDSA_MD_u_Prior_Interface.hpp"
#include "HDSA_Dense_Matrix.hpp"
#include "HDSA_Linear_Algebra.hpp"

namespace HDSA
{

  template <class RealT>
  class MD_OUU_Ensemble_Weighting_Matrix
  {

  private:
    const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> data_interface_;
    const HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT>> u_prior_interface_;
    int ens_size_;
    RealT max_marginal_var_percent_, min_cond_variance_percent_;
    bool assume_independent_;

    RealT reg_opt_, reg_opt_min_cond_var_percent_;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> W_s_, W_s_inv_, R_inv_, C_;

  public:
    MD_OUU_Ensemble_Weighting_Matrix(const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> &data_interface, const HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT>> &u_prior_interface, int ens_size, RealT max_marginal_var_percent = 1.0, RealT min_cond_variance_percent = 0.1, bool assume_independent = false)
        : data_interface_(data_interface), u_prior_interface_(u_prior_interface), ens_size_(ens_size), max_marginal_var_percent_(max_marginal_var_percent), min_cond_variance_percent_(min_cond_variance_percent), assume_independent_(assume_independent)
    {
      if (assume_independent_)
      {
        W_s_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(ens_size_, ens_size_);
        W_s_inv_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(ens_size_, ens_size_);
        R_inv_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(ens_size_, ens_size_);
        C_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(ens_size_, ens_size_);
        for (int s = 0; s < ens_size_; s++)
        {
          W_s_->Set_Entry(s, s, 1.0);
          W_s_inv_->Set_Entry(s, s, 1.0);
          R_inv_->Set_Entry(s, s, 1.0);
          C_->Set_Entry(s, s, 1.0);
        }
      }
      else
      {
        Compute_Matrices();
      }
    }

    void Compute_Matrices(void)
    {
      HDSA::Ptr<const HDSA::MultiVector<RealT>> D = data_interface_->Get_D();
      int N = D->Number_of_Vectors();

      // Assemble Eta
      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> Eta = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(ens_size_, ens_size_);
      for (int ell = 0; ell < N; ell++)
      {
        HDSA::Ptr<HDSA::Vector<RealT>> d = (*D)[ell];
        HDSA::Ensemble_Vector<RealT> d_ens = dynamic_cast<const HDSA::Ensemble_Vector<RealT> &>(*d);
        HDSA::Ptr<HDSA::Vector<RealT>> vec1 = d_ens[0]->Clone();
        HDSA::Ptr<HDSA::Vector<RealT>> vec2 = d_ens[0]->Clone();
        for (int s = 0; s < ens_size_; s++)
        {
          for (int k = 0; k < ens_size_; k++)
          {
            vec1->Zeros();
            vec2->Zeros();
            if (s == k)
            {
              u_prior_interface_->Apply_M_u(*vec1, *d_ens[s]);
              RealT val = (*Eta)(s, k) + vec1->Dot(*d_ens[s]) + 1.e-16;
              Eta->Set_Entry(s, k, val);
            }
            else
            {
              vec1->Set(*d_ens[s]);
              vec1->Scaled_Plus(-1.0, *d_ens[k]);
              u_prior_interface_->Apply_M_u(*vec2, *vec1);
              RealT val = (*Eta)(s, k) + vec1->Dot(*vec2) + 1.e-16;
              Eta->Set_Entry(s, k, val);
            }
          }
        }
      }
      Eta->Scale(1.0 / static_cast<RealT>(N));

      // Assemble A
      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> A = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(ens_size_, ens_size_);
      for (int s = 0; s < ens_size_; s++)
      {
        for (int k = 0; k < ens_size_; k++)
        {
          if (s == k)
          {
            A->Set_Entry(s, k, 1.0);
          }
          else
          {
            RealT val = std::sqrt((*Eta)(s, s)) * std::sqrt((*Eta)(k, k)) / (*Eta)(s, k);
            A->Set_Entry(s, k, val);
          }
        }
      }

      // Find highly correlated pair
      int row = 0;
      RealT best_val = std::numeric_limits<double>::infinity();
      for (int s = 0; s < ens_size_; s++)
      {
        RealT row_sum = 0.0;
        for (int k = 0; k < ens_size_; k++)
        {
          row_sum += (*Eta)(s, k);
        }
        if (row_sum < best_val)
        {
          row = s;
          best_val = row_sum;
        }
      }

      int col = 0;
      best_val = std::numeric_limits<double>::infinity();
      for (int k = 0; k < ens_size_; k++)
      {
        if (((*Eta)(row, k) < best_val) && (k != row))
        {
          col = k;
          best_val = (*Eta)(row, k);
        }
      }

      // Search over 0.1, 0.2, ..., 0.9
      RealT best_coeff = 0.0;
      RealT best_obj = std::numeric_limits<double>::infinity();
      for (int i = 1; i < 10; i++)
      {
        RealT coeff = static_cast<RealT>(i) / 10.0;
        RealT obj_val = Min_Cond_Var_Obj(*A, coeff, row, col);
        if (obj_val < best_obj)
        {
          best_coeff = coeff;
          best_obj = obj_val;
        }
      }

      // Search over [best_coeff-0.1,best_coeff+0.1]
      RealT lb = best_coeff - 0.1;
      for (int i = 0; i < 21; i++)
      {
        RealT coeff = lb + static_cast<RealT>(i) / 100.0;
        RealT obj_val = Min_Cond_Var_Obj(*A, coeff, row, col);
        if (obj_val < best_obj)
        {
          best_coeff = coeff;
          best_obj = obj_val;
        }
      }

      // Search over [best_coeff-0.01,best_coeff+0.01]
      lb = best_coeff - 0.01;
      for (int i = 0; i < 21; i++)
      {
        RealT coeff = lb + static_cast<RealT>(i) / 1000.0;
        RealT obj_val = Min_Cond_Var_Obj(*A, coeff, row, col);
        if (obj_val < best_obj)
        {
          best_coeff = coeff;
          best_obj = obj_val;
        }
      }

      reg_opt_ = best_coeff;
      W_s_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(ens_size_, ens_size_);
      W_s_inv_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(ens_size_, ens_size_);
      R_inv_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(ens_size_, ens_size_);
      C_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(ens_size_, ens_size_);
      Assemble_Matrices(*W_s_, *W_s_inv_, *R_inv_, *C_, *A, reg_opt_);

      reg_opt_min_cond_var_percent_ = ((*W_s_inv_)(row, row) - std::pow((*W_s_inv_)(row, col), 2.0) / (*W_s_inv_)(col, col)) / (*W_s_inv_)(row, row);

      if (std::abs(reg_opt_min_cond_var_percent_ - min_cond_variance_percent_) / min_cond_variance_percent_ > 1.e-2)
      {
        std::cout << "Error in MD_OUU_Ensemble_Weighting_Matrix: Failed to achieve conditional variance target" << std::endl;
      }
    }

    void Assemble_Matrices(HDSA::Dense_Matrix<RealT> &Ws, HDSA::Dense_Matrix<RealT> &Wsinv, HDSA::Dense_Matrix<RealT> &Rinv, HDSA::Dense_Matrix<RealT> &C, const HDSA::Dense_Matrix<RealT> &A, RealT reg_coeff)
    {
      // Assemble C
      for (int s = 0; s < ens_size_; s++)
      {
        for (int k = 0; k < ens_size_; k++)
        {
          if (s == k)
          {
            C.Set_Entry(s, k, 1.0);
          }
          else
          {
            C.Set_Entry(s, k, reg_coeff * A(s, k));
          }
        }
      }

      // Assemble Ws
      for (int s = 0; s < ens_size_; s++)
      {
        RealT running_val = 0.0;
        for (int k = 0; k < ens_size_; k++)
        {
          if (s != k)
          {
            Ws.Set_Entry(s, k, -2.0 * C(s, k));
            running_val += C(s, k);
          }
        }
        Ws.Set_Entry(s, s, C(s, s) + 2.0 * running_val);
      }

      // Factorize Ws
      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> R = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(ens_size_, ens_size_);
      HDSA::Linear_Algebra::Cholesky_Factorization<RealT>(Ws, *R);

      // Invert R
      for (int k = 0; k < ens_size_; k++)
      {
        HDSA::Ptr<HDSA::Dense_Matrix<RealT>> x = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(ens_size_, 1);
        HDSA::Ptr<HDSA::Dense_Matrix<RealT>> b = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(ens_size_, 1);
        b->Set_Entry(k, 0, 1.0);
        HDSA::Linear_Algebra::Upper_Tri_Solve<RealT>(*x, *b, *R);
        for (int s = 0; s < ens_size_; s++)
        {
          Rinv.Set_Entry(s, k, (*x)(s, 0));
        }
      }

      // Compute Ws inverse
      Rinv.Multiply(Wsinv, Rinv, false, true);

      // Rescale
      RealT diag_max = 0.0;
      for (int s = 0; s < ens_size_; s++)
      {
        diag_max = std::max(diag_max, Wsinv(s, s));
      }
      RealT scaling = max_marginal_var_percent_ / diag_max;
      Wsinv.Scale(scaling);
      Ws.Scale(1.0 / scaling);
      Rinv.Scale(std::sqrt(scaling));
    }

    RealT Min_Cond_Var_Obj(const HDSA::Dense_Matrix<RealT> &A, RealT reg_coeff, int i, int j)
    {
      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> Ws = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(ens_size_, ens_size_);
      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> Wsinv = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(ens_size_, ens_size_);
      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> Rinv = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(ens_size_, ens_size_);
      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> C = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(ens_size_, ens_size_);
      Assemble_Matrices(*Ws, *Wsinv, *Rinv, *C, A, reg_coeff);
      RealT val = (*Wsinv)(i, i) - std::pow((*Wsinv)(i, j), 2.0) / (*Wsinv)(j, j) - min_cond_variance_percent_ * (*Wsinv)(i, i);
      val = std::pow(val, 2.0);
      return val;
    }

    int Get_ens_size(void) const
    {
      return ens_size_;
    }

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> Get_W_s(void) const
    {
      return W_s_;
    }

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> Get_W_s_inv(void) const
    {
      return W_s_inv_;
    }

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> Get_R_inv(void) const
    {
      return R_inv_;
    }
  };

}

#endif
