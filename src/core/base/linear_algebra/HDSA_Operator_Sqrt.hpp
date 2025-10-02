#ifndef HDSA_OPERATOR_SQRT_HPP
#define HDSA_OPERATOR_SQRT_HPP

#include "HDSA_Dense_Matrix.hpp"
#include "HDSA_Linear_Algebra.hpp"

namespace HDSA
{

  template <class RealT>
  class Operator_Sqrt
  {
  private:
    int max_iter_;
    RealT tol_;

  public:
    Operator_Sqrt(int max_iter = 1000, RealT tol = 1.e-8) : max_iter_(max_iter), tol_(tol)
    {
    }

    virtual ~Operator_Sqrt()
    {
    }

    virtual void Apply_Operator(HDSA::Vector<RealT> &vec_out, const HDSA::Vector<RealT> &vec_in) const = 0;

    void Apply_Operator_Sqrt(HDSA::Vector<RealT> &vec_out, const HDSA::Vector<RealT> &vec_in)
    {
      RealT Norm_in = vec_in.Norm();

      std::vector<HDSA::Ptr<HDSA::Vector<RealT>>> V;
      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> T = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(max_iter_, max_iter_);
      RealT alpha = 0.0;
      RealT beta = 0.0;
      std::vector<RealT> rel_res;
      std::vector<RealT> ykp;
      std::vector<RealT> yk;

      HDSA::Ptr<HDSA::Vector<RealT>> v_j = vec_in.Clone();

      v_j->Set(vec_in);
      v_j->Scale(1.0 / Norm_in);

      for (int j = 0; j < max_iter_; j++)
      {
        HDSA::Ptr<HDSA::Vector<RealT>> v_push = vec_in.Clone();
        v_push->Set(*v_j);
        V.push_back(v_push);

        HDSA::Ptr<HDSA::Vector<RealT>> w_j = vec_in.Clone();
        Apply_Operator(*w_j, *v_j);
        alpha = w_j->Dot(*v_j);
        w_j->axpy(-alpha, *v_j);
        if (j > 0)
        {
          w_j->axpy(-beta, *V[j - 1]);
        }
        beta = w_j->Norm();

        v_j->Set(*w_j);
        v_j->Scale(1.0 / beta);

        for (int k = 0; k < j; k++)
        {
          RealT val = v_j->Dot(*V[j]);
          v_j->axpy(-val, *V[j]);
        }
        RealT val = v_j->Norm();
        v_j->Scale(1.0 / val);

        for (int k = 0; k < j; k++)
        {
          RealT val = v_j->Dot(*V[j]);
          v_j->axpy(-val, *V[j]);
        }
        val = v_j->Norm();
        v_j->Scale(1.0 / val);

        T->Set_Entry(j, j, alpha);
        T->Set_Entry(j, j + 1, beta);
        T->Set_Entry(j + 1, j, beta);

        HDSA::Ptr<HDSA::Dense_Matrix<RealT>> Tk = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(j + 1, j + 1);
        for (int i = 0; i < j + 1; i++)
        {
          for (int k = 0; k < j + 1; k++)
          {
            Tk->Set_Entry(i, k, (*T)(i, k));
          }
        }

        HDSA::Ptr<HDSA::Dense_Matrix<RealT>> V = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(j + 1, j + 1);
        HDSA::Ptr<HDSA::Dense_Matrix<RealT>> S = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(j + 1, 1);
        HDSA::Linear_Algebra::Symmetric_Eig_Decomposition<RealT>(*Tk, *V, *S);
        HDSA::Ptr<HDSA::Dense_Matrix<RealT>> e = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(j + 1, 1);
        e->Set_Entry(0, 0, Norm_in);
        HDSA::Ptr<HDSA::Dense_Matrix<RealT>> Ve = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(j + 1, 1);
        V->Multiply(*Ve, *e, true, false);
        yk.clear();
        yk.resize(j + 1);
        for (int i = 0; i < j + 1; i++)
        {
          for (int k = 0; k < j + 1; k++)
          {
            yk[i] += (*Ve)(k, 0) * std::sqrt((*S)(k, 0)) * (*V)(i, k);
          }
        }

        if (j == 0)
        {
          rel_res.push_back(std::sqrt(2));
        }
        else
        {
          RealT num = 0.0;
          RealT denom = 0.0;
          for (int k = 0; k < j; k++)
          {
            num += std::pow(yk[k] - ykp[k], 2.0);
            denom += std::pow(yk[k], 2.0);
          }
          num += std::pow(yk[j], 2.0);
          denom += std::pow(yk[j], 2.0);
          rel_res.push_back(std::sqrt(num / denom));
        }

        if (rel_res[j] < tol_)
        {
          break;
        }
        else
        {
          ykp.clear();
          ykp.resize(yk.size());
          for (int k = 0; k < ykp.size(); k++)
          {
            ykp[k] = yk[k];
          }
        }
      }

      vec_out.Zeros();
      for (int i = 0; i < yk.size(); i++)
      {
        vec_out.axpy(yk[i], *V[i]);
      }
    }
  };

}

#endif
