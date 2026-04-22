/***********************************************************************
 HdsaLib - A library for Hyper-differential Sensitivity Analysis
 
 Questions? Contact Joseph Hart (joshart@sandia.gov)
************************************************************************/

#ifndef ELLIPTIC_Z_PRIOR_INTERFACE_SIMOPT_TEST_PROBLEM_HPP
#define ELLIPTIC_Z_PRIOR_INTERFACE_SIMOPT_TEST_PROBLEM_HPP

#include "HDSA_MD_Elliptic_z_Prior_Interface.hpp"

template <class RealT>
class Elliptic_z_Prior_Interface_SimOptTestProb : public HDSA::MD_Elliptic_z_Prior_Interface<RealT>
{

private:
  int m_;                                    // Mesh resolution
  HDSA::Ptr<HDSA::Dense_Matrix<RealT>> S_;   // Stiffness matrix
  HDSA::Ptr<HDSA::Dense_Matrix<RealT>> M_;   // Mass matrix
  HDSA::Ptr<HDSA::Dense_Matrix<RealT>> E_z_; // State elliptic operator
  const HDSA::Ptr<HDSA::Random_Number_Generator<RealT>> random_number_generator_;

public:
  Elliptic_z_Prior_Interface_SimOptTestProb(RealT &alpha_z, RealT &beta_z, int &m, const HDSA::Ptr<HDSA::Random_Number_Generator<RealT>> &random_number_generator)
      : HDSA::MD_Elliptic_z_Prior_Interface<RealT>(alpha_z), random_number_generator_(random_number_generator)
  {
    m_ = m;

    RealT h = 1.0 / static_cast<RealT>(m_ - 1);

    S_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, m_);
    M_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, m_);

    S_->Set_Entry(0, 0, 1.0 / h);
    S_->Set_Entry(0, 1, -1.0 / h);
    for (int i = 1; i < m_ - 1; i++)
    {
      S_->Set_Entry(i, i, 2.0 / h);
      S_->Set_Entry(i, i - 1, -1.0 / h);
      S_->Set_Entry(i, i + 1, -1.0 / h);
    }
    S_->Set_Entry(m_ - 1, m_ - 2, -1.0 / h);
    S_->Set_Entry(m_ - 1, m_ - 1, 1.0 / h);

    M_->Set_Entry(0, 0, (1.0 / 3.0) * h);
    M_->Set_Entry(0, 1, (1.0 / 6.0) * h);
    for (int i = 1; i < m_ - 1; i++)
    {
      M_->Set_Entry(i, i, (2.0 / 3.0) * h);
      M_->Set_Entry(i, i - 1, (1.0 / 6.0) * h);
      M_->Set_Entry(i, i + 1, (1.0 / 6.0) * h);
    }
    M_->Set_Entry(m_ - 1, m_ - 2, (1.0 / 6.0) * h);
    M_->Set_Entry(m_ - 1, m_ - 1, (1.0 / 3.0) * h);

    E_z_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, m_);
    for (int i = 0; i < m_; i++)
    {
      for (int j = 0; j < m_; j++)
      {
        RealT val = beta_z * (*S_)(i, j) + (*M_)(i, j);
        E_z_->Set_Entry(i, j, val);
      }
    }
  }

  virtual ~Elliptic_z_Prior_Interface_SimOptTestProb()
  {
  }

  void Apply_E_z_Inverse(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in) const
  {
    HDSA::ROL_Vector<RealT> &z_out_rol = dynamic_cast<HDSA::ROL_Vector<RealT> &>(z_out);
    ROL::Ptr<std::vector<RealT>> z_out_std = dynamic_cast<ROL::StdVector<RealT> &>(*z_out_rol.rol_vec).getVector();

    const HDSA::ROL_Vector<RealT> &z_in_rol = dynamic_cast<const HDSA::ROL_Vector<RealT> &>(z_in);
    ROL::Ptr<const std::vector<RealT>> z_in_std = dynamic_cast<const ROL::StdVector<RealT> &>(*z_in_rol.rol_vec).getVector();

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> b = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    for (int k = 0; k < m_; k++)
    {
      b->Set_Entry(k, 0, (*z_in_std)[k]);
    }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> x = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*E_z_, *x, *b);
    for (int k = 0; k < m_; k++)
    {
      (*z_out_std)[k] = (*x)(k, 0);
    }
  }

  void Apply_E_z_Inverse_Transpose(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in) const
  {
    Apply_E_z_Inverse(z_out, z_in);
  }

  void Apply_M_z(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in) const
  {
    HDSA::ROL_Vector<RealT> &z_out_rol = dynamic_cast<HDSA::ROL_Vector<RealT> &>(z_out);
    ROL::Ptr<std::vector<RealT>> z_out_std = dynamic_cast<ROL::StdVector<RealT> &>(*z_out_rol.rol_vec).getVector();

    const HDSA::ROL_Vector<RealT> &z_in_rol = dynamic_cast<const HDSA::ROL_Vector<RealT> &>(z_in);
    ROL::Ptr<const std::vector<RealT>> z_in_std = dynamic_cast<const ROL::StdVector<RealT> &>(*z_in_rol.rol_vec).getVector();

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> b = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    for (int k = 0; k < m_; k++)
    {
      b->Set_Entry(k, 0, (*z_in_std)[k]);
    }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> x = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    M_->Multiply(*x, *b);
    for (int k = 0; k < m_; k++)
    {
      (*z_out_std)[k] = (*x)(k, 0);
    }
  }

  void Sample_with_Covariance_W_z_Inverse(HDSA::MultiVector<RealT> &samples) const
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> R = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, m_);
    HDSA::Linear_Algebra::Cholesky_Factorization(*M_, *R);

    int num_samples = samples.Number_of_Vectors();
    for (int i = 0; i < num_samples; i++)
    {

      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> b = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
      for (int k = 0; k < m_; k++)
      {
        b->Set_Entry(k, 0, random_number_generator_->Generate_Standard_Normal_Sample());
      }
      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> x = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);

      // Should be E_z^{-1}*R^T*b
      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> tmp = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
      R->Multiply(*tmp, *b, true);
      HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*E_z_, *x, *tmp);

      HDSA::ROL_Vector<RealT> &z_i = dynamic_cast<HDSA::ROL_Vector<RealT> &>(*samples[i]);
      ROL::Ptr<std::vector<RealT>> z_i_std = dynamic_cast<ROL::StdVector<RealT> &>(*z_i.rol_vec).getVector();
      for (int k = 0; k < m_; k++)
      {
        (*z_i_std)[k] = (*x)(k, 0);
      }
    }
  }

  void Apply_E_z(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in) const
  {
    HDSA::ROL_Vector<RealT> &z_out_rol = dynamic_cast<HDSA::ROL_Vector<RealT> &>(z_out);
    ROL::Ptr<std::vector<RealT>> z_out_std = dynamic_cast<ROL::StdVector<RealT> &>(*z_out_rol.rol_vec).getVector();

    const HDSA::ROL_Vector<RealT> &z_in_rol = dynamic_cast<const HDSA::ROL_Vector<RealT> &>(z_in);
    ROL::Ptr<const std::vector<RealT>> z_in_std = dynamic_cast<const ROL::StdVector<RealT> &>(*z_in_rol.rol_vec).getVector();

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> b = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    for (int k = 0; k < m_; k++)
    {
      b->Set_Entry(k, 0, (*z_in_std)[k]);
    }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> x = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    E_z_->Multiply(*x, *b);
    for (int k = 0; k < m_; k++)
    {
      (*z_out_std)[k] = (*x)(k, 0);
    }
  }

  void Apply_E_z_Transpose(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in) const
  {
    Apply_E_z(z_out, z_in);
  }

  void Apply_M_z_Inverse(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in) const
  {
    HDSA::ROL_Vector<RealT> &z_out_rol = dynamic_cast<HDSA::ROL_Vector<RealT> &>(z_out);
    ROL::Ptr<std::vector<RealT>> z_out_std = dynamic_cast<ROL::StdVector<RealT> &>(*z_out_rol.rol_vec).getVector();

    const HDSA::ROL_Vector<RealT> &z_in_rol = dynamic_cast<const HDSA::ROL_Vector<RealT> &>(z_in);
    ROL::Ptr<const std::vector<RealT>> z_in_std = dynamic_cast<const ROL::StdVector<RealT> &>(*z_in_rol.rol_vec).getVector();

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> b = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    for (int k = 0; k < m_; k++)
    {
      b->Set_Entry(k, 0, (*z_in_std)[k]);
    }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> x = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*M_, *x, *b);
    for (int k = 0; k < m_; k++)
    {
      (*z_out_std)[k] = (*x)(k, 0);
    }
  }
};

#endif
