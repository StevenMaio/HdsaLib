#ifndef HDSA_MD_Z_PRIOR_INTERFACE_SYNTHETIC_TEST_HESSIAN_GEVEP_HPP
#define HDSA_MD_Z_PRIOR_INTERFACE_SYNTHETIC_TEST_HESSIAN_GEVEP_HPP

template <class RealT>
class MD_z_Prior_Interface_synthetic_test_hessian_gevp : public HDSA::MD_z_Prior_Interface<RealT>
{

private:
  int m_;                                    // Mesh resolution
  HDSA::Ptr<HDSA::Dense_Matrix<RealT>> x_;   // Mesh nodes on [0,1]
  HDSA::Ptr<HDSA::Dense_Matrix<RealT>> S_;   // Stiffness matrix
  HDSA::Ptr<HDSA::Dense_Matrix<RealT>> M_;   // Mass matrix
  HDSA::Ptr<HDSA::Dense_Matrix<RealT>> W_z_; // Control weighting matrix
  const HDSA::Ptr<HDSA::Random_Number_Generator<RealT>> random_number_generator_;

public:
  MD_z_Prior_Interface_synthetic_test_hessian_gevp(const HDSA::Ptr<HDSA::Random_Number_Generator<RealT>> &random_number_generator) : random_number_generator_(random_number_generator)
  {
    m_ = 51;
    RealT h = 1.0 / static_cast<RealT>(m_ - 1);
    x_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    for (int k = 0; k < m_; k++)
    {
      x_->Replace_Element(k, 0, static_cast<RealT>(k) / static_cast<RealT>(m_ - 1));
    }

    S_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, m_);
    M_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, m_);

    S_->Replace_Element(0, 0, 1.0 / h);
    S_->Replace_Element(0, 1, -1.0 / h);
    for (int i = 1; i < m_ - 1; i++)
    {
      S_->Replace_Element(i, i, 2.0 / h);
      S_->Replace_Element(i, i - 1, -1.0 / h);
      S_->Replace_Element(i, i + 1, -1.0 / h);
    }
    S_->Replace_Element(m_ - 1, m_ - 2, -1.0 / h);
    S_->Replace_Element(m_ - 1, m_ - 1, 1.0 / h);

    M_->Replace_Element(0, 0, (1.0 / 3.0) * h);
    M_->Replace_Element(0, 1, (1.0 / 6.0) * h);
    for (int i = 1; i < m_ - 1; i++)
    {
      M_->Replace_Element(i, i, (2.0 / 3.0) * h);
      M_->Replace_Element(i, i - 1, (1.0 / 6.0) * h);
      M_->Replace_Element(i, i + 1, (1.0 / 6.0) * h);
    }
    M_->Replace_Element(m_ - 1, m_ - 2, (1.0 / 6.0) * h);
    M_->Replace_Element(m_ - 1, m_ - 1, (1.0 / 3.0) * h);

    W_z_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, m_);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> I = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, m_);
    for (int k = 0; k < m_; k++)
    {
      I->Replace_Element(k, k, 1.0);
    }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> Minv = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, m_);
    HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*M_, *Minv, *I);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> E_z = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, m_);
    for (int i = 0; i < m_; i++)
    {
      for (int j = 0; j < m_; j++)
      {
        RealT val = (2.0) * ((1.e-2) * (*S_)(i, j) + (*M_)(i, j));
        E_z->Replace_Element(i, j, val);
      }
    }

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> tmp = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, m_);
    Minv->Multiply(*tmp, *E_z);
    E_z->Multiply(*W_z_, *tmp);
  }

  virtual ~MD_z_Prior_Interface_synthetic_test_hessian_gevp()
  {
  }

  void Apply_M_z(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in) const
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> b = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    const Std_Vector<RealT> &z_in_std = dynamic_cast<const Std_Vector<RealT> &>(z_in);
    Std_Vector<RealT> &z_out_std = dynamic_cast<Std_Vector<RealT> &>(z_out);
    for (int k = 0; k < m_; k++)
    {
      b->Replace_Element(k, 0, z_in_std(k));
    }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> x = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);

    M_->Multiply(*x, *b);
    for (int k = 0; k < m_; k++)
    {
      z_out_std.Replace_Element(k, (*x)(k, 0));
    }
  }

  void Apply_W_z_Inverse(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in) const
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> b = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    const Std_Vector<RealT> &z_in_std = dynamic_cast<const Std_Vector<RealT> &>(z_in);
    Std_Vector<RealT> &z_out_std = dynamic_cast<Std_Vector<RealT> &>(z_out);
    for (int k = 0; k < m_; k++)
    {
      b->Replace_Element(k, 0, z_in_std(k));
    }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> x = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*W_z_, *x, *b);
    for (int k = 0; k < m_; k++)
    {
      z_out_std.Replace_Element(k, (*x)(k, 0));
    }
  }

  // Compute samples from a mean zero Gaussian with covariance W_z^{-1}
  virtual void Sample_with_Covariance_W_z_Inverse(HDSA::MultiVector<RealT> &samples) const
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> R = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, m_);
    HDSA::Linear_Algebra::Cholesky_Factorization(*W_z_, *R);

    int num_samples = samples.Number_of_Vectors();
    for (int i = 0; i < num_samples; i++)
    {
      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> b = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
      for (int k = 0; k < m_; k++)
      {
        b->Replace_Element(k, 0, random_number_generator_->Generate_Standard_Normal_Sample());
      }
      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> x = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
      HDSA::Linear_Algebra::Upper_Tri_Solve<RealT>(*x, *b, *R);
      Std_Vector<RealT> &vec_out_std = dynamic_cast<Std_Vector<RealT> &>(*samples[i]);
      for (int k = 0; k < m_; k++)
      {
        vec_out_std.Replace_Element(k, (*x)(k, 0));
      }
    }
  }

  void Apply_W_z(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in) const
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> b = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    const Std_Vector<RealT> &z_in_std = dynamic_cast<const Std_Vector<RealT> &>(z_in);
    Std_Vector<RealT> &z_out_std = dynamic_cast<Std_Vector<RealT> &>(z_out);
    for (int k = 0; k < m_; k++)
    {
      b->Replace_Element(k, 0, z_in_std(k));
    }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> x = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);

    W_z_->Multiply(*x, *b);
    for (int k = 0; k < m_; k++)
    {
      z_out_std.Replace_Element(k, (*x)(k, 0));
    }
  }
};

#endif
