/***********************************************************************
 HdsaLib - A library for Hyper-differential Sensitivity Analysis
 
 Questions? Contact Joseph Hart (joshart@sandia.gov)
************************************************************************/

#ifndef SIMOPT_TEST_PROBLEM_HPP
#define SIMOPT_TEST_PROBLEM_HPP

#include "ROL_Algorithm.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_CompositeStep.hpp"
#include "ROL_ConstraintStatusTest.hpp"
#include "ROL_Types.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_LAPACK.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_StdVector.hpp"
#include "ROL_Vector_SimOpt.hpp"
#include "ROL_Constraint_SimOpt.hpp"
#include "ROL_Objective_SimOpt.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <ctime>

template <class Real>
class Constraint_SimOptTestProb : public ROL::Constraint_SimOpt<Real>
{

private:
  int m_;

public:
  Constraint_SimOptTestProb(int m)
  {
    m_ = m;
  }

  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol)
  {
    ROL::Ptr<std::vector<Real>> cp =
        dynamic_cast<ROL::StdVector<Real> &>(c).getVector();
    ROL::Ptr<const std::vector<Real>> up =
        dynamic_cast<const ROL::StdVector<Real> &>(u).getVector();
    ROL::Ptr<const std::vector<Real>> zp =
        dynamic_cast<const ROL::StdVector<Real> &>(z).getVector();
    for (int i = 0; i < m_; i++)
    {
      (*cp)[i] = (*up)[i] - std::pow((*zp)[i], 3.0);
    }
  }

  void applyJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                       const ROL::Vector<Real> &z, Real &tol)
  {
    ROL::Ptr<std::vector<Real>> jvp =
        dynamic_cast<ROL::StdVector<Real> &>(jv).getVector();
    ROL::Ptr<const std::vector<Real>> vp =
        dynamic_cast<const ROL::StdVector<Real> &>(v).getVector();
    ROL::Ptr<const std::vector<Real>> up =
        dynamic_cast<const ROL::StdVector<Real> &>(u).getVector();
    ROL::Ptr<const std::vector<Real>> zp =
        dynamic_cast<const ROL::StdVector<Real> &>(z).getVector();
    for (int i = 0; i < m_; i++)
    {
      (*jvp)[i] = (*vp)[i];
    }
  }

  void applyInverseJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol)
  {
    ROL::Ptr<std::vector<Real>> jvp =
        dynamic_cast<ROL::StdVector<Real> &>(jv).getVector();
    ROL::Ptr<const std::vector<Real>> vp =
        dynamic_cast<const ROL::StdVector<Real> &>(v).getVector();
    ROL::Ptr<const std::vector<Real>> up =
        dynamic_cast<const ROL::StdVector<Real> &>(u).getVector();
    ROL::Ptr<const std::vector<Real>> zp =
        dynamic_cast<const ROL::StdVector<Real> &>(z).getVector();
    for (int i = 0; i < m_; i++)
    {
      (*jvp)[i] = (*vp)[i];
    }
  }

  void applyAdjointJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol)
  {
    ROL::Ptr<std::vector<Real>> jvp =
        dynamic_cast<ROL::StdVector<Real> &>(jv).getVector();
    ROL::Ptr<const std::vector<Real>> vp =
        dynamic_cast<const ROL::StdVector<Real> &>(v).getVector();
    ROL::Ptr<const std::vector<Real>> up =
        dynamic_cast<const ROL::StdVector<Real> &>(u).getVector();
    ROL::Ptr<const std::vector<Real>> zp =
        dynamic_cast<const ROL::StdVector<Real> &>(z).getVector();
    for (int i = 0; i < m_; i++)
    {
      (*jvp)[i] = (*vp)[i];
    }
  }

  void applyInverseAdjointJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol)
  {
    ROL::Ptr<std::vector<Real>> jvp =
        dynamic_cast<ROL::StdVector<Real> &>(jv).getVector();
    ROL::Ptr<const std::vector<Real>> vp =
        dynamic_cast<const ROL::StdVector<Real> &>(v).getVector();
    ROL::Ptr<const std::vector<Real>> up =
        dynamic_cast<const ROL::StdVector<Real> &>(u).getVector();
    ROL::Ptr<const std::vector<Real>> zp =
        dynamic_cast<const ROL::StdVector<Real> &>(z).getVector();
    for (int i = 0; i < m_; i++)
    {
      (*jvp)[i] = (*vp)[i];
    }
  }

  void applyJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol)
  {
    ROL::Ptr<std::vector<Real>> jvp =
        dynamic_cast<ROL::StdVector<Real> &>(jv).getVector();
    ROL::Ptr<const std::vector<Real>> vp =
        dynamic_cast<const ROL::StdVector<Real> &>(v).getVector();
    ROL::Ptr<const std::vector<Real>> up =
        dynamic_cast<const ROL::StdVector<Real> &>(u).getVector();
    ROL::Ptr<const std::vector<Real>> zp =
        dynamic_cast<const ROL::StdVector<Real> &>(z).getVector();
    jv.zero();
    for (int i = 0; i < m_; i++)
    {
      (*jvp)[i] = -3.0 * std::pow((*zp)[i], 2.0) * (*vp)[i];
    }
  }

  void applyAdjointJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol)
  {
    ROL::Ptr<std::vector<Real>> jvp =
        dynamic_cast<ROL::StdVector<Real> &>(jv).getVector();
    ROL::Ptr<const std::vector<Real>> vp =
        dynamic_cast<const ROL::StdVector<Real> &>(v).getVector();
    ROL::Ptr<const std::vector<Real>> up =
        dynamic_cast<const ROL::StdVector<Real> &>(u).getVector();
    ROL::Ptr<const std::vector<Real>> zp =
        dynamic_cast<const ROL::StdVector<Real> &>(z).getVector();
    jv.zero();
    for (int i = 0; i < m_; i++)
    {
      (*jvp)[i] = -3.0 * std::pow((*zp)[i], 2.0) * (*vp)[i];
    }
  }

  void applyAdjointHessian_11(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol)
  {
    ahwv.zero();
  }

  void applyAdjointHessian_12(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol)
  {
    ahwv.zero();
  }

  void applyAdjointHessian_21(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol)
  {
    ahwv.zero();
  }

  void applyAdjointHessian_22(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol)
  {
    ROL::Ptr<std::vector<Real>> ahwvp =
        dynamic_cast<ROL::StdVector<Real> &>(ahwv).getVector();
    ROL::Ptr<const std::vector<Real>> wp =
        dynamic_cast<const ROL::StdVector<Real> &>(w).getVector();
    ROL::Ptr<const std::vector<Real>> vp =
        dynamic_cast<const ROL::StdVector<Real> &>(v).getVector();
    ROL::Ptr<const std::vector<Real>> up =
        dynamic_cast<const ROL::StdVector<Real> &>(u).getVector();
    ROL::Ptr<const std::vector<Real>> zp =
        dynamic_cast<const ROL::StdVector<Real> &>(z).getVector();
    ahwv.zero();
    for (int i = 0; i < m_; i++)
    {
      (*ahwvp)[i] = -6.0 * (*wp)[i] * (*zp)[i] * (*vp)[i];
    }
  }
};

template <class Real>
class Objective_SimOptTestProb : public ROL::Objective_SimOpt<Real>
{
private:
  int m_;
  std::vector<Real> target_;

public:
  Objective_SimOptTestProb(int m)
  {
    m_ = m;
    target_.resize(m_);
    for (int i = 0; i < m_; i++)
    {
      target_[i] = std::pow(1.0 + (double)i / double(m_ - 1), 3.0);
    }
  }

  Real value(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol)
  {
    ROL::Ptr<const std::vector<Real>> up =
        dynamic_cast<const ROL::StdVector<Real> &>(u).getVector();
    ROL::Ptr<const std::vector<Real>> zp =
        dynamic_cast<const ROL::StdVector<Real> &>(z).getVector();
    Real val = 0.0;
    for (int i = 0; i < m_; i++)
    {
      val += std::pow((*up)[i] - target_[i], 2.0);
    }
    return 0.5 * val;
  }

  void gradient_1(ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol)
  {
    // Unwrap g
    ROL::Ptr<std::vector<Real>> gup =
        dynamic_cast<ROL::StdVector<Real> &>(g).getVector();
    // Unwrap x
    ROL::Ptr<const std::vector<Real>> up =
        dynamic_cast<const ROL::StdVector<Real> &>(u).getVector();
    ROL::Ptr<const std::vector<Real>> zp =
        dynamic_cast<const ROL::StdVector<Real> &>(z).getVector();
    // COMPUTE GRADIENT WRT U
    for (int i = 0; i < m_; i++)
    {
      (*gup)[i] = (*up)[i] - target_[i];
    }
  }

  void gradient_2(ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol)
  {
    g.zero();
  }

  void hessVec_11(ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
                  const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol)
  {
    ROL::Ptr<std::vector<Real>> hvup =
        dynamic_cast<ROL::StdVector<Real> &>(hv).getVector();
    // Unwrap v
    ROL::Ptr<const std::vector<Real>> vup =
        dynamic_cast<const ROL::StdVector<Real> &>(v).getVector();
    // COMPUTE GRADIENT WRT U
    for (int i = 0; i < m_; i++)
    {
      (*hvup)[i] = (*vup)[i];
    }
  }

  void hessVec_12(ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
                  const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol)
  {
    hv.zero();
  }

  void hessVec_21(ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
                  const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol)
  {
    hv.zero();
  }

  void hessVec_22(ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
                  const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol)
  {
    hv.zero();
  }
};

#endif
