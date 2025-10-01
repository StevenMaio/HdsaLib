#ifndef PDEOPT_STATE_COST_SHALLOW_ICE_HPP
#define PDEOPT_STATE_COST_SHALLOW_ICE_HPP

#include "pde_shallow_ice.hpp"

template <class Real>
class State_Cost_shallow_ice : public Objective_SimOpt_TS<Real>
{
private:
  std::vector<std::vector<Real>> data_; // data[k][i] is the kth time instance in the ith spatial node of the true state
  std::vector<Real> data_weight_;
  std::vector<int> data_weight_id_;
  ROL::Ptr<Tpetra::MultiVector<>> w_ptr_;
  ROL::Ptr<ROL::Vector<Real>> wp_;
  ROL::Ptr<Tpetra::MultiVector<>> noisy_data_ptr_;
  ROL::Ptr<ROL::Vector<Real>> noisy_data_;

  void Apply_Weight_Vector(ROL::Vector<Real> &vec)
  {
    wp_->zero();
    for (unsigned k = 0; k < data_weight_.size(); k++)
    {
      w_ptr_->replaceGlobalValue(data_weight_id_[k], 0, data_weight_[k]);
    }
    vec.applyBinary(ROL::Elementwise::Multiply<Real>(), *wp_);
  }

  void Apply_Sqrt_Weight_Vector(ROL::Vector<Real> &vec)
  {
    wp_->zero();
    for (unsigned k = 0; k < data_weight_.size(); k++)
    {
      w_ptr_->replaceGlobalValue(data_weight_id_[k], 0, std::sqrt(data_weight_[k]));
    }
    vec.applyBinary(ROL::Elementwise::Multiply<Real>(), *wp_);
  }

public:
  State_Cost_shallow_ice(std::vector<ROL::TimeStamp<Real>> &timeStamp, std::vector<std::vector<Real>> &data,
                         const std::vector<Real> &data_weight, const std::vector<int> &data_weight_id, ROL::Ptr<ROL::Vector<Real>> &up)
      : Objective_SimOpt_TS<Real>(timeStamp), data_(data), data_weight_(data_weight), data_weight_id_(data_weight_id)
  {
    wp_ = up->clone();
    w_ptr_ = dynamic_cast<ROL::TpetraMultiVector<Real> &>(*wp_).getVector();
    noisy_data_ = up->clone();
    noisy_data_ptr_ = dynamic_cast<ROL::TpetraMultiVector<Real> &>(*noisy_data_).getVector();
    wp_->zero();
    noisy_data_->zero();

    for (unsigned k = 0; k < data_weight_.size(); k++)
    {
      int index = data_weight_id_[k];
      noisy_data_ptr_->replaceGlobalValue(index, 0, data[0][index]);
    }
  }

  Real value(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol)
  {
    ROL::Ptr<ROL::Vector<Real>> vec = u.clone();
    vec->set(u);
    vec->axpy(-1.0, *noisy_data_);
    Apply_Sqrt_Weight_Vector(*vec);
    return 0.5 * (vec->dot(*vec));
  }

  void gradient_1(ROL::Vector<Real> &g, const ROL::Vector<Real> &u,
                  const ROL::Vector<Real> &z, Real &tol)
  {
    g.set(u);
    g.scale(1.0);
    g.axpy(-1.0, *noisy_data_);
    Apply_Weight_Vector(g);
  }

  void gradient_2(ROL::Vector<Real> &g, const ROL::Vector<Real> &u,
                  const ROL::Vector<Real> &z, Real &tol)
  {
    g.zero();
  }

  void hessVec_11(ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
                  const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol)
  {
    hv.set(v);
    Apply_Weight_Vector(hv);
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

  void update(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &t, bool flag = true, int iter = -1)
  {
    Objective_SimOpt_TS<Real>::update(u, z, t, flag, iter);
    Update_Noisy_Data();
  }

  void Update_Noisy_Data()
  {
    noisy_data_->zero();
    for (unsigned k = 0; k < data_weight_.size(); k++)
    {
      int index = data_weight_id_[k];
      noisy_data_ptr_->replaceGlobalValue(index, 0, data_[Objective_SimOpt_TS<Real>::current_TS_][index]);
    }
  }

}; // State_Cost

template <class Real>
class QoI_H1_shallow_ice : public QoI<Real>
{

  ROL::Ptr<FE<Real>> fe_;
  ROL::Ptr<FieldHelper<Real>> fieldHelper_;

public:
  QoI_H1_shallow_ice(const ROL::Ptr<FE<Real>> &fe, const ROL::Ptr<FieldHelper<Real>> &fieldHelper) : fe_(fe), fieldHelper_(fieldHelper) {}

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real>> &val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    // Get relevant dimensions
    int c = z_coeff->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    int d = fe_->gradN()->dimension(3);

    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);

    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);

    ROL::Ptr<Intrepid::FieldContainer<Real>> grad_eval;
    grad_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    fe_->evaluateGradient(grad_eval, Z[0]);

    ROL::Ptr<Intrepid::FieldContainer<Real>> gradNorm_eval;
    gradNorm_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);

    for (int i = 0; i < c; i++)
    {
      for (int j = 0; j < p; j++)
      {
        (*gradNorm_eval)(i, j) = (*grad_eval)(i, j, 0) * (*grad_eval)(i, j, 0) + (*grad_eval)(i, j, 1) * (*grad_eval)(i, j, 1);
      }
    }

    // Create array of 1's to integrate against
    ROL::Ptr<Intrepid::FieldContainer<Real>> const_weights = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    const_weights->initialize(1.0);
    fe_->computeIntegral(val, gradNorm_eval, const_weights);
    Intrepid::RealSpaceTools<Real>::scale(*val, static_cast<Real>(0.5));

    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real>> &grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> QoI_H1_shallow_ice::gradient_1 is zero.");
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real>> &grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    // Get relevant dimensions
    int c = z_coeff->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    int f = fe_->N()->dimension(1);
    int d = fe_->gradN()->dimension(3);

    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> G(3);
    for (int i = 0; i < 3; ++i)
    {
      G[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    }

    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);

    // Build local gradient of state tracking term
    ROL::Ptr<Intrepid::FieldContainer<Real>> grad_eval;
    grad_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    fe_->evaluateGradient(grad_eval, Z[0]);

    ROL::Ptr<Intrepid::FieldContainer<Real>> diffx;
    diffx = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real>> diffy;
    diffy = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);

    for (int i = 0; i < c; i++)
    {
      for (int j = 0; j < p; j++)
      {
        (*diffx)(i, j) = (*grad_eval)(i, j, 0);
        (*diffy)(i, j) = (*grad_eval)(i, j, 1);
      }
    }

    Intrepid::FunctionSpaceTools::integrate<Real>(*G[0],
                                                  *diffx,
                                                  *(fe_->DNDdetJ(0)),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::FunctionSpaceTools::integrate<Real>(*G[0],
                                                  *diffy,
                                                  *(fe_->DNDdetJ(1)),
                                                  Intrepid::COMP_CPP, true);
    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real>> &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> QoI_H1_shallow_ice::HessVec_11 is zero.");
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real>> &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> QoI_H1_shallow_ice::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>> &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> QoI_H1_shallow_ice::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>> &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    int c = v_coeff->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    int f = fe_->N()->dimension(1);
    int d = fe_->gradN()->dimension(3);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> H(3);
    for (int i = 0; i < 3; ++i)
    {
      H[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    }
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);

    ROL::Ptr<Intrepid::FieldContainer<Real>> gradV_eval;
    gradV_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    fe_->evaluateGradient(gradV_eval, V[0]);

    ROL::Ptr<Intrepid::FieldContainer<Real>> diffx;
    diffx = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real>> diffy;
    diffy = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);

    for (int i = 0; i < c; i++)
    {
      for (int j = 0; j < p; j++)
      {
        (*diffx)(i, j) = (*gradV_eval)(i, j, 0);
        (*diffy)(i, j) = (*gradV_eval)(i, j, 1);
      }
    }

    Intrepid::FunctionSpaceTools::integrate<Real>(*H[0],
                                                  *diffx,
                                                  *(fe_->DNDdetJ(0)),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::FunctionSpaceTools::integrate<Real>(*H[0],
                                                  *diffy,
                                                  *(fe_->DNDdetJ(1)),
                                                  Intrepid::COMP_CPP, true);

    fieldHelper_->combineFieldCoeff(hess, H);
  }

}; // QoI_H1_shallow_ice

template <class Real>
class QoI_L2_shallow_ice : public QoI<Real>
{

  ROL::Ptr<FE<Real>> fe_;
  ROL::Ptr<FieldHelper<Real>> fieldHelper_;

public:
  QoI_L2_shallow_ice(const ROL::Ptr<FE<Real>> &fe, const ROL::Ptr<FieldHelper<Real>> &fieldHelper) : fe_(fe), fieldHelper_(fieldHelper) {}

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real>> &val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    // Get relevant dimensions
    int c = z_coeff->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);

    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    ROL::Ptr<Intrepid::FieldContainer<Real>> z =
        ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(z, Z[0]);

    // Build local state tracking term
    fe_->computeIntegral(val, z, z);
    Intrepid::RealSpaceTools<Real>::scale(*val, static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real>> &grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> QoI_L2_shallow_ice::gradient_1 is zero.");
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real>> &grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    // Get relevant dimensions
    int c = z_coeff->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    int f = fe_->N()->dimension(1);
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> G(3);
    for (int i = 0; i < 3; ++i)
    {
      G[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    }

    // Build local gradient of state tracking term
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    ROL::Ptr<Intrepid::FieldContainer<Real>> z =
        ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(z, Z[0]);
    Intrepid::FunctionSpaceTools::integrate<Real>(*G[0],
                                                  *z,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real>> &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> QoI_L2_shallow_ice::HessVec_11 is zero.");
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real>> &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> QoI_L2_shallow_ice::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>> &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> QoI_L2_shallow_ice::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>> &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    int c = v_coeff->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    int f = fe_->N()->dimension(1);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> H(3);
    for (int i = 0; i < 3; ++i)
    {
      H[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    }
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    ROL::Ptr<Intrepid::FieldContainer<Real>> valV_eval =
        ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(valV_eval, V[0]);
    Intrepid::FunctionSpaceTools::integrate<Real>(*H[0],
                                                  *valV_eval,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    fieldHelper_->combineFieldCoeff(hess, H);
  }

}; // QoI_L2_shallow_ice

#endif
