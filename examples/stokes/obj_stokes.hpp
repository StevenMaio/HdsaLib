#ifndef PDEOPT_QOI_STOKES_HPP
#define PDEOPT_QOI_STOKES_HPP

#include "../../../PDE-OPT/TOOLS/qoi.hpp"
#include "pde_stokes.hpp"

template <class Real>
class QoI_Vertical_Velocity_Stokes : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real> > feVel_;
  const ROL::Ptr<FE<Real> > fePrs_;
  const ROL::Ptr<FieldHelper<Real> > fieldHelper_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > weight_;
  Real channelW_;
  Real stepW_;

  const Real eps_;

  Real weightFunc(const std::vector<Real> & x) const {
    Real in_domain = 0.0;
    if( (x[0] >= stepW_-eps_) && (x[0] <= channelW_-eps_) && (x[1] >= 0.2-eps_) && (x[1] <= 0.8+eps_) )
      {
	in_domain = 1.0;
      }
    return in_domain;
  }

public:
  QoI_Vertical_Velocity_Stokes(const ROL::Ptr<FE<Real> > &feVel,
				     const ROL::Ptr<FE<Real> > &fePrs,
				     const ROL::Ptr<FieldHelper<Real> > &fieldHelper,
				     Teuchos::ParameterList &parlist)
    : feVel_(feVel), fePrs_(fePrs), fieldHelper_(fieldHelper), eps_(std::sqrt(ROL::ROL_EPSILON<Real>())) {
    channelW_ = parlist.sublist("Geometry").get(    "Channel width", 8.0);
    stepW_    = parlist.sublist("Geometry").get(       "Step width", 1.0);
    int c = feVel_->cubPts()->dimension(0);
    int p = feVel_->cubPts()->dimension(1);
    int d = feVel_->cubPts()->dimension(2);
    std::vector<Real> pt(d);
    weight_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          pt[k] = (*feVel_->cubPts())(i,j,k);
        }
        (*weight_)(i,j) = weightFunc(pt);
      }
    }
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real> > & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = u_coeff->dimension(0);
    int p = feVel_->cubPts()->dimension(1);
    int d = feVel_->cubPts()->dimension(2);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Evaluate on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real> > valUY_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    feVel_->evaluateValue(valUY_eval, U[1]);
    // Multiply by weight
    ROL::Ptr<Intrepid::FieldContainer<Real> > weighted_valUY_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*weighted_valUY_eval,
                                                               *weight_,
                                                               *valUY_eval);
    // Compute L2 norm squared
    feVel_->computeIntegral(val,valUY_eval,weighted_valUY_eval,false);
    Intrepid::RealSpaceTools<Real>::scale(*val,static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = u_coeff->dimension(0);
    int fv = feVel_->N()->dimension(1);
    int fp = fePrs_->N()->dimension(1);
    int p = feVel_->cubPts()->dimension(1);
    int d = feVel_->cubPts()->dimension(2);

    Intrepid::FieldContainer<Real> velUX_grad(c, fv);
    Intrepid::FieldContainer<Real> velUY_grad(c, fv);
    Intrepid::FieldContainer<Real> presU_grad(c, fp);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > G;
    G.resize(fieldHelper_->numFields());
    G[0] = ROL::makePtrFromRef(velUX_grad);
    G[1] = ROL::makePtrFromRef(velUY_grad);
    G[2] = ROL::makePtrFromRef(presU_grad);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Evaluate on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real> > valUY_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    feVel_->evaluateValue(valUY_eval, U[1]);
    // Multiply by weight
    ROL::Ptr<Intrepid::FieldContainer<Real> > weighted_valUY_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*weighted_valUY_eval,
                                                               *weight_,
                                                               *valUY_eval);
    Intrepid::FunctionSpaceTools::integrate<Real>(*G[1],
                                                  *weighted_valUY_eval,
                                                  *(feVel_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Vertical_Velocity_Stokes::gradient_2 is zero.");
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    const int c  = feVel_->gradN()->dimension(0);
    const int fv = feVel_->gradN()->dimension(1);
    const int fp = fePrs_->gradN()->dimension(1);
    const int p = feVel_->cubPts()->dimension(1);
    // Initialize output grad
    Intrepid::FieldContainer<Real> velVX_grad(c, fv);
    Intrepid::FieldContainer<Real> velVY_grad(c, fv);
    Intrepid::FieldContainer<Real> presV_grad(c, fp);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > G;
    G.resize(fieldHelper_->numFields());
    G[0] = ROL::makePtrFromRef(velVX_grad);
    G[1] = ROL::makePtrFromRef(velVY_grad);
    G[2] = ROL::makePtrFromRef(presV_grad);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);

    ROL::Ptr<Intrepid::FieldContainer<Real> > valV_eval = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    feVel_->evaluateValue(valV_eval, V[1]);
    ROL::Ptr<Intrepid::FieldContainer<Real> > weighted_valV_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*weighted_valV_eval,
                                                               *weight_,
                                                               *valV_eval);
    Intrepid::FunctionSpaceTools::integrate<Real>(*G[1],
                                                  *weighted_valV_eval,
                                                  *(feVel_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    fieldHelper_->combineFieldCoeff(hess, G);   
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Vertical_Velocity_Stokes::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Vertical_Velocity_Stokes::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Vertical_Velocity_Stokes::HessVec_22 is zero.");
  }

}; // QoI_Vertical_Velocity_Stokes

template <class Real>
class QoI_L2Penalty_Stokes : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real> > feVel_;
  const ROL::Ptr<FE<Real> > fePrs_;
  const ROL::Ptr<FieldHelper<Real> > fieldHelper_;

public:
  QoI_L2Penalty_Stokes(const ROL::Ptr<FE<Real> > &feVel,
                             const ROL::Ptr<FE<Real> > &fePrs,
                             const ROL::Ptr<FieldHelper<Real> > &fieldHelper)
    : feVel_(feVel), fePrs_(fePrs), fieldHelper_(fieldHelper) {}

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real> > & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    const int c = feVel_->gradN()->dimension(0);
    int p = feVel_->cubPts()->dimension(1);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);

    // Build local state tracking term
    ROL::Ptr<Intrepid::FieldContainer<Real> > valX = ROL::makePtr<Intrepid::FieldContainer<Real> >(c);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valZX_eval = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    feVel_->evaluateValue(valZX_eval, Z[0]);
    feVel_->computeIntegral(valX,valZX_eval,valZX_eval);
    Intrepid::RealSpaceTools<Real>::scale(*valX,static_cast<Real>(0.5));

    // Build local state tracking term
    ROL::Ptr<Intrepid::FieldContainer<Real> > valY = ROL::makePtr<Intrepid::FieldContainer<Real> >(c);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valZY_eval = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    feVel_->evaluateValue(valZY_eval, Z[1]);
    feVel_->computeIntegral(valY,valZY_eval,valZY_eval);
    Intrepid::RealSpaceTools<Real>::scale(*valY,static_cast<Real>(0.5));

    // x-grad penalty
    ROL::Ptr<Intrepid::FieldContainer<Real> > valgradX = ROL::makePtr<Intrepid::FieldContainer<Real> >(c);
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradX_eval;
    gradX_eval = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p, 2);
    feVel_->evaluateGradient(gradX_eval, Z[0]);
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradXNorm_eval;
    gradXNorm_eval = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    for(int i = 0; i < c; i++)
      {
	for(int j = 0; j < p; j++)
	  {
	    (*gradXNorm_eval)(i,j) = (*gradX_eval)(i,j,0)*(*gradX_eval)(i,j,0) + (*gradX_eval)(i,j,1)*(*gradX_eval)(i,j,1);
	  }
      }
    // Create array of 1's to integrate against
    ROL::Ptr<Intrepid::FieldContainer<Real> > const_weights = ROL::makePtr<Intrepid::FieldContainer<Real> >(c,p);
    const_weights->initialize(1.0);
    feVel_->computeIntegral(valgradX,gradXNorm_eval,const_weights);
    Intrepid::RealSpaceTools<Real>::scale(*valgradX,static_cast<Real>(0.5));
    
    // y-grad penalty
    ROL::Ptr<Intrepid::FieldContainer<Real> > valgradY = ROL::makePtr<Intrepid::FieldContainer<Real> >(c);
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradY_eval;
    gradY_eval = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p, 2);
    feVel_->evaluateGradient(gradY_eval, Z[1]);
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradYNorm_eval;
    gradYNorm_eval = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    for(int i = 0; i < c; i++)
      {
	for(int j = 0; j < p; j++)
	  {
	    (*gradYNorm_eval)(i,j) = (*gradY_eval)(i,j,0)*(*gradY_eval)(i,j,0) + (*gradY_eval)(i,j,1)*(*gradY_eval)(i,j,1);
	  }
      }
    // Create array of 1's to integrate against
    feVel_->computeIntegral(valgradY,gradYNorm_eval,const_weights);
    Intrepid::RealSpaceTools<Real>::scale(*valgradY,static_cast<Real>(0.5));

    for(int i = 0; i < c; i++)
      {
	(*val)(i) = (*valX)(i) + (*valY)(i) + (*valgradX)(i) + (*valgradY)(i);
      }
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_Stokes::gradient_1 is zero.");
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {

    // Get relevant dimensions
    const int c  = feVel_->gradN()->dimension(0);
    int p = feVel_->cubPts()->dimension(1);
    const int fv = feVel_->gradN()->dimension(1);
    const int fp = fePrs_->gradN()->dimension(1);
    // Initialize output grad
    Intrepid::FieldContainer<Real> velZX_grad(c, fv);
    Intrepid::FieldContainer<Real> velZY_grad(c, fv);
    Intrepid::FieldContainer<Real> presZ_grad(c, fp);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > G;
    G.resize(fieldHelper_->numFields());
    G[0] = ROL::makePtrFromRef(velZX_grad);
    G[1] = ROL::makePtrFromRef(velZY_grad);
    G[2] = ROL::makePtrFromRef(presZ_grad);
    // Build local gradient of state tracking term
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valZX_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    feVel_->evaluateValue(valZX_eval, Z[0]);
    Intrepid::FunctionSpaceTools::integrate<Real>(*G[0],
                                                  *valZX_eval,
                                                  *(feVel_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valZY_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    feVel_->evaluateValue(valZY_eval, Z[1]);
    Intrepid::FunctionSpaceTools::integrate<Real>(*G[1],
                                                  *valZY_eval,
                                                  *(feVel_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);

    // Build local gradient of state tracking term
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradX_eval;
    gradX_eval = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p, 2);
    feVel_->evaluateGradient(gradX_eval, Z[0]);
    ROL::Ptr<Intrepid::FieldContainer<Real> > diffXx;
    diffXx = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > diffXy;
    diffXy = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    for(int i = 0; i < c; i++)
      {
	for(int j = 0; j < p; j++)
	  {
	    (*diffXx)(i,j) = (*gradX_eval)(i,j,0);
	    (*diffXy)(i,j) = (*gradX_eval)(i,j,1);
	  }
      }
    Intrepid::FunctionSpaceTools::integrate<Real>(*G[0],
                                                  *diffXx,
                                                  *(feVel_->DNDdetJ(0)),
                                                  Intrepid::COMP_CPP, true);
    Intrepid::FunctionSpaceTools::integrate<Real>(*G[0],
                                                  *diffXy,
                                                  *(feVel_->DNDdetJ(1)),
                                                  Intrepid::COMP_CPP, true);

    // Build local gradient of state tracking term
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradY_eval;
    gradY_eval = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p, 2);
    feVel_->evaluateGradient(gradY_eval, Z[1]);
    ROL::Ptr<Intrepid::FieldContainer<Real> > diffYx;
    diffYx = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > diffYy;
    diffYy = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    for(int i = 0; i < c; i++)
      {
	for(int j = 0; j < p; j++)
	  {
	    (*diffYx)(i,j) = (*gradY_eval)(i,j,0);
	    (*diffYy)(i,j) = (*gradY_eval)(i,j,1);
	  }
      }
    Intrepid::FunctionSpaceTools::integrate<Real>(*G[1],
                                                  *diffYx,
                                                  *(feVel_->DNDdetJ(0)),
                                                  Intrepid::COMP_CPP, true);
    Intrepid::FunctionSpaceTools::integrate<Real>(*G[1],
                                                  *diffYy,
                                                  *(feVel_->DNDdetJ(1)),
                                                  Intrepid::COMP_CPP, true);


    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_Stokes::HessVec_11 is zero.");
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_Stokes::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_Stokes::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    const int c  = feVel_->gradN()->dimension(0);
    const int fv = feVel_->gradN()->dimension(1);
    const int fp = fePrs_->gradN()->dimension(1);
    const int p = feVel_->cubPts()->dimension(1);
    // Initialize output grad
    Intrepid::FieldContainer<Real> velVX_grad(c, fv);
    Intrepid::FieldContainer<Real> velVY_grad(c, fv);
    Intrepid::FieldContainer<Real> presV_grad(c, fp);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > G;
    G.resize(fieldHelper_->numFields());
    G[0] = ROL::makePtrFromRef(velVX_grad);
    G[1] = ROL::makePtrFromRef(velVY_grad);
    G[2] = ROL::makePtrFromRef(presV_grad);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);

    ROL::Ptr<Intrepid::FieldContainer<Real> > valVX_eval = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    feVel_->evaluateValue(valVX_eval, V[0]);
    Intrepid::FunctionSpaceTools::integrate<Real>(*G[0],
                                                  *valVX_eval,
                                                  *(feVel_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valVY_eval = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    feVel_->evaluateValue(valVY_eval, V[1]);
    Intrepid::FunctionSpaceTools::integrate<Real>(*G[1],
                                                  *valVY_eval,
                                                  *(feVel_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);

    // Build local gradient of state tracking term
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradX_eval;
    gradX_eval = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p, 2);
    feVel_->evaluateGradient(gradX_eval, V[0]);
    ROL::Ptr<Intrepid::FieldContainer<Real> > diffXx;
    diffXx = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > diffXy;
    diffXy = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    for(int i = 0; i < c; i++)
      {
	for(int j = 0; j < p; j++)
	  {
	    (*diffXx)(i,j) = (*gradX_eval)(i,j,0);
	    (*diffXy)(i,j) = (*gradX_eval)(i,j,1);
	  }
      }
    Intrepid::FunctionSpaceTools::integrate<Real>(*G[0],
                                                  *diffXx,
                                                  *(feVel_->DNDdetJ(0)),
                                                  Intrepid::COMP_CPP, true);
    Intrepid::FunctionSpaceTools::integrate<Real>(*G[0],
                                                  *diffXy,
                                                  *(feVel_->DNDdetJ(1)),
                                                  Intrepid::COMP_CPP, true);

    // Build local gradient of state tracking term
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradY_eval;
    gradY_eval = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p, 2);
    feVel_->evaluateGradient(gradY_eval, V[1]);
    ROL::Ptr<Intrepid::FieldContainer<Real> > diffYx;
    diffYx = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > diffYy;
    diffYy = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    for(int i = 0; i < c; i++)
      {
	for(int j = 0; j < p; j++)
	  {
	    (*diffYx)(i,j) = (*gradY_eval)(i,j,0);
	    (*diffYy)(i,j) = (*gradY_eval)(i,j,1);
	  }
      }
    Intrepid::FunctionSpaceTools::integrate<Real>(*G[1],
                                                  *diffYx,
                                                  *(feVel_->DNDdetJ(0)),
                                                  Intrepid::COMP_CPP, true);
    Intrepid::FunctionSpaceTools::integrate<Real>(*G[1],
                                                  *diffYy,
                                                  *(feVel_->DNDdetJ(1)),
                                                  Intrepid::COMP_CPP, true);


    fieldHelper_->combineFieldCoeff(hess, G);
  }

}; // QoI_L2Penalty_Stokes

template <class Real>
class StdObjective_Stokes : public ROL::StdObjective<Real> {
private:
  Real alpha_;

public:
  StdObjective_Stokes(Teuchos::ParameterList &parlist) {
    alpha_    = parlist.sublist("Problem").get("Control penalty parameter",1.e-4);
  }

  Real value(const std::vector<Real> &x, Real &tol) {
    Real val = alpha_*x[1];
    val += x[0];
    return val;
  }

  void gradient(std::vector<Real> &g, const std::vector<Real> &x, Real &tol) {
    const Real one(1);
    g[0] = one;
    g[1] = alpha_;
  }

  void hessVec(std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol) {
    const Real zero(0);
    hv[0] = zero;
    hv[1] = zero;

  }

}; // OBJ_SCALAR

#endif
