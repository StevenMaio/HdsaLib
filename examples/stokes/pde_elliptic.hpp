#ifndef PDE_ELLIPTIC_HPP
#define PDE_ELLIPTIC_HPP

#include "../../../PDE-OPT/TOOLS/pde.hpp"
#include "../../../PDE-OPT/TOOLS/fe.hpp"
#include "../../../PDE-OPT/TOOLS/fieldhelper.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "ROL_Ptr.hpp"

template <class Real>
class PDE_Elliptic : public PDE<Real> {
private:
  // Finite element basis information
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > basisPtrVel_;
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > basisPtrPrs_;
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > basisPtrs_;
  // Cell cubature information
  ROL::Ptr<Intrepid::Cubature<Real> > cellCub_;
  ROL::Ptr<Intrepid::Cubature<Real> > bdryCub_;
  // Cell node information
  ROL::Ptr<Intrepid::FieldContainer<Real> > volCellNodes_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellNodes_;
  std::vector<std::vector<std::vector<int> > > bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<FE<Real> > feVel_;
  ROL::Ptr<FE<Real> > fePrs_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int> > fvidx_;
  std::vector<std::vector<int> > fpidx_;
  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellDofValues_;
  // Field pattern, offsets, etc.
  std::vector<std::vector<int> > fieldPattern_;  // local Field/DOF pattern; set from DOF manager 
  int numFields_;                                // number of fields (equations in the PDE)
  int numDofs_;                                  // total number of degrees of freedom for all (local) fields
  std::vector<int> offset_;                      // for each field, a counting offset
  std::vector<int> numFieldDofs_;                // for each field, number of degrees of freedom
  
  // Problem parameters.
  Real diff_;
  Real eye_;

  ROL::Ptr<FieldHelper<Real> > fieldHelper_;

public:
  PDE_Elliptic(Teuchos::ParameterList &parlist, Real diff, Real eye) {
    // Finite element fields.
    basisPtrVel_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real> >>();
    basisPtrPrs_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >>();
    basisPtrs_.clear();
    basisPtrs_.push_back(basisPtrVel_);  // Velocity X
    basisPtrs_.push_back(basisPtrVel_);  // Velocity Y
    basisPtrs_.push_back(basisPtrPrs_);  // Pressure
    // Quadrature rules.
    shards::CellTopology cellType = basisPtrs_[0]->getBaseCellTopology();        // get the cell type from any basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                           // create cubature factory
    int cubDegree = parlist.sublist("Problem").get("Cubature Degree", 4);        // set cubature degree, e.g., 4
    cellCub_ = cubFactory.create(cellType, cubDegree);                           // create default cubature

    int d = cellType.getDimension();
    shards::CellTopology bdryCellType = cellType.getCellTopologyData(d-1, 0);
    int bdryCubDegree = parlist.sublist("Problem").get("Boundary Cubature Degree",4); // set cubature degree, e.g., 4
    bdryCub_ = cubFactory.create(bdryCellType, bdryCubDegree);

    // Other problem parameters.
    diff_ = diff;
    eye_ = eye;

    numDofs_ = 0;
    numFields_ = basisPtrs_.size();
    offset_.resize(numFields_);
    numFieldDofs_.resize(numFields_);
    for (int i=0; i<numFields_; ++i) {
      if (i==0) {
        offset_[i]  = 0;
      }
      else {
        offset_[i]  = offset_[i-1] + basisPtrs_[i-1]->getCardinality();
      }
      numFieldDofs_[i] = basisPtrs_[i]->getCardinality();
      numDofs_ += numFieldDofs_[i];
    }
  }

  void residual(ROL::Ptr<Intrepid::FieldContainer<Real> > & res,
                const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    int c  = u_coeff->dimension(0);
    int p  = cellCub_->getNumPoints();
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize residuals.
    Intrepid::FieldContainer<Real> velX_res(c, fv);
    Intrepid::FieldContainer<Real> velY_res(c, fv);
    Intrepid::FieldContainer<Real> pres_res(c, fp);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > R;
    R.resize(numFields_);
    R[0] = ROL::makePtrFromRef(velX_res);
    R[1] = ROL::makePtrFromRef(velY_res);
    R[2] = ROL::makePtrFromRef(pres_res);

    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);

    // First component
    ROL::Ptr<Intrepid::FieldContainer<Real> > U0_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    feVel_->evaluateValue(U0_eval, U[0]);
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradU0_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p, d);
    feVel_->evaluateGradient(gradU0_eval, U[0]);
    // Integrate U * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*R[0],
						  *U0_eval,
						  *(feVel_->NdetJ()),
						  Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*R[0],eye_);
    ROL::Ptr<Intrepid::FieldContainer<Real> > laplace_term0 = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, fv);
    Intrepid::FunctionSpaceTools::integrate<Real>(*laplace_term0, *gradU0_eval, *(feVel_->gradNdetJ()), Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*laplace_term0,diff_);
    Intrepid::RealSpaceTools<Real>::add(*R[0],*laplace_term0);

    // Second component
    ROL::Ptr<Intrepid::FieldContainer<Real> > U1_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    feVel_->evaluateValue(U1_eval, U[1]);
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradU1_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p, d);
    feVel_->evaluateGradient(gradU1_eval, U[1]);
    // Integrate U * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*R[1],
						  *U1_eval,
						  *(feVel_->NdetJ()),
						  Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*R[1],eye_);
    ROL::Ptr<Intrepid::FieldContainer<Real> > laplace_term1 = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, fv);
    Intrepid::FunctionSpaceTools::integrate<Real>(*laplace_term1, *gradU1_eval, *(feVel_->gradNdetJ()), Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*laplace_term1,diff_);
    Intrepid::RealSpaceTools<Real>::add(*R[1],*laplace_term1);

    // Third component
    ROL::Ptr<Intrepid::FieldContainer<Real> > U2_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    fePrs_->evaluateValue(U2_eval, U[2]);
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradU2_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p, d);
    fePrs_->evaluateGradient(gradU2_eval, U[2]);
    // Integrate U * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*R[2],
						  *U2_eval,
						  *(fePrs_->NdetJ()),
						  Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*R[2],eye_);
    ROL::Ptr<Intrepid::FieldContainer<Real> > laplace_term2 = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, fp);
    Intrepid::FunctionSpaceTools::integrate<Real>(*laplace_term2, *gradU2_eval, *(fePrs_->gradNdetJ()), Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*laplace_term2,diff_);
    Intrepid::RealSpaceTools<Real>::add(*R[2],*laplace_term2);

    // Combine the residuals.
    fieldHelper_->combineFieldCoeff(res, R);
  }

  void Jacobian_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    int c  = u_coeff->dimension(0);
    int p  = cellCub_->getNumPoints();
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize residuals.
    Intrepid::FieldContainer<Real> velXvelX_jac(c, fv, fv);
    Intrepid::FieldContainer<Real> velXvelY_jac(c, fv, fv);
    Intrepid::FieldContainer<Real> velXpres_jac(c, fv, fp);
    Intrepid::FieldContainer<Real> velYvelX_jac(c, fv, fv);
    Intrepid::FieldContainer<Real> velYvelY_jac(c, fv, fv);
    Intrepid::FieldContainer<Real> velYpres_jac(c, fv, fp);
    Intrepid::FieldContainer<Real> presvelX_jac(c, fp, fv);
    Intrepid::FieldContainer<Real> presvelY_jac(c, fp, fv);
    Intrepid::FieldContainer<Real> prespres_jac(c, fp, fp);
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > J;
    J.resize(numFields_);
    J[0].resize(numFields_);
    J[1].resize(numFields_);
    J[2].resize(numFields_);
    J[0][0] = ROL::makePtrFromRef(velXvelX_jac); J[0][1] = ROL::makePtrFromRef(velXvelY_jac); J[0][2] = ROL::makePtrFromRef(velXpres_jac);  
    J[1][0] = ROL::makePtrFromRef(velYvelX_jac); J[1][1] = ROL::makePtrFromRef(velYvelY_jac); J[1][2] = ROL::makePtrFromRef(velYpres_jac);  
    J[2][0] = ROL::makePtrFromRef(presvelX_jac); J[2][1] = ROL::makePtrFromRef(presvelY_jac); J[2][2] = ROL::makePtrFromRef(prespres_jac);  

    // First component
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][0],
                                                  *(feVel_->N()),
                                                  *(feVel_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*J[0][0],eye_);
    ROL::Ptr<Intrepid::FieldContainer<Real> > laplace_term0 = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, fv, fv);
    Intrepid::FunctionSpaceTools::integrate<Real>(*laplace_term0, *(feVel_->gradN()), *(feVel_->gradNdetJ()), Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*laplace_term0,diff_);
    Intrepid::RealSpaceTools<Real>::add(*J[0][0],*laplace_term0);

    // Second component
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[1][1],
                                                  *(feVel_->N()),
                                                  *(feVel_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*J[1][1],eye_);
    ROL::Ptr<Intrepid::FieldContainer<Real> > laplace_term1 = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, fv, fv);
    Intrepid::FunctionSpaceTools::integrate<Real>(*laplace_term1, *(feVel_->gradN()), *(feVel_->gradNdetJ()), Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*laplace_term1,diff_);
    Intrepid::RealSpaceTools<Real>::add(*J[1][1],*laplace_term1);

   // Second component
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[2][2],
                                                  *(fePrs_->N()),
                                                  *(fePrs_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*J[2][2],eye_);
    ROL::Ptr<Intrepid::FieldContainer<Real> > laplace_term2 = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, fp, fp);
    Intrepid::FunctionSpaceTools::integrate<Real>(*laplace_term2, *(fePrs_->gradN()), *(fePrs_->gradNdetJ()), Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*laplace_term2,diff_);
    Intrepid::RealSpaceTools<Real>::add(*J[2][2],*laplace_term2);
 
    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);
  }


  void Jacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    int c  = u_coeff->dimension(0);
    int p  = cellCub_->getNumPoints();
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int d  = cellCub_->getDimension();

    // Initialize residuals.
    Intrepid::FieldContainer<Real> velXvelX_jac(c, fv, fv);
    Intrepid::FieldContainer<Real> velXvelY_jac(c, fv, fv);
    Intrepid::FieldContainer<Real> velXpres_jac(c, fv, fp);
    Intrepid::FieldContainer<Real> velYvelX_jac(c, fv, fv);
    Intrepid::FieldContainer<Real> velYvelY_jac(c, fv, fv);
    Intrepid::FieldContainer<Real> velYpres_jac(c, fv, fp);
    Intrepid::FieldContainer<Real> presvelX_jac(c, fp, fv);
    Intrepid::FieldContainer<Real> presvelY_jac(c, fp, fv);
    Intrepid::FieldContainer<Real> prespres_jac(c, fp, fp);
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > J;
    J.resize(numFields_);
    J[0].resize(numFields_);
    J[1].resize(numFields_);
    J[2].resize(numFields_);
    J[0][0] = ROL::makePtrFromRef(velXvelX_jac); J[0][1] = ROL::makePtrFromRef(velXvelY_jac); J[0][2] = ROL::makePtrFromRef(velXpres_jac);
    J[1][0] = ROL::makePtrFromRef(velYvelX_jac); J[1][1] = ROL::makePtrFromRef(velYvelY_jac); J[1][2] = ROL::makePtrFromRef(velYpres_jac);
    J[2][0] = ROL::makePtrFromRef(presvelX_jac); J[2][1] = ROL::makePtrFromRef(presvelY_jac); J[2][2] = ROL::makePtrFromRef(prespres_jac);

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);
  }

  void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Elliptic::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Elliptic::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Elliptic::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Elliptic::Hessian_22): Hessian is zero.");
  }

  void RieszMap_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & riesz) {
    // Retrieve dimensions.
    int c  = feVel_->N()->dimension(0);
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
 
    // Initialize residuals.
    Intrepid::FieldContainer<Real> velXvelX_jac(c, fv, fv);
    Intrepid::FieldContainer<Real> velXvelY_jac(c, fv, fv);
    Intrepid::FieldContainer<Real> velXpres_jac(c, fv, fp);
    Intrepid::FieldContainer<Real> velYvelX_jac(c, fv, fv);
    Intrepid::FieldContainer<Real> velYvelY_jac(c, fv, fv);
    Intrepid::FieldContainer<Real> velYpres_jac(c, fv, fp);
    Intrepid::FieldContainer<Real> presvelX_jac(c, fp, fv);
    Intrepid::FieldContainer<Real> presvelY_jac(c, fp, fv);
    Intrepid::FieldContainer<Real> prespres_jac(c, fp, fp);
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > J;
    J.resize(numFields_);
    J[0].resize(numFields_);
    J[1].resize(numFields_);
    J[2].resize(numFields_);
    J[0][0] = ROL::makePtrFromRef(velXvelX_jac); J[0][1] = ROL::makePtrFromRef(velXvelY_jac); J[0][2] = ROL::makePtrFromRef(velXpres_jac);  
    J[1][0] = ROL::makePtrFromRef(velYvelX_jac); J[1][1] = ROL::makePtrFromRef(velYvelY_jac); J[1][2] = ROL::makePtrFromRef(velYpres_jac);  
    J[2][0] = ROL::makePtrFromRef(presvelX_jac); J[2][1] = ROL::makePtrFromRef(presvelY_jac); J[2][2] = ROL::makePtrFromRef(prespres_jac);  

    *(J[0][0]) = *(feVel_->stiffMat());
    Intrepid::RealSpaceTools<Real>::add(*(J[0][0]),*(feVel_->massMat()));
    *(J[1][1]) = *(feVel_->stiffMat());
    Intrepid::RealSpaceTools<Real>::add(*(J[1][1]),*(feVel_->massMat()));
    *(J[2][2]) = *(fePrs_->massMat());

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(riesz, J);
  }

  void RieszMap_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & riesz) {
    // Retrieve dimensions.
    int c  = feVel_->N()->dimension(0);
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
 
    // Initialize residuals.
    Intrepid::FieldContainer<Real> velXvelX_jac(c, fv, fv);
    Intrepid::FieldContainer<Real> velXvelY_jac(c, fv, fv);
    Intrepid::FieldContainer<Real> velXpres_jac(c, fv, fp);
    Intrepid::FieldContainer<Real> velYvelX_jac(c, fv, fv);
    Intrepid::FieldContainer<Real> velYvelY_jac(c, fv, fv);
    Intrepid::FieldContainer<Real> velYpres_jac(c, fv, fp);
    Intrepid::FieldContainer<Real> presvelX_jac(c, fp, fv);
    Intrepid::FieldContainer<Real> presvelY_jac(c, fp, fv);
    Intrepid::FieldContainer<Real> prespres_jac(c, fp, fp);
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > J;
    J.resize(numFields_);
    J[0].resize(numFields_);
    J[1].resize(numFields_);
    J[2].resize(numFields_);
    J[0][0] = ROL::makePtrFromRef(velXvelX_jac); J[0][1] = ROL::makePtrFromRef(velXvelY_jac); J[0][2] = ROL::makePtrFromRef(velXpres_jac);  
    J[1][0] = ROL::makePtrFromRef(velYvelX_jac); J[1][1] = ROL::makePtrFromRef(velYvelY_jac); J[1][2] = ROL::makePtrFromRef(velYpres_jac);  
    J[2][0] = ROL::makePtrFromRef(presvelX_jac); J[2][1] = ROL::makePtrFromRef(presvelY_jac); J[2][2] = ROL::makePtrFromRef(prespres_jac);  

    *(J[0][0]) = *(feVel_->stiffMat());
    Intrepid::RealSpaceTools<Real>::add(*(J[0][0]),*(feVel_->massMat()));
    *(J[1][1]) = *(feVel_->stiffMat());
    Intrepid::RealSpaceTools<Real>::add(*(J[1][1]),*(feVel_->massMat()));
    *(J[2][2]) = *(fePrs_->massMat());

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(riesz, J);
  }

  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > getFields() {
    return basisPtrs_;
  }

  void setCellNodes(const ROL::Ptr<Intrepid::FieldContainer<Real> > &volCellNodes,
                    const std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int> > > &bdryCellLocIds) {
    volCellNodes_ = volCellNodes;
    bdryCellNodes_ = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    // Finite element definition.
    feVel_ = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtrVel_,cellCub_);
    fePrs_ = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtrPrs_,cellCub_);
    fvidx_ = feVel_->getBoundaryDofs();
    fpidx_ = fePrs_->getBoundaryDofs();
  }

  void setFieldPattern(const std::vector<std::vector<int> > & fieldPattern) {
    fieldPattern_ = fieldPattern;
    fieldHelper_ = ROL::makePtr<FieldHelper<Real>>(numFields_, numDofs_, numFieldDofs_, fieldPattern_);
  }

  const ROL::Ptr<FE<Real> > getVelocityFE(void) const {
    return feVel_;
  }

  const ROL::Ptr<FE<Real> > getPressureFE(void) const {
    return fePrs_;
  }

  const ROL::Ptr<FieldHelper<Real> > getFieldHelper(void) const {
    return fieldHelper_;
  }

}; // PDE_Elliptic

#endif
