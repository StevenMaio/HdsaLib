#ifndef PDE_DARCY_ADV_DIFF_HPP
#define PDE_DARCY_ADV_DIFF_HPP

#include "../../../PDE-OPT/TOOLS/pde.hpp"
#include "../../../PDE-OPT/TOOLS/fe.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "ROL_Ptr.hpp"

template <class Real>
class PDE_darcy_adv_diff : public PDE<Real> {
private:
  // Finite element basis information
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > basisPtrPrs_;
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > basisPtrCntm_;
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > basisPtrs_;
  // Cell cubature information
  ROL::Ptr<Intrepid::Cubature<Real> > cellCub_;
  ROL::Ptr<Intrepid::Cubature<Real> > bdryCub_;
  // Cell node information
  ROL::Ptr<Intrepid::FieldContainer<Real> > volCellNodes_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellNodes_;
  std::vector<std::vector<std::vector<int> > > bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<FE<Real> > fePrs_;
  ROL::Ptr<FE<Real> > feCntm_;
  std::vector<std::vector<ROL::Ptr<FE<Real> > > > fePrsBdry_;
  std::vector<std::vector<ROL::Ptr<FE<Real> > > > feCntmBdry_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int> > fpidx_;
  std::vector<std::vector<int> > fcidx_;

  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellPDofValues_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellCDofValues_;

  // Field pattern, offsets, etc.
  std::vector<std::vector<int> > fieldPattern_;  // local Field/DOF pattern; set from DOF manager 
  int numFields_;                                // number of fields (equations in the PDE)
  int numDofs_;                                  // total number of degrees of freedom for all (local) fields
  std::vector<int> offset_;                      // for each field, a counting offset
  std::vector<int> numFieldDofs_;                // for each field, number of degrees of freedom

  ROL::Ptr<FieldHelper<Real> > fieldHelper_;

  Real T_;

  Real a;
  int L;
  std::vector<Real> uncertain_basis_grid;

public:
  PDE_darcy_adv_diff(Teuchos::ParameterList &parlist) {
    // Finite element fields -- NOT DIMENSION INDEPENDENT!
    basisPtrPrs_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >>();
    basisPtrCntm_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >>();
    // Volume quadrature rules.
    shards::CellTopology cellType = basisPtrPrs_->getBaseCellTopology();         // get the cell type from any basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                           // create cubature factory
    int cubDegree = parlist.sublist("Problem").get("Cubature Degree", 4);        // set cubature degree, e.g., 4
    cellCub_ = cubFactory.create(cellType, cubDegree);                           // create default cubature
    // Boundary quadrature rules.
    int d = cellType.getDimension();
    shards::CellTopology bdryCellType = cellType.getCellTopologyData(d-1, 0);
    int bdryCubDegree = parlist.sublist("Problem").get("Boundary Cubature Degree",4); // set cubature degree, e.g., 4
    bdryCub_ = cubFactory.create(bdryCellType, bdryCubDegree);
    // Fill finite element basis container
    basisPtrs_.clear();
    basisPtrs_.push_back(basisPtrPrs_); // Pressure
    basisPtrs_.push_back(basisPtrCntm_); // Contaminant

    numDofs_ = 0;
    numFields_ = basisPtrs_.size();
    offset_.resize(numFields_);
    numFieldDofs_.resize(numFields_);
    for (int i=0; i<numFields_; ++i) {
      if (i==0) {
        offset_[i] = 0;
      }
      else {
        offset_[i] = offset_[i-1] + basisPtrs_[i-1]->getCardinality();
      }
      numFieldDofs_[i] = basisPtrs_[i]->getCardinality();
      numDofs_ += numFieldDofs_[i];
    }

    L = parlist.sublist("Problem").get("Number of Uncertainty Basis Function", 10);
    a = parlist.sublist("Problem").get("Noise Level", .2);
    uncertain_basis_grid.resize(L+1);
    for(int i = 0; i < L+1; i++)
      {
	uncertain_basis_grid[i] = static_cast<Real>(i)/static_cast<Real>(L);
      }

    T_  = parlist.sublist("Time Discretization").get("End Time",1.0);
  }

   void residual(ROL::Ptr<Intrepid::FieldContainer<Real> > & res,
                const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    
     // Retrieve dimensions.
    int c  = u_coeff->dimension(0);
    int p  = cellCub_->getNumPoints();
    int fp = basisPtrPrs_->getCardinality();
    int fc = basisPtrCntm_->getCardinality();
    int d = fePrs_->gradN()->dimension(3);
 
    // Initialize residuals.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > R(2);
    R[0]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    R[1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fc);

    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);

    // COMPUTE PDE COEFFICIENTS
    ROL::Ptr<Intrepid::FieldContainer<Real> > kappa
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > epsilon
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > rhs
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > Prs_source
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    Intrepid::RealSpaceTools<Real>::scale(*Prs_source, 0.0);
    computeCoefficients(epsilon,rhs);

    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    fePrs_->evaluateValue(kappa, Z[0]);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*kappa)(i,j) = std::exp((*kappa)(i,j));
      }
    }

    // Evaluate/interpolate gradient of finite element fields on cells.
    // Comptue grad(Prs)
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradPrs_eval;
    gradPrs_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    fePrs_->evaluateGradient(gradPrs_eval, U[0]);
    // Compute grad(Cntm)
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradCntm_eval;
    gradCntm_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    feCntm_->evaluateGradient(gradCntm_eval, U[1]);
    // Compute Cntm
    ROL::Ptr<Intrepid::FieldContainer<Real> > Cntm_eval;
    Cntm_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    feCntm_->evaluateValue(Cntm_eval, U[1]);


    // Pressure Equation

    // COMPUTE PERMEABILITY TERM
    // Multiply kappa * grad(Prs)
    Intrepid::FieldContainer<Real> kappa_gradPrs(c, p, d);
    Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(kappa_gradPrs,
                                                               *kappa,
                                                               *gradPrs_eval);
    // Integrate (kappa * grad(Prs)) . grad(N)
    Intrepid::FunctionSpaceTools::integrate<Real>(*R[0],
                                                  kappa_gradPrs,
                                                  *(fePrs_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, false);

    Intrepid::FunctionSpaceTools::integrate<Real>(*R[0],
						  *Prs_source,
						  *(fePrs_->NdetJ()),
						  Intrepid::COMP_CPP, true);

    // Contaminant Equation

    // COMPUTE DIFFUSION TERM
    // Multiply epsilon * grad(Cntm)
    Intrepid::FieldContainer<Real> epsilon_gradCntm(c, p, d);
    Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(epsilon_gradCntm,
                                                               *epsilon,
                                                               *gradCntm_eval);
    // Integrate (epsilon * grad(Ctnm)) . grad(N)
    Intrepid::FunctionSpaceTools::integrate<Real>(*R[1],
                                                  epsilon_gradCntm,
                                                  *(feCntm_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, false);

    // Compute advection terms
    Intrepid::FieldContainer<Real> kappa_gradPrs_gradCntm_eval(c, p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(kappa_gradPrs_gradCntm_eval,
    							    kappa_gradPrs,
    							    *gradCntm_eval);
    Intrepid::RealSpaceTools<Real>::scale(kappa_gradPrs_gradCntm_eval, -1.0);
    Intrepid::FunctionSpaceTools::integrate<Real>(*R[1],
                                                  kappa_gradPrs_gradCntm_eval,
                                                  (*feCntm_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);

    Intrepid::RealSpaceTools<Real>::scale(*Prs_source, -1.0);
    Intrepid::FieldContainer<Real> Prs_source_Cntm_eval(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(Prs_source_Cntm_eval,*Prs_source,*Cntm_eval);
    
    Intrepid::FunctionSpaceTools::integrate<Real>(*R[1],
                                                  Prs_source_Cntm_eval,
                                                  (*feCntm_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);

    // ADD RHS TO RESIDUAL
    Intrepid::FunctionSpaceTools::integrate<Real>(*R[1],
                                                  *rhs,
                                                  (*feCntm_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);
    
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
    int fp = basisPtrPrs_->getCardinality();
    int fc = basisPtrCntm_->getCardinality();
    int d = fePrs_->gradN()->dimension(3);

    // Initialize jacobians.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > J(2);
    for (int i = 0; i < 2; ++i) {
      J[i].resize(2);
    }
    J[0][0] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fp);
    J[1][0] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fc, fp);
    J[0][1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fc);
    J[1][1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fc, fc);

    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    
    // COMPUTE PDE COEFFICIENTS
    ROL::Ptr<Intrepid::FieldContainer<Real> > kappa
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > epsilon
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > rhs
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > Prs_source
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    Intrepid::RealSpaceTools<Real>::scale(*Prs_source, 0.0);
    computeCoefficients(epsilon,rhs);

    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    fePrs_->evaluateValue(kappa, Z[0]);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*kappa)(i,j) = std::exp((*kappa)(i,j));
      }
    }

    // Evaluate/interpolate gradient of finite element fields on cells.
    // Comptue grad(Prs)
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradPrs_eval;
    gradPrs_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    fePrs_->evaluateGradient(gradPrs_eval, U[0]);
    // Compute grad(Cntm)
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradCntm_eval;
    gradCntm_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    feCntm_->evaluateGradient(gradCntm_eval, U[1]);


    // Pressure Equation

    // COMPUTE DIFFUSION TERM
    // Multiply kappa * grad(N)
    Intrepid::FieldContainer<Real> kappa_gradN(c, fp, p, d);
    Intrepid::FunctionSpaceTools::tensorMultiplyDataField<Real>(kappa_gradN,
                                                                *kappa,
                                                                *(fePrs_->gradN()));
    // Integrate (kappa * grad(N)) . grad(N)
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][0],
                                                  kappa_gradN,
                                                  *(fePrs_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, false);
    // Contaminant Equation

    // COMPUTE DIFFUSION TERM
    // Multiply epsilon * grad(N)
    Intrepid::FieldContainer<Real> epsilon_gradN(c, fc, p, d);
    Intrepid::FunctionSpaceTools::tensorMultiplyDataField<Real>(epsilon_gradN,
                                                                *epsilon,
                                                                *(feCntm_->gradN()));
    // Integrate (epsilon * grad(N)) . grad(N)
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[1][1],
                                                  epsilon_gradN,
                                                  *(feCntm_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, false);

    // Comptue source terms
    Intrepid::RealSpaceTools<Real>::scale(*Prs_source, -1.0);

    Intrepid::FieldContainer<Real> Prs_source_N(c, fp, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(Prs_source_N,
                                                                *Prs_source,
                                                                *(feCntm_->N()));
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[1][1],
                                                  Prs_source_N,
                                                  *(feCntm_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);

    ROL::Ptr<Intrepid::FieldContainer<Real> > kappa_gradPrs = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p, d);
    Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(*kappa_gradPrs,
							       *kappa,
							       *gradPrs_eval);
    Intrepid::RealSpaceTools<Real>::scale(*kappa_gradPrs,static_cast<Real>(-1.0));
    Intrepid::FieldContainer<Real> kappa_gradPrs_gradN(c, fc, p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(kappa_gradPrs_gradN,
                                                             *kappa_gradPrs,
                                                             *(feCntm_->gradN()));
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[1][1],
                                                  *(feCntm_->NdetJ()),
                                                  kappa_gradPrs_gradN,
                                                  Intrepid::COMP_CPP, true);

    // Integrate (kappa * grad(N) . grad(Cntm)) * N
    Intrepid::FieldContainer<Real> kappa_gradN_gradCntm(c, fc, p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(kappa_gradN_gradCntm,
                                                             *gradCntm_eval,
							     kappa_gradN);
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[1][0],
                                                  *(feCntm_->NdetJ()),
                                                  kappa_gradN_gradCntm,
                                                  Intrepid::COMP_CPP, false);

    Intrepid::RealSpaceTools<Real>::scale(*J[1][0],static_cast<Real>(-1.0));

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
    int fp = basisPtrPrs_->getCardinality();
    int fc = basisPtrCntm_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize jacobians.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > J(2);
    for (int i = 0; i < 2; ++i) {
      J[i].resize(2);
    }
    J[0][0]     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fp);
    J[1][0]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fc, fp);
    J[0][1]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fc);
    J[1][1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fc, fc);

    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U, Z;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    fieldHelper_->splitFieldCoeff(Z, z_coeff);

    // Comptue grad(Prs)
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradPrs_eval;
    gradPrs_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    fePrs_->evaluateGradient(gradPrs_eval, U[0]);
    // Compute grad(Cntm)
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradCntm_eval;
    gradCntm_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    feCntm_->evaluateGradient(gradCntm_eval, U[1]);

    ROL::Ptr<Intrepid::FieldContainer<Real> > kappa =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    fePrs_->evaluateValue(kappa, Z[0]);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*kappa)(i,j) = std::exp((*kappa)(i,j));
      }
    }

    Intrepid::FieldContainer<Real> gradPrs_gradN_eval(c, fp, p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(gradPrs_gradN_eval,*gradPrs_eval,*(feCntm_->gradN()));
    Intrepid::FieldContainer<Real> gradPrs_gradN_kappa_eval(c, fp, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(gradPrs_gradN_kappa_eval,
								*kappa,
								gradPrs_gradN_eval);
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][0],
                                                  gradPrs_gradN_kappa_eval,
                                                  *(feCntm_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);



    // Cross term
    Intrepid::FieldContainer<Real> gradPrs_gradCntm_eval(c, p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(gradPrs_gradCntm_eval,*gradPrs_eval,*gradCntm_eval);
    Intrepid::RealSpaceTools<Real>::scale(gradPrs_gradCntm_eval,-1.0);
    Intrepid::FieldContainer<Real> gradPrs_gradCntm_N_eval(c, fp, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(gradPrs_gradCntm_N_eval,
								gradPrs_gradCntm_eval,
								*(feCntm_->N()));
    Intrepid::FieldContainer<Real> gradPrs_gradCntm_N_kappa_eval(c, fp, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(gradPrs_gradCntm_N_kappa_eval,
								*kappa,
							        gradPrs_gradCntm_N_eval);
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[1][0],
						  gradPrs_gradCntm_N_kappa_eval,
						  *(feCntm_->NdetJ()),
						  Intrepid::COMP_CPP, false);
    
    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);
  }

  void Jacobian_3(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_darcy_adv_diff::Jacobian_3): Jacobian is zero.");
  }

  void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> (PDE_darcy_adv_diff::Hessian_11): Not implemented.");
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> (PDE_darcy_adv_diff::Hessian_12): Not implemented.");
  }

  void Hessian_13(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_darcy_adv_diff::Hessian_13): Hessian is zero.");
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> (PDE_darcy_adv_diff::Hessian_21): Not implemented.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> (PDE_darcy_adv_diff::Hessian_22): Not implemented.");
  }

  void Hessian_23(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_darcy_adv_diff::Hessian_23): Hessian is zero.");
  }

  void Hessian_31(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_darcy_adv_diff::Hessian_31): Hessian is zero.");
  }

  void Hessian_32(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_darcy_adv_diff::Hessian_32): Hessian is zero.");
  }

  void Hessian_33(std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_darcy_adv_diff::Hessian_33): Hessian is zero.");
  }

   void RieszMap_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & riesz) {
     throw Exception::NotImplemented(">>> (PDE_darcy_adv_diff::RieszMap_1): Not implemented.");
  }

  void RieszMap_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & riesz) {
    throw Exception::NotImplemented(">>> (PDE_darcy_adv_diff::RieszMap_2): Not implemented.");
  }

  void setCellNodes(const ROL::Ptr<Intrepid::FieldContainer<Real> > &volCellNodes,
                    const std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int> > > &bdryCellLocIds) {
    volCellNodes_ = volCellNodes;
    bdryCellNodes_ = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    // Finite element definition.
    fePrs_ = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtrPrs_,cellCub_);
    feCntm_ = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtrCntm_,cellCub_);
    // Get boundary degrees of freedom.
    fpidx_ = fePrs_->getBoundaryDofs();
    fcidx_ = feCntm_->getBoundaryDofs();
    // Construct boundary FEs
    const int numSideSets = bdryCellNodes.size();
    fePrsBdry_.resize(numSideSets);
    feCntmBdry_.resize(numSideSets);
    for (int i = 0; i < numSideSets; ++i) {
      int numLocSides = bdryCellNodes[i].size();
      fePrsBdry_[i].resize(numLocSides);
      feCntmBdry_[i].resize(numLocSides);
      for (int j = 0; j < numLocSides; ++j) {
        if (bdryCellNodes[i][j] != ROL::nullPtr) {
          fePrsBdry_[i][j] = ROL::makePtr<FE<Real>>(bdryCellNodes[i][j],basisPtrPrs_,bdryCub_,j);
          feCntmBdry_[i][j] = ROL::makePtr<FE<Real>>(bdryCellNodes[i][j],basisPtrCntm_,bdryCub_,j);
        }
      }
    }

  }

  void setFieldPattern(const std::vector<std::vector<int> > & fieldPattern) {
    fieldPattern_ = fieldPattern;
    fieldHelper_ = ROL::makePtr<FieldHelper<Real>>(numFields_, numDofs_, numFieldDofs_, fieldPattern_);
  }


  const ROL::Ptr<FE<Real> > getPressureFE(void) const {
    return fePrs_;
  }

  const ROL::Ptr<FE<Real> > getContaminantFE(void) const {
    return feCntm_;
  }

  const std::vector<std::vector<ROL::Ptr<FE<Real> > > > getPressureBdryFE(void) const {
    return fePrsBdry_;
  }

  const std::vector<std::vector<ROL::Ptr<FE<Real> > > > getContaminantBdryFE(void) const {
    return feCntmBdry_;
  }

  const std::vector<std::vector<std::vector<int> > > getBdryCellLocIds(void) const {
    return bdryCellLocIds_;
  }

  const ROL::Ptr<FieldHelper<Real> > getFieldHelper(void) const {
    return fieldHelper_;
  }


  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > getFields() {
    return basisPtrs_;
  }

private:
  void computeCoefficients(ROL::Ptr<Intrepid::FieldContainer<Real> > &epsilon,
                           ROL::Ptr<Intrepid::FieldContainer<Real> > &rhs) const {
    int c = fePrs_->gradN()->dimension(0);
    int p = fePrs_->gradN()->dimension(2);
    int d = fePrs_->gradN()->dimension(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k) {
          pt[k] = (*fePrs_->cubPts())(i,j,k);
        }
        // Compute diffusion epsilon
        (*epsilon)(i,j) = evaluateDiffusion(pt);
        // Compute forcing term f
        (*rhs)(i,j) = -evaluateRHS(pt);
      }
    }    

  }
  
  Real evaluateDiffusion(const std::vector<Real> & x) const {
    Real val = 0.025;

    const std::vector<Real> param = PDE<Real>::getParameter();
    if((int)param.size()>0)
      {
    	Real val_uncertain = 0.0;
    	int offset = 16*9;
    	int count = 0;
        val_uncertain = param[offset];
    	val_uncertain = a*val_uncertain;
    	val_uncertain = 1+val_uncertain;
    	val = val*val_uncertain;
      }


    return val;
  }

  Real evaluateRHS(const std::vector<Real> & x) const {
    int N = 16;
    const std::vector<Real> xl = {0.2, 0.2, 0.2, 0.2, 0.4, 0.4, 0.4, 0.4, 0.6, 0.6, 0.6, 0.6, 0.8, 0.8, 0.8, 0.8};
    const std::vector<Real> yl = {0.2, 0.4, 0.6, 0.8, 0.2, 0.4, 0.6, 0.8, 0.2, 0.4, 0.6, 0.8, 0.2, 0.4, 0.6, 0.8};
    Real val = 0.0;

    for(int k = 0; k < N; k++)
      {
	val += 10.0*std::exp(-100.0 * ( (x[0]-xl[k])*(x[0]-xl[k])+(x[1]-yl[k])*(x[1]-yl[k]) ) );
      }

    const std::vector<Real> param = PDE<Real>::getParameter();
    if((int)param.size()>0)
      {
	Real val_uncertain = 0.0;
	Real basis_fun = 0.0;
	for(int k = 0; k < N; k++)
	  {
	    if(std::abs(x[0]-xl[k])<.1 && std::abs(x[1]-yl[k])<.1)
	      {
		int count = 0;
		for(int j = 0; j < 3; j++)
		  {
		    for(int i = 0; i < 3; i++)
		      {
			basis_fun = Parameter_Basis_Fun_Box_Eval(x[0]-xl[k],i);
			basis_fun *= Parameter_Basis_Fun_Box_Eval(x[1]-yl[k],j);
			val_uncertain += param[9*k+count]*basis_fun;
			count += 1;
		      }
		  }
	      }
	  }
    	val_uncertain = a*val_uncertain;
    	val_uncertain = 1+val_uncertain;
    	val = val*val_uncertain;
      }

    return val;
  }

  Real Parameter_Basis_Fun_Box_Eval(const Real & x, int i) const
  {
    // Evaluates ith 1D basis function on a box with grid of 3x3 nodes at point x
    Real val = 0.0;
    Real dist = std::abs(x-static_cast<Real>(i-1)*.2/4.0);
    if( dist < .2/4.0 )
      {
	val = 1.0 - dist*4.0/.2;
      }
    return val;
  }

}; // PDE_darcy_adv_diff

#endif
