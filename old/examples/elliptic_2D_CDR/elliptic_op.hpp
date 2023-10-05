#ifndef ELLIPTIC_OP_HPP
#define ELLIPTIC_OP_HPP

#include "../../../PDE-OPT/TOOLS/pde.hpp"
#include "../../../PDE-OPT/TOOLS/fe.hpp"

#include "Intrepid_HGRAD_LINE_Cn_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C2_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "ROL_Ptr.hpp"

template <class Real>
class Elliptic_Op : public PDE<Real> {
private:
  // Finite element basis information
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > basisPtr_;
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > basisPtrs_;
  // Cell cubature information
  ROL::Ptr<Intrepid::Cubature<Real> > cellCub_;
  // Cell node information
  ROL::Ptr<Intrepid::FieldContainer<Real> > volCellNodes_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellNodes_;
  std::vector<std::vector<std::vector<int> > > bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<FE<Real> > fe_vol_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int> > fidx_;
  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellDofValues_;

  Real epsilon_;
  bool useStateRiesz_;
  bool useControlRiesz_;

public:
  Elliptic_Op(Teuchos::ParameterList &parlist, Real epsilon) {
    // Finite element fields.
    int cubDegree  = parlist.sublist("Problem").get("Cubature Degree",4);
    basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >>();
    basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_);
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();
    Intrepid::DefaultCubatureFactory<Real> cubFactory;
    cellCub_ = cubFactory.create(cellType, cubDegree);
    // Problem data.
    useStateRiesz_   = parlist.sublist("Problem").get("Use State Riesz Map", true);      // use Riesz map for state variables?
    useControlRiesz_ = parlist.sublist("Problem").get("Use Control Riesz Map", true);    // use Riesz map for control variables?
    epsilon_ = epsilon; 
  }

   void residual(ROL::Ptr<Intrepid::FieldContainer<Real> > & res,
                const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    int c = u_coeff->dimension(0);
    int p = cellCub_->getNumPoints();
    int f = basisPtr_->getCardinality();
    int d = cellCub_->getDimension();
    // INITIALIZE RESIDUAL
    res = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // COMPUTE STIFFNESS TERM
    ROL::Ptr<Intrepid::FieldContainer<Real> > U_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_vol_->evaluateValue(U_eval, u_coeff);
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradU_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    
    ROL::Ptr<Intrepid::FieldContainer<Real> > laplace_term = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    Intrepid::FunctionSpaceTools::integrate<Real>(*laplace_term, *gradU_eval, *(fe_vol_->gradNdetJ()), Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*laplace_term,epsilon_);

    // ADD REACTION TERM TO RESIDUAL
    ROL::Ptr<Intrepid::FieldContainer<Real> > eye_term = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // Integrate U * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*eye_term,
                                                  *U_eval,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    
    Intrepid::RealSpaceTools<Real>::scale(*res,0.0);
    Intrepid::RealSpaceTools<Real>::add(*res,*laplace_term);
    Intrepid::RealSpaceTools<Real>::add(*res,*eye_term);

  }

  void Jacobian_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    int c = u_coeff->dimension(0);
    int p = cellCub_->getNumPoints();
    int f = basisPtr_->getCardinality();
    int d = cellCub_->getDimension();
    // INITILAIZE JACOBIAN
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    // COMPUTE STIFFNESS TERM
    ROL::Ptr<Intrepid::FieldContainer<Real> > laplace_term = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    Intrepid::FunctionSpaceTools::integrate<Real>(*laplace_term, *(fe_vol_->gradN()), *(fe_vol_->gradNdetJ()), Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*laplace_term,epsilon_);

    // ADD REACTION TERM TO JACOBIAN
    ROL::Ptr<Intrepid::FieldContainer<Real> > eye_term = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    // Integrate N * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*eye_term,
                                                  *(fe_vol_->N()),
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);

    Intrepid::RealSpaceTools<Real>::scale(*jac,0.0);
    Intrepid::RealSpaceTools<Real>::add(*jac,*laplace_term);
    Intrepid::RealSpaceTools<Real>::add(*jac,*eye_term);

  }   

  void Jacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    int c = u_coeff->dimension(0);
    int p = cellCub_->getNumPoints();
    int f = basisPtr_->getCardinality();
    int d = cellCub_->getDimension();
    // INITILAIZE JACOBIAN
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
  }


  void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Poisson::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Poisson::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Poisson::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Poisson::Hessian_22): Hessian is zero.");
  }

  void RieszMap_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & riesz) {
    // Optionally disable Riesz map ...
    if (!useStateRiesz_) {
      throw Exception::NotImplemented(">>> (PDE_Poisson::RieszMap_1): Not implemented.");
    }

    // ...otherwise ...

    int c = fe_vol_->N()->dimension(0);
    int f = fe_vol_->N()->dimension(1);
    // INITILAIZE JACOBIAN
    riesz = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    *riesz = *fe_vol_->stiffMat();
    Intrepid::RealSpaceTools<Real>::add(*riesz,*(fe_vol_->massMat()));
  }

  void RieszMap_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & riesz) {
    // Optionally disable Riesz map ...
    if (!useControlRiesz_) {
      throw Exception::NotImplemented(">>> (PDE_Poisson::RieszMap_2): Not implemented.");
    }

    // ...otherwise ...

    int c = fe_vol_->N()->dimension(0);
    int f = fe_vol_->N()->dimension(1);
    // INITILAIZE JACOBIAN
    riesz = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    *riesz = *fe_vol_->massMat();
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
    fe_vol_ = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtr_,cellCub_);
    // Set local boundary DOFs.
    fidx_ = fe_vol_->getBoundaryDofs();
    // Compute Dirichlet values at DOFs.
    int d = basisPtr_->getBaseCellTopology().getDimension();
    int numSidesets = bdryCellLocIds_.size();
    bdryCellDofValues_.resize(numSidesets);
  }

  const ROL::Ptr<FE<Real> > getFE(void) const {
    return fe_vol_;
  }

}; // Elliptic_Op

#endif
