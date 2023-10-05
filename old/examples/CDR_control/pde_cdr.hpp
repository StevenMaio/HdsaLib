#ifndef PDE_CDR_HPP
#define PDE_CDR_HPP

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
class PDE_CDR : public PDE<Real> {
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

  bool useStateRiesz_;
  bool useControlRiesz_;

  Real diff_coeff_, adv_coeff_, cubic_react_coeff_, linear_react_coeff_;

  Real DirichletFunc(const std::vector<Real> & coords, int sideset, int locSideId) const {
    return 0;
  }

  Real evaluateRHS(const std::vector<Real> &x) const {
    return 0.0;
  }

  Real evaluatexVelocityField(const std::vector<Real> &x) const {
    Real pi(M_PI);
    return std::cos(2.0*pi*x[0]);
  }

  Real evaluateyVelocityField(const std::vector<Real> &x) const {
    Real pi(M_PI);
    return 1.0 + std::pow(std::cos(2.0*pi*x[1]),2.0);
  }

  void computeRHS(ROL::Ptr<Intrepid::FieldContainer<Real> > &rhs) const {
    // GET DIMENSIONS
    int c = fe_vol_->gradN()->dimension(0);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k) {
          pt[k] = (*fe_vol_->cubPts())(i,j,k);
        }
        // Compute forcing term f
        (*rhs)(i,j) = -evaluateRHS(pt);
      }
    }
  }

  void computeVelocityField(ROL::Ptr<Intrepid::FieldContainer<Real> > & V) const {
    // GET DIMENSIONS
    int c = fe_vol_->gradN()->dimension(0);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k) {
          pt[k] = (*fe_vol_->cubPts())(i,j,k);
        }
        // Compute velocity field
        (*V)(i,j,0) = evaluatexVelocityField(pt);
	(*V)(i,j,1) = evaluateyVelocityField(pt);
      }
    }
  }

public:
  PDE_CDR(Teuchos::ParameterList &parlist) {
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

    diff_coeff_ = parlist.sublist("Problem").get("Diffusion Coefficient", 1.0);
    adv_coeff_ = parlist.sublist("Problem").get("Advection Coefficient", 1.0);
    cubic_react_coeff_ = parlist.sublist("Problem").get("Cubic Reaction Coefficient", 1.0);
    linear_react_coeff_ = parlist.sublist("Problem").get("Linear Reaction Coefficient", 1.0);
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
    
    ROL::Ptr<Intrepid::FieldContainer<Real> > diff_term = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    Intrepid::FunctionSpaceTools::integrate<Real>(*diff_term, *gradU_eval, *(fe_vol_->gradNdetJ()), Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*diff_term,diff_coeff_);

    // ADD ADVECTION TERM TO RESIDUAL
    ROL::Ptr<Intrepid::FieldContainer<Real> > V
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    computeVelocityField(V);
    // Multiply V . grad(U)
    Intrepid::FieldContainer<Real> V_gradU(c, p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(V_gradU,
                                                            *V,
                                                            *gradU_eval);
    ROL::Ptr<Intrepid::FieldContainer<Real> > adv_term = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // Integrate (V . grad(U)) * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*adv_term,
                                                  V_gradU,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*adv_term,adv_coeff_);

    // ADD REACTION TERM TO RESIDUAL
    ROL::Ptr<Intrepid::FieldContainer<Real> > react_term = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    ROL::Ptr<Intrepid::FieldContainer<Real> > U_eval3 = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for(int i = 0; i < c; i++)
      {
	for(int j = 0; j < p; j++)
	  {
	    (*U_eval3)(i,j) = cubic_react_coeff_*std::pow((*U_eval)(i,j),3.0) + linear_react_coeff_*(*U_eval)(i,j);
	  } 
      }
    // Integrate U * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*react_term,
                                                  *U_eval3,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*react_term,-1.0);

    // COMPUTE RHS
    ROL::Ptr<Intrepid::FieldContainer<Real> > rhs
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    computeRHS(rhs);
    ROL::Ptr<Intrepid::FieldContainer<Real> > rhs_term = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    Intrepid::FunctionSpaceTools::integrate<Real>(*rhs_term, *rhs, *(fe_vol_->NdetJ()), Intrepid::COMP_CPP, false);

    // ADD CONTROL TERM TO RESIDUAL
    ROL::Ptr<Intrepid::FieldContainer<Real> > valZ_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_vol_->evaluateValue(valZ_eval, z_coeff);
    Intrepid::RealSpaceTools<Real>::scale(*valZ_eval,static_cast<Real>(-1));
    ROL::Ptr<Intrepid::FieldContainer<Real> > source_term = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    Intrepid::FunctionSpaceTools::integrate<Real>(*source_term, *valZ_eval, *(fe_vol_->NdetJ()), Intrepid::COMP_CPP, false);
    
    Intrepid::RealSpaceTools<Real>::scale(*res,0.0);
    Intrepid::RealSpaceTools<Real>::add(*res,*diff_term);
    Intrepid::RealSpaceTools<Real>::add(*res,*adv_term);
    Intrepid::RealSpaceTools<Real>::add(*res,*react_term);
    Intrepid::RealSpaceTools<Real>::add(*res,*rhs_term);
    Intrepid::RealSpaceTools<Real>::add(*res,*source_term);

    // APPLY DIRICHLET CONDITIONS
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < 1; ++i) {
        int numLocalSideIds = bdryCellLocIds_[i].size();
        for (int j = 0; j < numLocalSideIds; ++j) {
          int numCellsSide = bdryCellLocIds_[i][j].size();
          int numBdryDofs = fidx_[j].size();
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[i][j][k];
            for (int l = 0; l < numBdryDofs; ++l) {
              (*res)(cidx,fidx_[j][l]) = (*u_coeff)(cidx,fidx_[j][l]) - (*bdryCellDofValues_[i][j])(k,fidx_[j][l]);
            }
          }
        }
      }
    }
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
    ROL::Ptr<Intrepid::FieldContainer<Real> > U_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_vol_->evaluateValue(U_eval, u_coeff);

    ROL::Ptr<Intrepid::FieldContainer<Real> > diff_term = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    Intrepid::FunctionSpaceTools::integrate<Real>(*diff_term, *(fe_vol_->gradN()), *(fe_vol_->gradNdetJ()), Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*diff_term,diff_coeff_);

    // ADD ADVECTION TERM TO JACOBIAN
    // Multiply V . grad(N)
    ROL::Ptr<Intrepid::FieldContainer<Real> > V
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    computeVelocityField(V);
    Intrepid::FieldContainer<Real> V_gradN(c, f, p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(V_gradN,
                                                             *V,
                                                             *(fe_vol_->gradN()));
    // Integrate (V . grad(U)) * N
    ROL::Ptr<Intrepid::FieldContainer<Real> > adv_term = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    Intrepid::FunctionSpaceTools::integrate<Real>(*adv_term,
                                                  *(fe_vol_->NdetJ()),
                                                  V_gradN,
                                                  Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*adv_term,adv_coeff_);

    // ADD REACTION TERM TO RESIDUAL
    ROL::Ptr<Intrepid::FieldContainer<Real> > react_term = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    ROL::Ptr<Intrepid::FieldContainer<Real> > U_eval3_diff = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for(int i = 0; i < c; i++)
      {
	for(int j = 0; j < p; j++)
	  {
	    (*U_eval3_diff)(i,j) = cubic_react_coeff_*3.0*std::pow((*U_eval)(i,j),2.0) + linear_react_coeff_;
	  } 
      }
    ROL::Ptr<Intrepid::FieldContainer<Real> > U_eval3_diff_N = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*U_eval3_diff_N,
                                                                *U_eval3_diff,
                                                                *(fe_vol_->N()));
    // Integrate U * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*react_term,
                                                  *U_eval3_diff_N,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*react_term,-1.0);


    Intrepid::RealSpaceTools<Real>::scale(*jac,0.0);
    Intrepid::RealSpaceTools<Real>::add(*jac,*diff_term);
    Intrepid::RealSpaceTools<Real>::add(*jac,*adv_term);
    Intrepid::RealSpaceTools<Real>::add(*jac,*react_term);

    // APPLY DIRICHLET CONDITIONS
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < 1; ++i) {
        int numLocalSideIds = bdryCellLocIds_[i].size();
        for (int j = 0; j < numLocalSideIds; ++j) {
          int numCellsSide = bdryCellLocIds_[i][j].size();
          int numBdryDofs = fidx_[j].size();
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[i][j][k];
            for (int l = 0; l < numBdryDofs; ++l) {
              for (int m = 0; m < f; ++m) {
                (*jac)(cidx,fidx_[j][l],m) = static_cast<Real>(0);
              }
              (*jac)(cidx,fidx_[j][l],fidx_[j][l]) = static_cast<Real>(1);
            }
          }
        }
      }
    }
  }

  void Jacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    if ( z_coeff != ROL::nullPtr ) {
      // GET DIMENSIONS
      int c = u_coeff->dimension(0);
      int f = basisPtr_->getCardinality();
      // INITIALIZE JACOBIAN
      jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);

      Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
						    *(fe_vol_->N()),
						    *(fe_vol_->NdetJ()),
						    Intrepid::COMP_CPP,false);
      
      Intrepid::RealSpaceTools<Real>::scale(*jac,-1.0);
      
      // APPLY DIRICHLET CONDITIONS
      int numSideSets = bdryCellLocIds_.size();
      if (numSideSets > 0) {
        for (int i = 0; i < 1; ++i) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m = 0; m < f; ++m) {
                  (*jac)(cidx,fidx_[j][l],m) = static_cast<Real>(0);
                }
              }
            }
          }
        }
      }
    }
    else {
      throw Exception::Zero(">>> (PDE_CDR::Jacobian_2): Jacobian is zero.");
    }
  }

  void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    int c = u_coeff->dimension(0);
    int p = cellCub_->getNumPoints();
    int f = basisPtr_->getCardinality();
    int d = cellCub_->getDimension();
    // INITILAIZE JACOBIAN
    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);

    ROL::Ptr<Intrepid::FieldContainer<Real> > L = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    for(int i = 0; i < c; i++)
      {
	for(int j = 0; j < f; j++)
	  {
	    (*L)(i,j) = (*l_coeff)(i,j);
	  }
      }

    // APPLY DIRICHLET CONDITIONS
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < 1; ++i) {
        int numLocalSideIds = bdryCellLocIds_[i].size();
        for (int j = 0; j < numLocalSideIds; ++j) {
          int numCellsSide = bdryCellLocIds_[i][j].size();
          int numBdryDofs = fidx_[j].size();
          for (int k = 0; k < numCellsSide; ++k) {
            int cidx = bdryCellLocIds_[i][j][k];
            for (int l = 0; l < numBdryDofs; ++l) {
              (*L)(cidx,fidx_[j][l]) = 0.0;
            }
          }
        }
      }
    }
   
    ROL::Ptr<Intrepid::FieldContainer<Real> > U_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_vol_->evaluateValue(U_eval, u_coeff);

    ROL::Ptr<Intrepid::FieldContainer<Real> > L_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_vol_->evaluateValue(L_eval, L);

    // ADD REACTION TERM TO RESIDUAL
    ROL::Ptr<Intrepid::FieldContainer<Real> > react_term = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    ROL::Ptr<Intrepid::FieldContainer<Real> > U_eval3_diff2 = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for(int i = 0; i < c; i++)
      {
	for(int j = 0; j < p; j++)
	  {
	    (*U_eval3_diff2)(i,j) = cubic_react_coeff_*6.0*(*U_eval)(i,j)*(*L_eval)(i,j);
	  } 
      }
    ROL::Ptr<Intrepid::FieldContainer<Real> > U_eval3_diff2_N = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*U_eval3_diff2_N,
                                                                *U_eval3_diff2,
                                                                *(fe_vol_->N()));
    // Integrate U * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*react_term,
                                                  *U_eval3_diff2_N,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*react_term,-1.0);

    Intrepid::RealSpaceTools<Real>::add(*hess,*react_term);
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_CDR::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_CDR::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_CDR::Hessian_22): Hessian is zero.");
  }

  void RieszMap_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & riesz) {
    // Optionally disable Riesz map ...
    if (!useStateRiesz_) {
      throw Exception::NotImplemented(">>> (PDE_CDR::RieszMap_1): Not implemented.");
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
      throw Exception::NotImplemented(">>> (PDE_CDR::RieszMap_2): Not implemented.");
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
    for (int i=0; i<numSidesets; ++i) {
      int numLocSides = bdryCellLocIds_[i].size();
      bdryCellDofValues_[i].resize(numLocSides);
      for (int j=0; j<numLocSides; ++j) {
        int c = bdryCellLocIds_[i][j].size();
        int f = basisPtr_->getCardinality();
        bdryCellDofValues_[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
        ROL::Ptr<Intrepid::FieldContainer<Real> > coords =
          ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, d);
        if (c > 0) {
          fe_vol_->computeDofCoords(coords, bdryCellNodes_[i][j]);
        }
        for (int k=0; k<c; ++k) {
          for (int l=0; l<f; ++l) {
            std::vector<Real> dofpoint(d);
            for (int m=0; m<d; ++m) {
              dofpoint[m] = (*coords)(k, l, m);
            }
            (*bdryCellDofValues_[i][j])(k, l) = DirichletFunc(dofpoint, i, j);
          }
        }
      }
    }
  }

  const ROL::Ptr<FE<Real> > getFE(void) const {
    return fe_vol_;
  }

}; // PDE_CDR

#endif
