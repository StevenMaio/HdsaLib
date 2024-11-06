#ifndef PDE_DARCY_FLOW_HPP
#define PDE_DARCY_FLOW_HPP

#include "../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/pde.hpp"
#include "../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/fe.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "ROL_Ptr.hpp"

// Defines the porous media flow PDE 
// -\nabla \cdot (\kappa * \nabla u) = z
// where u is the state, z is the source, and \kappa is the permeability

template <class Real>
class PDE_darcy_flow : public PDE<Real> {
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

  Real a;
  int L;
  std::vector<Real> uncertain_basis_grid;
 
public:
  PDE_darcy_flow(Teuchos::ParameterList &parlist) {
    // Finite element fields.
    int basisOrder = parlist.sublist("PDE Poisson").get("Basis Order",1);
    if (basisOrder == 1) {
      basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >>();
    }
    else if (basisOrder == 2) {
      basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real> >>();
    }
    basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_);
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();        // get the cell type from any basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                       // create cubature factory
    int cubDegree = parlist.sublist("PDE Poisson").get("Cubature Degree",2); // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create(cellType, cubDegree);                       // create default cubature
    
    L = parlist.sublist("Problem").get("Number of Uncertainty Basis Function", 10);
    a = parlist.sublist("Problem").get("Noise Level", .2);
    uncertain_basis_grid.resize(L+1);
    for(int i = 0; i < L+1; i++)
      {
	uncertain_basis_grid[i] = static_cast<Real>(i)/static_cast<Real>(L);
      }

  }

  void residual(ROL::Ptr<Intrepid::FieldContainer<Real> > & res,
                const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);
    // INITIALIZE RESIDUAL
    res = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);
    // COMPUTE PDE COEFFICIENTS
    ROL::Ptr<Intrepid::FieldContainer<Real> > rhs
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real>> kappa =
        ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_vol_->evaluateValue(kappa, z_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*kappa)(i,j) = std::exp((*kappa)(i,j));
      }
    }

    // Compute grad(U)
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradU_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    // Multiply kappa * grad(U)
    Intrepid::FieldContainer<Real> kappa_gradU(c, p, d);
    Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(kappa_gradU,
                                                               *kappa,
                                                               *gradU_eval);
    // Integrate (kappa * grad(U)) . grad(N)
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  kappa_gradU,
                                                  *(fe_vol_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, false);

    // Add source term to residual
    computeSource(rhs);
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
						  *rhs,
						  *(fe_vol_->NdetJ()),
						  Intrepid::COMP_CPP, true);

    // APPLY DIRICHLET CONDITIONS
    int numLocalSideIds = bdryCellLocIds_[0].size();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[0][j].size();
      int numBdryDofs = fidx_[j].size();
      for (int k = 0; k < numCellsSide; ++k) {
        int cidx = bdryCellLocIds_[0][j][k];
        for (int l = 0; l < numBdryDofs; ++l) {
          (*res)(cidx,fidx_[j][l])
            = (*u_coeff)(cidx,fidx_[j][l]) - (*bdryCellDofValues_[0][j])(k,fidx_[j][l]);
        }
      }
    }
  }

  void Jacobian_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);
    // INITILAIZE JACOBIAN
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    // COMPUTE PDE COEFFICIENTS
    ROL::Ptr<Intrepid::FieldContainer<Real>> kappa =
        ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_vol_->evaluateValue(kappa, z_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*kappa)(i,j) = std::exp((*kappa)(i,j));
      }
    }
    // Multiply kappa * grad(N)
    Intrepid::FieldContainer<Real> kappa_gradN(c, f, p, d);
    Intrepid::FunctionSpaceTools::tensorMultiplyDataField<Real>(kappa_gradN,
                                                                *kappa,
                                                                *(fe_vol_->gradN()));
    // Integrate (kappa * grad(N)) . grad(N)
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  kappa_gradN,
                                                  *(fe_vol_->gradNdetJ()),
                                                  Intrepid::COMP_CPP, false);
  
    // APPLY DIRICHLET CONDITIONS
    int numLocalSideIds = bdryCellLocIds_[0].size();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[0][j].size();
      int numBdryDofs = fidx_[j].size();
      for (int k = 0; k < numCellsSide; ++k) {
        int cidx = bdryCellLocIds_[0][j][k];
        for (int l = 0; l < numBdryDofs; ++l) {
          //std::cout << "\n   j=" << j << "  l=" << l << "  " << fidx[j][l];
          for (int m = 0; m < f; ++m) {
            (*jac)(cidx,fidx_[j][l],m) = static_cast<Real>(0);
          }
          (*jac)(cidx,fidx_[j][l],fidx_[j][l]) = static_cast<Real>(1);
        }
      }
    }
  }
  
  void Jacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);
    // INITILAIZE JACOBIAN
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);

    // Comptue grad(Prs)
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradPrs_eval;
    gradPrs_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    fe_vol_->evaluateGradient(gradPrs_eval, u_coeff);

    ROL::Ptr<Intrepid::FieldContainer<Real> > kappa =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    fe_vol_->evaluateValue(kappa, z_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*kappa)(i,j) = std::exp((*kappa)(i,j));
      }
    }

    Intrepid::FieldContainer<Real> gradPrs_gradN_eval(c, f, p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(gradPrs_gradN_eval,*gradPrs_eval,*(fe_vol_->gradN()));
    Intrepid::FieldContainer<Real> gradPrs_gradN_kappa_eval(c, f, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(gradPrs_gradN_kappa_eval,
								*kappa,
								gradPrs_gradN_eval);
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  gradPrs_gradN_kappa_eval,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);

    // APPLY DIRICHLET CONDITIONS
    int numLocalSideIds = bdryCellLocIds_[0].size();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[0][j].size();
      int numBdryDofs = fidx_[j].size();
      for (int k = 0; k < numCellsSide; ++k) {
        int cidx = bdryCellLocIds_[0][j][k];
        for (int l = 0; l < numBdryDofs; ++l) {
          //std::cout << "\n   j=" << j << "  l=" << l << "  " << fidx[j][l];
          for (int m = 0; m < f; ++m) {
            (*jac)(cidx,fidx_[j][l],m) = static_cast<Real>(0);
          }
          (*jac)(cidx,fidx_[j][l],fidx_[j][l]) = static_cast<Real>(0);
        }
      }
    }

  }

 void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_darcy_flow::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);

    ROL::Ptr<Intrepid::FieldContainer<Real> > lambda_coeff = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);
    for (int i = 0; i < c; i++)
    {
      for (int j = 0; j < f; j++)
      {
        (*lambda_coeff)(i,j) = (*l_coeff)(i,j);
      }
    }
    // APPLY DIRICHLET CONDITIONS
    int numLocalSideIds = bdryCellLocIds_[0].size();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[0][j].size();
      int numBdryDofs = fidx_[j].size();
      for (int k = 0; k < numCellsSide; ++k) {
        int cidx = bdryCellLocIds_[0][j][k];
        for (int l = 0; l < numBdryDofs; ++l) {
          (*lambda_coeff)(cidx,fidx_[j][l]) = 0.0;
        }
      }
    }

    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradL_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    fe_vol_->evaluateGradient(gradL_eval, lambda_coeff);

     ROL::Ptr<Intrepid::FieldContainer<Real> > kappa =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    fe_vol_->evaluateValue(kappa, z_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*kappa)(i,j) = std::exp((*kappa)(i,j));
      }
    }

    Intrepid::FieldContainer<Real> gradL_gradN_eval(c, f, p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(gradL_gradN_eval,*gradL_eval,*(fe_vol_->gradN()));
    Intrepid::FieldContainer<Real> gradL_gradN_kappa_eval(c, f, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(gradL_gradN_kappa_eval,
								*kappa,
								gradL_gradN_eval);
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  *(fe_vol_->NdetJ()),
                                                  gradL_gradN_kappa_eval,
                                                  Intrepid::COMP_CPP, false);

  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);

    ROL::Ptr<Intrepid::FieldContainer<Real> > lambda_coeff = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);
    for (int i = 0; i < c; i++)
    {
      for (int j = 0; j < f; j++)
      {
        (*lambda_coeff)(i,j) = (*l_coeff)(i,j);
      }
    }
    // APPLY DIRICHLET CONDITIONS
    int numLocalSideIds = bdryCellLocIds_[0].size();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[0][j].size();
      int numBdryDofs = fidx_[j].size();
      for (int k = 0; k < numCellsSide; ++k) {
        int cidx = bdryCellLocIds_[0][j][k];
        for (int l = 0; l < numBdryDofs; ++l) {
          (*lambda_coeff)(cidx,fidx_[j][l]) = 0.0;
        }
      }
    }

    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradL_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    fe_vol_->evaluateGradient(gradL_eval, lambda_coeff);

     ROL::Ptr<Intrepid::FieldContainer<Real> > kappa =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    fe_vol_->evaluateValue(kappa, z_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*kappa)(i,j) = std::exp((*kappa)(i,j));
      }
    }

    Intrepid::FieldContainer<Real> gradL_gradN_eval(c, f, p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(gradL_gradN_eval,*gradL_eval,*(fe_vol_->gradN()));
    Intrepid::FieldContainer<Real> gradL_gradN_kappa_eval(c, f, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(gradL_gradN_kappa_eval,
								*kappa,
								gradL_gradN_eval);
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  gradL_gradN_kappa_eval,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);

    ROL::Ptr<Intrepid::FieldContainer<Real> > lambda_coeff = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);
    for (int i = 0; i < c; i++)
    {
      for (int j = 0; j < f; j++)
      {
        (*lambda_coeff)(i,j) = (*l_coeff)(i,j);
      }
    }
    // APPLY DIRICHLET CONDITIONS
    int numLocalSideIds = bdryCellLocIds_[0].size();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[0][j].size();
      int numBdryDofs = fidx_[j].size();
      for (int k = 0; k < numCellsSide; ++k) {
        int cidx = bdryCellLocIds_[0][j][k];
        for (int l = 0; l < numBdryDofs; ++l) {
          (*lambda_coeff)(cidx,fidx_[j][l]) = 0.0;
        }
      }
    }

    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradL_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    fe_vol_->evaluateGradient(gradL_eval, lambda_coeff);

    // Comptue grad(Prs)
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradPrs_eval;
    gradPrs_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    fe_vol_->evaluateGradient(gradPrs_eval, u_coeff);

    ROL::Ptr<Intrepid::FieldContainer<Real> > kappa =
      ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    fe_vol_->evaluateValue(kappa, z_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*kappa)(i,j) = std::exp((*kappa)(i,j));
      }
    }

    Intrepid::FieldContainer<Real> gradPrs_gradL_eval(c, p);
    Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(gradPrs_gradL_eval,*gradPrs_eval,*gradL_eval);
    Intrepid::FieldContainer<Real> gradPrs_gradL_kappa_eval(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(gradPrs_gradL_kappa_eval,
								*kappa,
								gradPrs_gradL_eval);

    Intrepid::FieldContainer<Real> gradPrs_gradL_kappa_N_eval(c, p, f);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(gradPrs_gradL_kappa_N_eval,
								gradPrs_gradL_kappa_eval,
                *(fe_vol_->N()));

    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  gradPrs_gradL_kappa_N_eval,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

  void RieszMap_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & riesz) {
    // GET DIMENSIONS
    int c = fe_vol_->N()->dimension(0);
    int f = fe_vol_->N()->dimension(1);
    // INITIALIZE RIESZ
    riesz = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    *riesz = *fe_vol_->stiffMat();
    Intrepid::RealSpaceTools<Real>::add(*riesz,*(fe_vol_->massMat()));
  }

  void RieszMap_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & riesz) {
    // GET DIMENSIONS
    int c = fe_vol_->N()->dimension(0);
    int f = fe_vol_->N()->dimension(1);
    // INITIALIZE RIESZ
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
            (*bdryCellDofValues_[i][j])(k, l) = evaluateDirichlet(dofpoint, i, j);
          }
        }
      }
    }

  }

  const ROL::Ptr<FE<Real> > getFE(void) const {
    return fe_vol_;
  }

private:

private:
  void computeSource(ROL::Ptr<Intrepid::FieldContainer<Real> > &rhs) const {
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

  Real evaluateRHS(const std::vector<Real> & x) const {
    Real val = 0.0;
    return val;
  }

  Real evaluateDirichlet(const std::vector<Real> & coords, int sideset, int locSideId) const {
    Real val = 0.0;
    if(coords[1] == 1.0)
      {
	      val = 1.0;
      }
    else
      {
	      val = 0.0;
      }
    return val;
  }



}; // PDE_darcy_flow

#endif
