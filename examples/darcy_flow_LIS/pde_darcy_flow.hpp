#ifndef PDE_DARCY_FLOW_HPP
#define PDE_DARCY_FLOW_HPP

#include "../../../PDE-OPT/TOOLS/pde.hpp"
#include "../../../PDE-OPT/TOOLS/fe.hpp"

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
    ROL::Ptr<Intrepid::FieldContainer<Real> > Prs_source
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > kappa
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    computePermeability(kappa);
    // COMPUTE PERMEABILITY TERM
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

    // Add source and control terms to residual
    int size = z_param->size();
    for (int i = 0; i < size; ++i) 
      {
	ROL::Ptr<Intrepid::FieldContainer<Real> > ctrl = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
	computeControlOperator(ctrl,(*z_param)[i],i);
	Intrepid::RealSpaceTools<Real>::add(*Prs_source, *ctrl);
      }
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
						  *Prs_source,
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
    ROL::Ptr<Intrepid::FieldContainer<Real> > kappa
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    computePermeability(kappa);
    // COMPUTE PERMEABILITY TERM
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
    throw Exception::Zero(">>> (PDE_darcy_flow::Jacobian_2): Jacobian is zero.");
  }

  void Jacobian_3(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
		  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // GET DIMENSIONS
    int c  = u_coeff->dimension(0);
    int p  = cellCub_->getNumPoints();
    int f = basisPtr_->getCardinality();

    // ADD CONTROL TERM TO RESIDUAL
    int size = z_param->size();
    jac.resize(size);

    for (int i = 0; i < size; ++i) 
      {
	ROL::Ptr<Intrepid::FieldContainer<Real> > ctrl = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
	computeControlOperator(ctrl,1.0,i);
	jac[i] = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);
	Intrepid::FunctionSpaceTools::integrate<Real>(*jac[i],
						      *ctrl,
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
  	      (*jac[i])(cidx,fidx_[j][l]) = static_cast<Real>(0);
  	    }
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
    throw Exception::Zero(">>> (PDE_darcy_flow::Hessian_12): Hessian is zero.");
  }

  void Hessian_13(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_darcy_flow::Hessian_13): Hessian is zero.");
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_darcy_flow::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_darcy_flow::Hessian_22): Hessian is zero.");
  }

  void Hessian_23(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_darcy_flow::Hessian_23): Hessian is zero.");
  }

  void Hessian_31(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_darcy_flow::Hessian_31): Hessian is zero.");
  }

  void Hessian_32(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_darcy_flow::Hessian_32): Hessian is zero.");
  }

  void Hessian_33(std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_darcy_flow::Hessian_33): Hessian is zero.");
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
  Real evaluateControlOperator(const std::vector<Real> &x, const int i) const {
    const std::vector<Real> xl = {0.2, 0.2, 0.5, 0.5, 0.8, 0.8};
    const std::vector<Real> yl = {0.5, 0.8, 0.5, 0.8, 0.5, 0.8};
    return 1000.0*std::exp(-1000.0*( (x[0]-xl[i])*(x[0]-xl[i]) + (x[1]-yl[i])*(x[1]-yl[i]) ) );
  }
  
  void computeControlOperator(ROL::Ptr<Intrepid::FieldContainer<Real> > &ctrl,
                              const Real z, const int I) const {
    // GET DIMENSIONS
    int c = fe_vol_->gradN()->dimension(0);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);
    std::vector<Real> pt(d), adv(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k) {
          pt[k] = (*fe_vol_->cubPts())(i,j,k);
        }
        // Compute control operator
        (*ctrl)(i,j) = -z*evaluateControlOperator(pt,I);
      }
    }
  }

  Real evaluateDirichlet(const std::vector<Real> & coords, int sideset, int locSideId) const {
    Real val = 0.0;
    if(coords[1] == 1.0)
      {
	val = 5.0;
      }
    else
      {
	val = 8.0;
      }
    return val;
  }

  Real evaluatePermeability(const std::vector<Real> &x) const {
    int T = 5;
    std::vector<Real> amp = std::vector<Real>(T);
    std::vector<Real> sdx = std::vector<Real>(T);
    std::vector<Real> sdy = std::vector<Real>(T);
    std::vector<Real> c = std::vector<Real>(T);
    std::vector<Real> h = std::vector<Real>(T);
    std::vector<Real> k= std::vector<Real>(T);

    amp[0] = -2.0; amp[1] = -2.2; amp[2] = 2.1; amp[3] = 2.0; amp[4] = 1.9;
    sdx[0] = 1.0; sdx[1] = 2.0; sdx[2] = 17.0; sdx[3] = 7.0; sdx[4] = 29.0;
    sdy[0] = 50.0; sdy[1] = 43.0; sdy[2] = 34.0; sdy[3] = 46.0; sdy[4] = 58.0;
    c[0] = 10.0; c[1] = 6.0; c[2] = 8.0; c[3] = 6.0; c[4] = 11.0;
    h[0] = .6; h[1] = .45; h[2] = .27; h[3] = .53; h[4] = .78;
    k[0] = .75; k[1] = .15; k[2] = .43; k[3] = .55; k[4] = .62;

    Real val = 0.0;
    for(int i = 0; i < T; i++)
      {
	val += amp[i]*std::exp(-1.0*( sdx[i]*(x[0]-h[i])*(x[0]-h[i]) + c[i]*(x[0]-h[i])*(x[1]-k[i]) + sdy[i]*(x[1]-k[i])*(x[1]-k[i]) ) );
      }
    
    val = std::pow(10.0,val);
    
    const std::vector<Real> param = PDE<Real>::getParameter();
    if((int)param.size()>0)
      {
	Real val_uncertain = 0.0;
	int offset = 0;
	int count = 0;
	Real basis_fun = 0.0;
	for(int j = 0; j < L+1; j++)
	  {
	    for(int i = 0; i < L+1; i++)
	      {
		basis_fun = std::exp(-10.0*static_cast<Real>(L+1)*(x[0]-uncertain_basis_grid[i])*(x[0]-uncertain_basis_grid[i]));
		basis_fun *= std::exp(-10.0*static_cast<Real>(L+1)*(x[1]-uncertain_basis_grid[j])*(x[1]-uncertain_basis_grid[j]));
		val_uncertain += param[offset + count]*basis_fun;
		count = count+1;
	      }
	  }
	val_uncertain = a*val_uncertain;
	val_uncertain = 1+val_uncertain;
	val = val*val_uncertain;
      }
    
    return val;
  }

  void computePermeability(ROL::Ptr<Intrepid::FieldContainer<Real> > &kappa) const {
    // GET DIMENSIONS
    int c = fe_vol_->gradN()->dimension(0);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);
    std::vector<Real> pt(d), adv(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k) {
          pt[k] = (*fe_vol_->cubPts())(i,j,k);
        }
        // Compute permeability kappa
        (*kappa)(i,j) = evaluatePermeability(pt);
      }
    }
  }


}; // PDE_darcy_flow

#endif
