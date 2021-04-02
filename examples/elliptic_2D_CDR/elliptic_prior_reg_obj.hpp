#ifndef ELLIPTIC_PRIOR_REGULARIZATION_OBJECTIVE_HPP
#define ELLIPTIC_PRIOR_REGULARIZATION_OBJECTIVE_HPP

#include "ROL_Objective_SimOpt.hpp"
#include "../../../PDE-OPT/TOOLS/qoi.hpp"

template <class Real>
class PDE_Reg_Op : public PDE<Real> {
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

  Real gamma_, alpha_;
  bool dirichlet_;

public:
  PDE_Reg_Op(Teuchos::ParameterList &parlist) {
    // Finite element fields.
    int cubDegree  = parlist.sublist("Problem").get("Cubature Degree",4);
    basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >>();
    basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_);
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();
    Intrepid::DefaultCubatureFactory<Real> cubFactory;
    cellCub_ = cubFactory.create(cellType, cubDegree);
    gamma_ = parlist.sublist("Problem").get("gamma",1.0);
    alpha_ = parlist.sublist("Problem").get("alpha",1.0);
    dirichlet_ = parlist.sublist("Problem").get("Dirichlet Prior",true);
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
    Intrepid::RealSpaceTools<Real>::scale(*laplace_term,gamma_);

    // ADD REACTION TERM TO RESIDUAL
    ROL::Ptr<Intrepid::FieldContainer<Real> > eye_term = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // Integrate U * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*eye_term,
                                                  *U_eval,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*eye_term,alpha_);
    
    Intrepid::RealSpaceTools<Real>::scale(*res,0.0);
    Intrepid::RealSpaceTools<Real>::add(*res,*laplace_term);
    Intrepid::RealSpaceTools<Real>::add(*res,*eye_term);

    if(dirichlet_)
      {
	// APPLY DIRICHLET CONDITIONS
	int numSideSets = bdryCellLocIds_.size();
	if (numSideSets > 0) {
	  for (int i = 0; i < numSideSets; ++i) {
	    int numLocalSideIds = bdryCellLocIds_[i].size();
	    for (int j = 0; j < numLocalSideIds; ++j) {
	      int numCellsSide = bdryCellLocIds_[i][j].size();
	      int numBdryDofs = fidx_[j].size();
	      for (int k = 0; k < numCellsSide; ++k) {
		int cidx = bdryCellLocIds_[i][j][k];
		for (int l = 0; l < numBdryDofs; ++l) {
		  (*res)(cidx,fidx_[j][l]) = (1.e-3)*( (*u_coeff)(cidx,fidx_[j][l]) - (*bdryCellDofValues_[i][j])(k,fidx_[j][l]) );
		}
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
    ROL::Ptr<Intrepid::FieldContainer<Real> > laplace_term = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    Intrepid::FunctionSpaceTools::integrate<Real>(*laplace_term, *(fe_vol_->gradN()), *(fe_vol_->gradNdetJ()), Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*laplace_term,gamma_);

    // ADD REACTION TERM TO JACOBIAN
    ROL::Ptr<Intrepid::FieldContainer<Real> > eye_term = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    // Integrate N * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*eye_term,
                                                  *(fe_vol_->N()),
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*eye_term,alpha_);

    Intrepid::RealSpaceTools<Real>::scale(*jac,0.0);
    Intrepid::RealSpaceTools<Real>::add(*jac,*laplace_term);
    Intrepid::RealSpaceTools<Real>::add(*jac,*eye_term);

    if(dirichlet_)
      {
	// APPLY DIRICHLET CONDITIONS
	int numSideSets = bdryCellLocIds_.size();
	if (numSideSets > 0) {
	  for (int i = 0; i < numSideSets; ++i) {
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
		  (*jac)(cidx,fidx_[j][l],fidx_[j][l]) = (1.e-3)*static_cast<Real>(1);
		}
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
    // GET DIMENSIONS
    int c = u_coeff->dimension(0);
    int p = cellCub_->getNumPoints();
    int f = basisPtr_->getCardinality();
    int d = cellCub_->getDimension();
    // INITILAIZE JACOBIAN
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    //throw Exception::Zero(">>> (PDE_Reg_Op::Jacobian_2): Jacobian is zero.");
  }

  void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Reg_Op::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Reg_Op::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Reg_Op::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Reg_Op::Hessian_22): Hessian is zero.");
  }

  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > getFields() {
    return basisPtrs_;
  }

  Real DirichletFunc(const std::vector<Real> & coords, int sideset, int locSideId) const {
    return 0;
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

}; // PDE_Reg_Op

template <class Real>
class PDE_Mass_Mat : public PDE<Real> {
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

public:
  PDE_Mass_Mat(Teuchos::ParameterList &parlist) {
    // Finite element fields.
    int cubDegree  = parlist.sublist("Problem").get("Cubature Degree",4);
    basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >>();
    basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_);
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();
    Intrepid::DefaultCubatureFactory<Real> cubFactory;
    cellCub_ = cubFactory.create(cellType, cubDegree);
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
    
    // Integrate U * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *U_eval,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
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
    // Integrate N * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  *(fe_vol_->N()),
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
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
    throw Exception::Zero(">>> (PDE_Reg_Op::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Reg_Op::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Reg_Op::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Reg_Op::Hessian_22): Hessian is zero.");
  }

  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > getFields() {
    return basisPtrs_;
  }

  Real DirichletFunc(const std::vector<Real> & coords, int sideset, int locSideId) const {
    return 0;
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

}; // PDE_Mass_Mat

template<class Real>
void Map_Index_to_Coords(int k, std::vector<Real> & coords, int d, int nx, int ny, Real nx_float, Real ny_float)
{
  coords[0] = static_cast<Real>( (k%(nx+1)) )*(1.0/nx_float);
  coords[1] = static_cast<Real>( std::floor( static_cast<Real>(k)/(nx_float+1.0) ) )*(1.0/ny_float);
}

template<class Real>
Real Mean_Eval(Real x, Real y)
{
  Real val = 16.0*x*(1.0-x)*y*(1.0-y);
  return val;
}

template<class Real>
void Set_Prior_Mean(HDSA::Ptr<Tpetra::MultiVector<> > & z_ptr, const Teuchos::RCP<Teuchos::ParameterList> & parlist)
{
  int nx = parlist->sublist("Geometry").get("NX",0);
  int ny = parlist->sublist("Geometry").get("NY",0);
  int dim = (nx+1)*(ny+1);
  Real nx_float = static_cast<Real>(nx);
  Real ny_float = static_cast<Real>(ny);
  std::vector<Real> coords = std::vector<Real>(2,0.0);
  for(int k = 0; k < dim; k++)
    {
      Map_Index_to_Coords<Real>(k,coords,dim,nx,ny,nx_float,ny_float);
      z_ptr->replaceGlobalValue(k,0,Mean_Eval<Real>(coords[0],coords[1]));
    }
}

template <class Real>
class Elliptic_Prior_Regularization_Objective : public ROL::Objective_SimOpt<Real> {

private:
  HDSA::Ptr<ROL::Constraint_SimOpt<Real> > con_;
  HDSA::Ptr<ROL::Constraint_SimOpt<Real> > mass_mat_con_;
  HDSA::Ptr<ROL::Vector<Real> > prior_mean_;

public:
  Elliptic_Prior_Regularization_Objective(const HDSA::Ptr<const HDSA::Comm<int> > & comm, const HDSA::Ptr<HDSA::ParameterList > & parlist, HDSA::Ptr<std::ostream> & outStream, bool construct_operators = false)
    : ROL::Objective_SimOpt<Real>()
  { 
    HDSA::Ptr<MeshManager<Real> > meshMgr = HDSA::makePtr<MeshManager_Rectangle<Real> >(*parlist);

    HDSA::Ptr<PDE<Real> > elliptic_pde = HDSA::makePtr<PDE_Reg_Op<Real> >(*parlist);
    con_ = HDSA::makePtr<Linear_PDE_Constraint<Real> >(elliptic_pde,meshMgr,comm->Get_Teuchos_Communicator(),*parlist,*outStream,true);

    HDSA::Ptr<PDE<Real> > mass_mat_pde = HDSA::makePtr<PDE_Mass_Mat<Real> >(*parlist);
    mass_mat_con_ = HDSA::makePtr<Linear_PDE_Constraint<Real> >(mass_mat_pde,meshMgr,comm->Get_Teuchos_Communicator(),*parlist,*outStream,true);

    // Cast the constraint and get the assembler.
    HDSA::Ptr<Linear_PDE_Constraint<Real> > pdecon = HDSA::dynamicPtrCast<Linear_PDE_Constraint<Real> >(con_);
    HDSA::Ptr<Assembler<Real> > assembler = pdecon->getAssembler();
    HDSA::Ptr<Tpetra::MultiVector<> > prior_mean_ptr  = assembler->createControlVector();   prior_mean_ptr->putScalar(0.0);
    prior_mean_  = HDSA::makePtr<PDE_PrimalOptVector<Real> >(prior_mean_ptr,elliptic_pde,assembler);

    Set_Prior_Mean<Real>(prior_mean_ptr,parlist);

    if(construct_operators)
      {
	HDSA::Ptr<ROL::Vector<Real> > zp = prior_mean_->clone();
	HDSA::Ptr<ROL::Vector<Real> > Kz = zp->clone();
	HDSA::Ptr<ROL::Vector<Real> > Mz = zp->clone();
	HDSA::Ptr<ROL::Vector<Real> > tmp1 = zp->clone();
	int n = zp->dimension();

	std::vector<std::vector<Real> > K;
	std::vector<std::vector<Real> > M;
	K.resize(n);
	M.resize(n);
	for(int i = 0; i < n; i++)
	  {
	    K[i].resize(n);
	    M[i].resize(n);
	  }

	Real tol = 1.e-8;
	for(int j = 0; j < n; j++)
	  {
	    std::cout << "Constructing the " << j+1 << "th out of " << n << " columns for the prior operators" << std::endl;
	    zp->set(*Kz->basis(j));
	    con_->applyJacobian_1(*Kz,*zp,*tmp1,*tmp1,tol);
	    mass_mat_con_->applyJacobian_1(*Mz,*zp,*tmp1,*tmp1,tol);
	    for(int i = 0; i < n; i++)
	      {
		K[i][j] = Kz->dot(*zp->basis(i));
		M[i][j] = Mz->dot(*zp->basis(i));
	      }
	  }
	std::cout << "Writing the operators to a file" << std::endl;
	std::string name;
	std::ofstream fout;
	name = "Elliptic_Operator.txt";
	fout.open(name);
	for(int i = 0; i < n; i++)
	  {
	    for(int j = 0; j < n; j++)
	      {
		fout << std::setprecision(16) << K[i][j] << "  ";
	      }
	    fout << "  " << std::endl;
	  }
	fout.close();
	name = "Mass_Mat.txt";
	fout.open(name);
	for(int i = 0; i < n; i++)
	  {
	    for(int j = 0; j < n; j++)
	      {
		fout << std::setprecision(16) << M[i][j] << "  ";
	      }
	    fout << "  " << std::endl;
	  }
	fout.close();
	name = "Prior_Mean.txt";
	fout.open(name);
	for(int i = 0; i < n; i++)
	  {
	    fout << std::setprecision(16) << prior_mean_->dot(*zp->basis(i)) << std::endl;
	  }
	fout.close();
      }
  }
  
  Real value(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<ROL::Vector<Real> > tmp = z.clone();
    tmp->set(z);
    tmp->axpy(-1.0,*prior_mean_);
    ROL::Ptr<ROL::Vector<Real> > Az = z.clone();
    con_->applyJacobian_1(*Az,*tmp,u,u,tol);
    ROL::Ptr<ROL::Vector<Real> > MinvAz = z.clone();
    mass_mat_con_->applyInverseJacobian_1(*MinvAz,*Az,u,u,tol);
    ROL::Ptr<ROL::Vector<Real> > AtMinvAz = z.clone();
    con_->applyAdjointJacobian_1(*AtMinvAz,*MinvAz,u,u,tol);
    Real val = 0.5*AtMinvAz->dot(*tmp);
    return val;
  }

  void gradient_1(ROL::Vector<Real> &g, const ROL::Vector<Real> &u,
                  const ROL::Vector<Real> &z, Real &tol ) {
    g.zero();
  }

  void gradient_2(ROL::Vector<Real> &g, const ROL::Vector<Real> &u,
                  const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<ROL::Vector<Real> > tmp = z.clone();
    tmp->set(z);
    tmp->axpy(-1.0,*prior_mean_);
    ROL::Ptr<ROL::Vector<Real> > Az = z.clone();
    con_->applyJacobian_1(*Az,*tmp,u,u,tol);
    ROL::Ptr<ROL::Vector<Real> > MinvAz = z.clone();
    mass_mat_con_->applyInverseJacobian_1(*MinvAz,*Az,u,u,tol);
    con_->applyAdjointJacobian_1(g,*MinvAz,u,u,tol);
  }

  void hessVec_11( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
             const ROL::Vector<Real> &u,  const ROL::Vector<Real> &z, Real &tol ) {
    hv.zero();
  }

  void hessVec_12( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    hv.zero();
  }

  void hessVec_21( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    hv.zero();
  }

  void hessVec_22( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
             const ROL::Vector<Real> &u,  const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<ROL::Vector<Real> > Av = z.clone();
    con_->applyJacobian_1(*Av,v,u,u,tol);
    ROL::Ptr<ROL::Vector<Real> > MinvAv = z.clone();
    mass_mat_con_->applyInverseJacobian_1(*MinvAv,*Av,u,u,tol);
    con_->applyAdjointJacobian_1(hv,*MinvAv,u,u,tol);
  } 

}; // Elliptic_Prior_Regularization_Objective

#endif
