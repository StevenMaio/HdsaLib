#ifndef ELLIPTIC_PRIOR_REGULARIZATION_OBJECTIVE_HPP
#define ELLIPTIC_PRIOR_REGULARIZATION_OBJECTIVE_HPP

#include "ROL_Objective_SimOpt.hpp"
#include "../../../PDE-OPT/TOOLS/linearpdeconstraint.hpp"

template <class Real>
class PDE_Reg_Op : public PDE<Real> {
private:
  // Finite element basis information
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > basisPtr_;
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > basisPtrs_;
  // Cell cubature information
  ROL::Ptr<Intrepid::Cubature<Real> > cellCub_;
  ROL::Ptr<Intrepid::Cubature<Real> > bdryCub_;
  // Cell node information
  ROL::Ptr<Intrepid::FieldContainer<Real> > volCellNodes_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellNodes_;
  std::vector<std::vector<std::vector<int> > > bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<FE<Real> > fe_;
  std::vector<std::vector<ROL::Ptr<FE<Real> > > > feBdry_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int> > fidx_;

  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellDofValues_;

  // Field pattern, offsets, etc.
  std::vector<std::vector<int> > fieldPattern_;  // local Field/DOF pattern; set from DOF manager 
  int numFields_;                                // number of fields (equations in the PDE)
  int numDofs_;                                  // total number of degrees of freedom for all (local) fields
  std::vector<int> offset_;                      // for each field, a counting offset
  std::vector<int> numFieldDofs_;                // for each field, number of degrees of freedom

  ROL::Ptr<FieldHelper<Real> > fieldHelper_;

  Real gamma_, alpha_;
  bool dirichlet_;

public:
  PDE_Reg_Op(Teuchos::ParameterList &parlist) {
    // Finite element fields -- NOT DIMENSION INDEPENDENT!
    basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >>();
    // Volume quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();         // get the cell type from any basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                           // create cubature factory
    int cubDegree = parlist.sublist("Problem").get("Cubature Degree", 2);        // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create(cellType, cubDegree);                           // create default cubature
    // Boundary quadrature rules.
    int d = cellType.getDimension();
    shards::CellTopology bdryCellType = cellType.getCellTopologyData(d-1, 0);
    int bdryCubDegree = parlist.sublist("Problem").get("Boundary Cubature Degree",2); // set cubature degree, e.g., 2
    bdryCub_ = cubFactory.create(bdryCellType, bdryCubDegree);
    basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_); basisPtrs_.push_back(basisPtr_); basisPtrs_.push_back(basisPtr_);

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

    // Initialize residuals.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > R(3);
    R[0]   = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);
    R[1] = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);
    R[2] = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);

    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);

    // COMPUTE STIFFNESS TERM
    ROL::Ptr<Intrepid::FieldContainer<Real> > U_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(U_eval, U[0]);
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradU_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    fe_->evaluateGradient(gradU_eval, U[0]);
    
    ROL::Ptr<Intrepid::FieldContainer<Real> > laplace_term = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    Intrepid::FunctionSpaceTools::integrate<Real>(*laplace_term, *gradU_eval, *(fe_->gradNdetJ()), Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*laplace_term,gamma_);

    // ADD REACTION TERM TO RESIDUAL
    ROL::Ptr<Intrepid::FieldContainer<Real> > eye_term = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // Integrate U * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*eye_term,
                                                  *U_eval,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*eye_term,alpha_);
    
    Intrepid::RealSpaceTools<Real>::scale(*R[0],0.0);
    Intrepid::RealSpaceTools<Real>::add(*R[0],*laplace_term);
    Intrepid::RealSpaceTools<Real>::add(*R[0],*eye_term);

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
		  (*R[0])(cidx,fidx_[j][l]) = (1.e-3)*( (*U[0])(cidx,fidx_[j][l]) - (*bdryCellDofValues_[i][j])(k,fidx_[j][l]) );
		}
	      }
	    }
	  }
	}
      }

    // Combine the residuals.
    fieldHelper_->combineFieldCoeff(res, R);

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
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > J(3);
    for (int i = 0; i < 3; ++i) {
      J[i].resize(3);
      for (int j = 0; j < 3; j++)
	{
	  J[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
	}
    }

    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);

    // COMPUTE STIFFNESS TERM
    ROL::Ptr<Intrepid::FieldContainer<Real> > laplace_term = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    Intrepid::FunctionSpaceTools::integrate<Real>(*laplace_term, *(fe_->gradN()), *(fe_->gradNdetJ()), Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*laplace_term,gamma_);

    // ADD REACTION TERM TO JACOBIAN
    ROL::Ptr<Intrepid::FieldContainer<Real> > eye_term = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    // Integrate N * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*eye_term,
                                                  *(fe_->N()),
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*eye_term,alpha_);

    Intrepid::RealSpaceTools<Real>::scale(*J[0][0],0.0);
    Intrepid::RealSpaceTools<Real>::add(*J[0][0],*laplace_term);
    Intrepid::RealSpaceTools<Real>::add(*J[0][0],*eye_term);

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
		    (*J[0][0])(cidx,fidx_[j][l],m) = static_cast<Real>(0);
		  }
		  (*J[0][0])(cidx,fidx_[j][l],fidx_[j][l]) = (1.e-3)*static_cast<Real>(1);
		}
	      }
	    }
	  }
	}
      }

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);

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
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > J(3);
    for (int i = 0; i < 3; ++i) {
      J[i].resize(3);
      for (int j = 0; j < 3; j++)
	{
	  J[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
	}
    }  
    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);
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

  void setFieldPattern(const std::vector<std::vector<int> > & fieldPattern) {
    fieldPattern_ = fieldPattern;
    fieldHelper_ = ROL::makePtr<FieldHelper<Real>>(numFields_, numDofs_, numFieldDofs_, fieldPattern_);
  }

 void setCellNodes(const ROL::Ptr<Intrepid::FieldContainer<Real> > &volCellNodes,
                    const std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int> > > &bdryCellLocIds) {
    volCellNodes_ = volCellNodes;
    bdryCellNodes_ = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    // Finite element definition.
    fe_ = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtr_,cellCub_);
    // Set local boundary DOFs.
    fidx_ = fe_->getBoundaryDofs();
    // Compute Dirichlet values at DOFs.
    int d = basisPtr_->getBaseCellTopology().getDimension();
    int numSidesets = bdryCellLocIds_.size();
    bdryCellDofValues_.resize(numSidesets);
    feBdry_.resize(numSidesets); 
    for (int i=0; i<numSidesets; ++i) {
      int numLocSides = bdryCellLocIds_[i].size();
      bdryCellDofValues_[i].resize(numLocSides);
      feBdry_[i].resize(numLocSides);
      for (int j=0; j<numLocSides; ++j) {
        int c = bdryCellLocIds_[i][j].size();
        int f = basisPtr_->getCardinality();
	if(c > 0)
	  {
	    feBdry_[i][j] = ROL::makePtr<FE<Real> >(bdryCellNodes_[i][j],basisPtr_,bdryCub_,j);
	  }
	bdryCellDofValues_[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
        ROL::Ptr<Intrepid::FieldContainer<Real> > coords =
          ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, d);
        if (c > 0) {
          fe_->computeDofCoords(coords, bdryCellNodes_[i][j]);
        }
        for (int k=0; k<c; ++k) {
          for (int l=0; l<f; ++l) {
            std::vector<Real> dofpoint(d);
            for (int m=0; m<d; ++m) {
              dofpoint[m] = (*coords)(k, l, m);
            }
            (*bdryCellDofValues_[i][j])(k, l) = 0.0;
          }
        }
      }
    }

  }

  ROL::Ptr<Intrepid::FieldContainer<Real> > getBoundaryCoeff(const Intrepid::FieldContainer<Real> & cell_coeff, int sideSet, int cell) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[sideSet][cell];
    const int numCellsSide = bdryCellLocId.size();
    const int f = basisPtr_->getCardinality();
    
    ROL::Ptr<Intrepid::FieldContainer<Real > > bdry_coeff = 
      ROL::makePtr<Intrepid::FieldContainer<Real > >(numCellsSide, f);
    for (int i = 0; i < numCellsSide; ++i) {
      for (int j = 0; j < f; ++j) {
	(*bdry_coeff)(i, j) = cell_coeff(bdryCellLocId[i], j);
      }
    }
    return bdry_coeff;
  }

  const ROL::Ptr<FE<Real> > getFE(void) const {
    return fe_;
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
  ROL::Ptr<Intrepid::Cubature<Real> > bdryCub_;
  // Cell node information
  ROL::Ptr<Intrepid::FieldContainer<Real> > volCellNodes_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellNodes_;
  std::vector<std::vector<std::vector<int> > > bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<FE<Real> > fe_;
  std::vector<std::vector<ROL::Ptr<FE<Real> > > > feBdry_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int> > fidx_;

  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellDofValues_;

  // Field pattern, offsets, etc.
  std::vector<std::vector<int> > fieldPattern_;  // local Field/DOF pattern; set from DOF manager 
  int numFields_;                                // number of fields (equations in the PDE)
  int numDofs_;                                  // total number of degrees of freedom for all (local) fields
  std::vector<int> offset_;                      // for each field, a counting offset
  std::vector<int> numFieldDofs_;                // for each field, number of degrees of freedom

  ROL::Ptr<FieldHelper<Real> > fieldHelper_;

public:
  PDE_Mass_Mat(Teuchos::ParameterList &parlist) {
 // Finite element fields -- NOT DIMENSION INDEPENDENT!
    basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >>();
    // Volume quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();         // get the cell type from any basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                           // create cubature factory
    int cubDegree = parlist.sublist("Problem").get("Cubature Degree", 2);        // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create(cellType, cubDegree);                           // create default cubature
    // Boundary quadrature rules.
    int d = cellType.getDimension();
    shards::CellTopology bdryCellType = cellType.getCellTopologyData(d-1, 0);
    int bdryCubDegree = parlist.sublist("Problem").get("Boundary Cubature Degree",2); // set cubature degree, e.g., 2
    bdryCub_ = cubFactory.create(bdryCellType, bdryCubDegree);
    basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_); basisPtrs_.push_back(basisPtr_); basisPtrs_.push_back(basisPtr_);

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

    // Initialize residuals.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > R(3);
    R[0]   = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);
    R[1] = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);
    R[2] = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);

    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // COMPUTE STIFFNESS TERM
    ROL::Ptr<Intrepid::FieldContainer<Real> > U0_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(U0_eval, U[0]);
    // Integrate U * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*R[0],
                                                  *U0_eval,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    ROL::Ptr<Intrepid::FieldContainer<Real> > U1_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(U1_eval, U[1]);
    // Integrate U * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*R[1],
                                                  *U1_eval,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    ROL::Ptr<Intrepid::FieldContainer<Real> > U2_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(U2_eval, U[2]);
    // Integrate U * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*R[2],
                                                  *U2_eval,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);

  // Combine the residuals.
    fieldHelper_->combineFieldCoeff(res, R);
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
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > J(3);
    for (int i = 0; i < 3; ++i) {
      J[i].resize(3);
      for (int j = 0; j < 3; j++)
	{
	  J[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
	}
    }

    // Integrate N * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][0],
                                                  *(fe_->N()),
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);

    
    Intrepid::RealSpaceTools<Real>::scale(*J[1][1],0.0);
    Intrepid::RealSpaceTools<Real>::add(*J[1][1],*J[0][0]);
    Intrepid::RealSpaceTools<Real>::scale(*J[2][2],0.0);
    Intrepid::RealSpaceTools<Real>::add(*J[2][2],*J[0][0]);

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);
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
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > J(3);
    for (int i = 0; i < 3; ++i) {
      J[i].resize(3);
      for (int j = 0; j < 3; j++)
	{
	  J[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
	}
    }

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);
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

  void setFieldPattern(const std::vector<std::vector<int> > & fieldPattern) {
    fieldPattern_ = fieldPattern;
    fieldHelper_ = ROL::makePtr<FieldHelper<Real>>(numFields_, numDofs_, numFieldDofs_, fieldPattern_);
  }

 void setCellNodes(const ROL::Ptr<Intrepid::FieldContainer<Real> > &volCellNodes,
                    const std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int> > > &bdryCellLocIds) {
    volCellNodes_ = volCellNodes;
    bdryCellNodes_ = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    // Finite element definition.
    fe_ = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtr_,cellCub_);
    // Set local boundary DOFs.
    fidx_ = fe_->getBoundaryDofs();
    // Compute Dirichlet values at DOFs.
    int d = basisPtr_->getBaseCellTopology().getDimension();
    int numSidesets = bdryCellLocIds_.size();
    bdryCellDofValues_.resize(numSidesets);
    feBdry_.resize(numSidesets); 
    for (int i=0; i<numSidesets; ++i) {
      int numLocSides = bdryCellLocIds_[i].size();
      bdryCellDofValues_[i].resize(numLocSides);
      feBdry_[i].resize(numLocSides);
      for (int j=0; j<numLocSides; ++j) {
        int c = bdryCellLocIds_[i][j].size();
        int f = basisPtr_->getCardinality();
	if(c > 0)
	  {
	    feBdry_[i][j] = ROL::makePtr<FE<Real> >(bdryCellNodes_[i][j],basisPtr_,bdryCub_,j);
	  }
	bdryCellDofValues_[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
        ROL::Ptr<Intrepid::FieldContainer<Real> > coords =
          ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, d);
        if (c > 0) {
          fe_->computeDofCoords(coords, bdryCellNodes_[i][j]);
        }
        for (int k=0; k<c; ++k) {
          for (int l=0; l<f; ++l) {
            std::vector<Real> dofpoint(d);
            for (int m=0; m<d; ++m) {
              dofpoint[m] = (*coords)(k, l, m);
            }
            (*bdryCellDofValues_[i][j])(k, l) = 0.0;
          }
        }
      }
    }

  }

  ROL::Ptr<Intrepid::FieldContainer<Real> > getBoundaryCoeff(const Intrepid::FieldContainer<Real> & cell_coeff, int sideSet, int cell) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[sideSet][cell];
    const int numCellsSide = bdryCellLocId.size();
    const int f = basisPtr_->getCardinality();
    
    ROL::Ptr<Intrepid::FieldContainer<Real > > bdry_coeff = 
      ROL::makePtr<Intrepid::FieldContainer<Real > >(numCellsSide, f);
    for (int i = 0; i < numCellsSide; ++i) {
      for (int j = 0; j < f; ++j) {
	(*bdry_coeff)(i, j) = cell_coeff(bdryCellLocId[i], j);
      }
    }
    return bdry_coeff;
  }

  const ROL::Ptr<FE<Real> > getFE(void) const {
    return fe_;
  }

}; // PDE_Mass_Mat

// template<class Real>
// void Map_Index_to_Coords(int k, std::vector<Real> & coords, int d, int nx, int ny, Real nx_float, Real ny_float)
// {
//   coords[0] = static_cast<Real>( (k%(nx+1)) )*(1.0/nx_float);
//   coords[1] = static_cast<Real>( std::floor( static_cast<Real>(k)/(nx_float+1.0) ) )*(1.0/ny_float);
// }

// template<class Real>
// Real Mean_Eval(Real x, Real y)
// {
//   Real val = x*(1.0-x)*y*(1.0-y);
//   return val;
// }

template<class Real>
void Set_Prior_Mean(HDSA::Ptr<Tpetra::MultiVector<> > & z_ptr, const Teuchos::RCP<Teuchos::ParameterList> & parlist)
{
  int num_coeff_load = parlist->sublist("Problem").get("Number of Coefficients in Loaded Fields", 10); 
  std::vector<Real> prior_mean_coeff = std::vector<Real>(num_coeff_load);
  // read in data
  std::ifstream in("matlab_prior_mean.txt");          
  // read the elements in the file into a vector  
  // test file open   
  if (in) 
    {   
      for(int j = 0; j < num_coeff_load; j++)
	{
	  in >> prior_mean_coeff[j];
	}
    }
  else
    {
      std::cout << "Error loading the data from matlab_prior_mean.txt" << std::endl;
    }  
  
  for(int k = 0; k < num_coeff_load; k++)
    {
      z_ptr->replaceGlobalValue(3*k,0,prior_mean_coeff[k]);
    }

  // int nx = parlist->sublist("Geometry").get("NX",0);
  // int ny = parlist->sublist("Geometry").get("NY",0);
  // int dim = (nx+1)*(ny+1);
  // Real nx_float = static_cast<Real>(nx);
  // Real ny_float = static_cast<Real>(ny);
  // std::vector<Real> coords = std::vector<Real>(2,0.0);
  // for(int k = 0; k < dim; k++)
  //   {
  //     Map_Index_to_Coords<Real>(3*k,coords,dim,nx,ny,nx_float,ny_float);
  //     z_ptr->replaceGlobalValue(k,0,Mean_Eval<Real>(coords[0],coords[1]));
  //   }
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
	    zp->set(*Kz->basis(j));
	    con_->applyJacobian_1(*Kz,*zp,*tmp1,*tmp1,tol);
	    mass_mat_con_->applyJacobian_1(*Mz,*zp,*tmp1,*tmp1,tol);
	    for(int i = 0; i < n; i++)
	      {
		K[i][j] = Kz->dot(*zp->basis(i));
		M[i][j] = Mz->dot(*zp->basis(i));
	      }
	  }
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
