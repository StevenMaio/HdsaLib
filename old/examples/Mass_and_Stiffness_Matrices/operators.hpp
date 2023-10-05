#ifndef ELLIPTIC_OPERATORS_HPP
#define ELLIPTIC_OPERATORS_HPP

#include "../../../PDE-OPT/TOOLS/linearpdeconstraint.hpp"

template <class Real>
class PDE_Stiffness_Op : public PDE<Real> {
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
  PDE_Stiffness_Op(Teuchos::ParameterList &parlist) {
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
    basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_);

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
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > R(1);
    R[0]   = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);

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
    
    Intrepid::RealSpaceTools<Real>::scale(*R[0],0.0);
    Intrepid::RealSpaceTools<Real>::add(*R[0],*laplace_term);

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
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > J(1);
    for (int i = 0; i < 1; ++i) {
      J[i].resize(1);
      for (int j = 0; j < 1; j++)
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

    Intrepid::RealSpaceTools<Real>::scale(*J[0][0],0.0);
    Intrepid::RealSpaceTools<Real>::add(*J[0][0],*laplace_term);

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
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > J(1);
    for (int i = 0; i < 1; ++i) {
      J[i].resize(1);
      for (int j = 0; j < 1; j++)
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
    throw Exception::Zero(">>> (PDE_Stiffness_Op::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Stiffness_Op::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Stiffness_Op::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Stiffness_Op::Hessian_22): Hessian is zero.");
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

}; // PDE_Stiffness_Op

template <class Real>
class PDE_Mass_Op : public PDE<Real> {
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
  PDE_Mass_Op(Teuchos::ParameterList &parlist) {
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
    basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_);

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
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > R(1);
    R[0]   = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);
  
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
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > J(1);
    for (int i = 0; i < 1; ++i) {
      J[i].resize(1);
      for (int j = 0; j < 1; j++)
	{
	  J[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
	}
    }

    // Integrate N * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][0],
                                                  *(fe_->N()),
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);

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
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > J(1);
    for (int i = 0; i < 1; ++i) {
      J[i].resize(1);
      for (int j = 0; j < 1; j++)
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
    throw Exception::Zero(">>> (PDE_Stiffness_Op::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Stiffness_Op::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Stiffness_Op::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Stiffness_Op::Hessian_22): Hessian is zero.");
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


#endif
