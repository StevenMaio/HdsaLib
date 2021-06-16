#ifndef CONTROL_MASS_MAT_HPP
#define CONTROL_MASS_MAT_HPP

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
class Control_Mass_Mat : public PDE<Real> {
private:
  // Finite element basis information
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > basisPtrVel_;
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > basisPtrPrs_;
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > basisPtrThr_;
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
  ROL::Ptr<FE<Real> > feThr_;
  std::vector<std::vector<ROL::Ptr<FE<Real> > > > feVelBdry_;
  std::vector<std::vector<ROL::Ptr<FE<Real> > > > fePrsBdry_;
  std::vector<std::vector<ROL::Ptr<FE<Real> > > > feThrBdry_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int> > fvidx_;
  std::vector<std::vector<int> > fpidx_;
  std::vector<std::vector<int> > fhidx_;
  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellVDofValues_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellTDofValues_;
  // Field pattern, offsets, etc.
  std::vector<std::vector<int> > fieldPattern_;  // local Field/DOF pattern; set from DOF manager 
  int numFields_;                                // number of fields (equations in the PDE)
  int numDofs_;                                  // total number of degrees of freedom for all (local) fields
  std::vector<int> offset_;                      // for each field, a counting offset
  std::vector<int> numFieldDofs_;                // for each field, number of degrees of freedom

  ROL::Ptr<FieldHelper<Real> > fieldHelper_;

  ROL::Ptr<Intrepid::FieldContainer<Real> > getThermalBoundaryCoeff(
      const Intrepid::FieldContainer<Real> & cell_coeff,
      int sideSet, int cell) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[sideSet][cell];
    const int numCellsSide = bdryCellLocId.size();
    const int f = basisPtrThr_->getCardinality();
    
    ROL::Ptr<Intrepid::FieldContainer<Real > > bdry_coeff
      = ROL::makePtr<Intrepid::FieldContainer<Real >>(numCellsSide, f);
    for (int i = 0; i < numCellsSide; ++i) {
      for (int j = 0; j < f; ++j) {
        (*bdry_coeff)(i, j) = cell_coeff(bdryCellLocId[i], j);
      }
    }
    return bdry_coeff;
  }

public:
  Control_Mass_Mat(Teuchos::ParameterList &parlist) {
    // Finite element fields -- NOT DIMENSION INDEPENDENT!
    basisPtrVel_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real> >>();
    basisPtrPrs_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >>();
    basisPtrThr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real> >>();
    // Volume quadrature rules.
    shards::CellTopology cellType = basisPtrVel_->getBaseCellTopology();         // get the cell type from any basis
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
    for (int i = 0; i < d; ++i) {
      basisPtrs_.push_back(basisPtrVel_); // Velocity
    }
    basisPtrs_.push_back(basisPtrPrs_); // Pressure
    basisPtrs_.push_back(basisPtrThr_); // Heat

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
    // Retrieve dimensions.
    int c  = u_coeff->dimension(0);
    int p  = cellCub_->getNumPoints();
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int fh = basisPtrThr_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize residuals.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > R(d+2);
    for (int i = 0; i < d; ++i) {
      R[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
    }
    R[d]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    R[d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh);

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
    int fh = basisPtrThr_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize jacobians.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > J(d+2);
    for (int i = 0; i < d+2; ++i) {
      J[i].resize(d+2);
    }
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        J[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, fv);
      }
      J[d][i]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fv);
      J[i][d]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, fp);
      J[d+1][i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh, fv);
      J[i][d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, fh);
    }
    J[d][d]     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fp);
    J[d+1][d]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh, fp);
    J[d][d+1]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fh);
    J[d+1][d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh, fh);

    // APPLY BOUNDARY CONDITIONS
    int numSideSets = bdryCellLocIds_.size();
    const int numCubPerSide = bdryCub_->getNumPoints();
    if (numSideSets > 0) {
      // ROBIN CONDITIONS
      for (int i = 0; i < numSideSets; ++i) {
        // Control boundaries
        if (i==1 || i==2) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            ROL::Ptr<Intrepid::FieldContainer<Real> > robinJac;
            if (numCellsSide) {
              // Evaluate boundary integral
	      ROL::Ptr<Intrepid::FieldContainer<Real> > thr_bc = ROL::makePtr<Intrepid::FieldContainer<Real> >(numCellsSide, fh, fh);
              Intrepid::FunctionSpaceTools::integrate<Real>(*thr_bc,
                                                            *(feThrBdry_[i][j]->N()),
							    *(feThrBdry_[i][j]->NdetJ()),
							    Intrepid::COMP_CPP, false);
              for (int k = 0; k < numCellsSide; ++k) {
                int cidx = bdryCellLocIds_[i][j][k];
                for (int l = 0; l < fh; ++l) {
                  for (int m = 0; m < fh; ++m) {
                    (*J[d+1][d+1])(cidx,l,m) += (*thr_bc)(k,l,m);
                  }
                }
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
    // Retrieve dimensions.
    int c  = u_coeff->dimension(0);
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int fh = basisPtrThr_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize jacobians.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > J(d+2);
    for (int i = 0; i < d+2; ++i) {
      J[i].resize(d+2);
    }
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        J[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, fv);
      }
      J[d][i]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fv);
      J[i][d]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, fp);
      J[d+1][i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh, fv);
      J[i][d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, fh);
    }
    J[d][d]     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fp);
    J[d+1][d]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh, fp);
    J[d][d+1]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fh);
    J[d+1][d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh, fh);

    // APPLY DIRICHLET CONDITIONS
    int numSideSets = bdryCellLocIds_.size();
    const int numCubPerSide = bdryCub_->getNumPoints();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        // Control boundaries
        if (i==1 || i==2) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            if (numCellsSide) {
              // Evaluate boundary integral
              Intrepid::FieldContainer<Real> robinJac(numCellsSide, fh, fh);
              Intrepid::FunctionSpaceTools::integrate<Real>(robinJac,
                                                            *(feThrBdry_[i][j]->N()),
                                                            *(feThrBdry_[i][j]->NdetJ()),
                                                            Intrepid::COMP_CPP, false);
              for (int k = 0; k < numCellsSide; ++k) {
                int cidx = bdryCellLocIds_[i][j][k];
                for (int l = 0; l < fh; ++l) {
                  for (int m = 0; m < fh; ++m) {
                    (*J[d+1][d+1])(cidx,l,m) += robinJac(k,l,m);
                  }
                }
              }
            }
          }
        }
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
    throw Exception::Zero(">>> (PDE_Control_Mass_Mat::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Control_Mass_Mat::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Control_Mass_Mat::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Control_Mass_Mat::Hessian_22): Hessian is zero.");
  }

  void RieszMap_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & riesz) {
    // Retrieve dimensions.
    int c  = feVel_->N()->dimension(0);
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int fh = basisPtrThr_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize jacobians.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > J(d+2);
    for (int i = 0; i < d+2; ++i) {
      J[i].resize(d+2);
    }
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        J[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, fv);
      }
      J[d][i]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fv);
      J[i][d]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, fp);
      J[d+1][i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh, fv);
      J[i][d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, fh);
    }
    J[d][d]     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fp);
    J[d+1][d]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh, fp);
    J[d][d+1]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fh);
    J[d+1][d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh, fh);

    for (int i = 0; i < d; ++i) {
      *(J[i][i]) = *(feVel_->stiffMat());
      Intrepid::RealSpaceTools<Real>::add(*(J[i][i]),*(feVel_->massMat()));
    }
    *(J[d][d]) = *(fePrs_->massMat());
    *(J[d+1][d+1]) = *(feThr_->stiffMat());
    Intrepid::RealSpaceTools<Real>::add(*(J[d+1][d+1]),*(feThr_->massMat()));

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(riesz, J);
  }

  void RieszMap_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & riesz) {
    // Retrieve dimensions.
    int c  = feVel_->N()->dimension(0);
    int fv = basisPtrVel_->getCardinality();
    int fp = basisPtrPrs_->getCardinality();
    int fh = basisPtrThr_->getCardinality();
    int d  = cellCub_->getDimension();
 
    // Initialize jacobians.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > J(d+2);
    for (int i = 0; i < d+2; ++i) {
      J[i].resize(d+2);
    }
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        J[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, fv);
      }
      J[d][i]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fv);
      J[i][d]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, fp);
      J[d+1][i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh, fv);
      J[i][d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, fh);
    }
    J[d][d]     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fp);
    J[d+1][d]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh, fp);
    J[d][d+1]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp, fh);
    J[d+1][d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh, fh);


    for (int i = 0; i < d; ++i) {
      *(J[i][i]) = *(feVel_->massMat());
    }
    *(J[d][d]) = *(fePrs_->massMat());
    *(J[d+1][d+1]) = *(feThr_->massMat());

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
    feThr_ = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtrThr_,cellCub_);
    // Get boundary degrees of freedom.
    fvidx_ = feVel_->getBoundaryDofs();
    fpidx_ = fePrs_->getBoundaryDofs();
    fhidx_ = feThr_->getBoundaryDofs();
    // Construct boundary FEs
    const int numSideSets = bdryCellNodes.size();
    feVelBdry_.resize(numSideSets);
    fePrsBdry_.resize(numSideSets);
    feThrBdry_.resize(numSideSets);
    for (int i = 0; i < numSideSets; ++i) {
      int numLocSides = bdryCellNodes[i].size();
      feVelBdry_[i].resize(numLocSides);
      fePrsBdry_[i].resize(numLocSides);
      feThrBdry_[i].resize(numLocSides);
      for (int j = 0; j < numLocSides; ++j) {
        if (bdryCellNodes[i][j] != ROL::nullPtr) {
          feVelBdry_[i][j] = ROL::makePtr<FE<Real>>(bdryCellNodes[i][j],basisPtrVel_,bdryCub_,j);
          fePrsBdry_[i][j] = ROL::makePtr<FE<Real>>(bdryCellNodes[i][j],basisPtrPrs_,bdryCub_,j);
          feThrBdry_[i][j] = ROL::makePtr<FE<Real>>(bdryCellNodes[i][j],basisPtrThr_,bdryCub_,j);
        }
      }
    }
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

  const ROL::Ptr<FE<Real> > getThermalFE(void) const {
    return feThr_;
  }

  const std::vector<std::vector<ROL::Ptr<FE<Real> > > > getVelocityBdryFE(void) const {
    return feVelBdry_;
  }

  const std::vector<std::vector<ROL::Ptr<FE<Real> > > > getPressureBdryFE(void) const {
    return fePrsBdry_;
  }

  const std::vector<std::vector<ROL::Ptr<FE<Real> > > > getThermalBdryFE(void) const {
    return feThrBdry_;
  }

  const std::vector<std::vector<std::vector<int> > > getBdryCellLocIds(void) const {
    return bdryCellLocIds_;
  }

  const ROL::Ptr<FieldHelper<Real> > getFieldHelper(void) const {
    return fieldHelper_;
  }

}; // PDE_Control_Mass_Mat

#endif
