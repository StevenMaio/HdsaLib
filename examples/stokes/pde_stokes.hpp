#ifndef PDE_STOKES_HPP
#define PDE_STOKES_HPP

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
class PDE_Stokes : public PDE<Real> {
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
  std::vector<ROL::Ptr<FE<Real> > > feVelBdry_;
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
  Real viscosity_;
  Real gravity_;
  Real advection_coeff_;
  Real channelW_;
  Real stepW_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > ctrlWeight_;  
  Real initiate_control_;

  ROL::Ptr<FieldHelper<Real> > fieldHelper_;

  Real DirichletFunc(const std::vector<Real> & coords, int sideset, int locSideId, int dir) const {
    Real val(0);
    if ((sideset==4) && (dir==0)) { // Left upper x-velocity
        val = 8.0*(coords[1]-0.5)*(1.0-coords[1]);
    }
    if ((sideset==1) && (dir==0)) { // Right lower x-velocity
      val = 12.0*coords[1]*(1.0-coords[1]);
    }
    if ((sideset==2) && (dir==0)) { // Right upper x-velocity
      val = 12.0*coords[1]*(1.0-coords[1]);
    }
   if ((sideset==3) && (dir==0)) { // Top x-velocity
     val = coords[0]/channelW_;
    }
   if ((sideset==0) && (dir==0)) { // Bottom x-velocity
     val = coords[0]/channelW_;
    }
    if ((sideset==4) && (dir==1)) { // Left upper y-velocity
        val = 0.0;
    }
   if ((sideset==5) && (dir==1)) { // Middle y-velocity
      val = 0.0;
    }
    if ((sideset==1) && (dir==1)) { // Right lower y-velocity
      val = -1.5 + coords[1];
    }
    if ((sideset==2) && (dir==1)) { // Right upper y-velocity
      val = -1.5 + coords[1];
    }
   if ((sideset==3) && (dir==1)) { // Top y-velocity
      val = -0.5*coords[0]/channelW_;
    }
   if ((sideset==0) && (dir==1)) { // Bottom y-velocity
      val = -1.5*coords[0]/channelW_;
    }
   if ((sideset==4) && (dir==2)) { // Left upper pressure
      val = 10.0;
    }
   if ((sideset==5) && (dir==2)) { // middle pressure
      val = 10.0;
    }
   if ((sideset==6) && (dir==2)) { // Left lower pressure
      val = 10.0;
    }
   if ((sideset==1) && (dir==2)) { // Right lower pressure
      val = 2.0;
    }
   if ((sideset==2) && (dir==2)) { // Right upper pressure
      val = 2.0;
    }
    return val;
  }

  void computeDirichlet(void) {
    // Compute Dirichlet values at DOFs.
    int d = basisPtrVel_->getBaseCellTopology().getDimension();
    int numSidesets = bdryCellLocIds_.size();
    bdryCellDofValues_.resize(numSidesets);
    for (int i=0; i<numSidesets; ++i) {
      int numLocSides = bdryCellLocIds_[i].size();
      bdryCellDofValues_[i].resize(numLocSides);
      for (int j=0; j<numLocSides; ++j) {
        int c = bdryCellLocIds_[i][j].size();
        int f = basisPtrVel_->getCardinality();
        bdryCellDofValues_[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, d+1);
        ROL::Ptr<Intrepid::FieldContainer<Real> > coords =
          ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, d);
        if (c > 0) {
          feVel_->computeDofCoords(coords, bdryCellNodes_[i][j]);
        }
        for (int k=0; k<c; ++k) {
          for (int l=0; l<f; ++l) {
            std::vector<Real> dofpoint(d);
            //std::cout << "Sideset " << i << " LocalSide " << j << "  Cell " << k << "  Field " << l << "  Coord ";
            for (int m=0; m<d; ++m) {
              dofpoint[m] = (*coords)(k, l, m);
              //std::cout << dofpoint[m] << "  ";
            }

            for (int m=0; m<d+1; ++m) {
              (*bdryCellDofValues_[i][j])(k, l, m) = DirichletFunc(dofpoint, i, j, m);
              //std::cout << "  " << m << "-Value " << DirichletFunc(dofpoint, i, j, m);
            }
            //std::cout << std::endl;
          }
        }
      }
    }
  }

  // void computeDirichlet(void) {
  //   // Compute Dirichlet values at DOFs.
  //   int d = basisPtrVel_->getBaseCellTopology().getDimension();
  //   int numSidesets = bdryCellLocIds_.size();
  //   bdryCellDofValues_.resize(numSidesets);
  //   for (int i=0; i<numSidesets; ++i) {
  //     int numLocSides = bdryCellLocIds_[i].size();
  //     bdryCellDofValues_[i].resize(numLocSides);
  //     for (int j=0; j<numLocSides; ++j) {
  //       int c = bdryCellLocIds_[i][j].size();
  //       int f = basisPtrVel_->getCardinality();
  //       bdryCellDofValues_[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, d);
  //       ROL::Ptr<Intrepid::FieldContainer<Real> > coords =
  //         ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, d);
  //       if (c > 0) {
  //         feVel_->computeDofCoords(coords, bdryCellNodes_[i][j]);
  //       }
  //       for (int k=0; k<c; ++k) {
  //         for (int l=0; l<f; ++l) {
  //           std::vector<Real> dofpoint(d);
  //           //std::cout << "Sideset " << i << " LocalSide " << j << "  Cell " << k << "  Field " << l << "  Coord ";
  //           for (int m=0; m<d; ++m) {
  //             dofpoint[m] = (*coords)(k, l, m);
  //             //std::cout << dofpoint[m] << "  ";
  //           }

  //           for (int m=0; m<d; ++m) {
  //             (*bdryCellDofValues_[i][j])(k, l, m) = DirichletFunc(dofpoint, i, j, m);
  //             //std::cout << "  " << m << "-Value " << DirichletFunc(dofpoint, i, j, m);
  //           }
  //           //std::cout << std::endl;
  //         }
  //       }
  //     }
  //   }
  // }

  Real viscosityFunc(const std::vector<Real> & coords) const {
    return viscosity_;
  }

  void computeViscosity(ROL::Ptr<Intrepid::FieldContainer<Real> > &nu) const {
    int c = feVel_->gradN()->dimension(0);
    int p = feVel_->gradN()->dimension(2);
    int d = feVel_->gradN()->dimension(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k) {
          pt[k] = (*feVel_->cubPts())(i,j,k);
        }
        // Compute spatially dependent viscosity
        (*nu)(i,j) = viscosityFunc(pt);
      }
    }
  }

  void computeControlWeight(void) {
    int c = feVel_->cubPts()->dimension(0);
    int p = feVel_->cubPts()->dimension(1);
    int d = feVel_->cubPts()->dimension(2);

    ctrlWeight_ = ROL::makePtr<Intrepid::FieldContainer<Real> >(c,p);
 
    const Real zero(0), one(1);
    bool inside(false);
    std::vector<Real> x(d);
    for (int i = 0; i < c; ++i) {
      inside = false;
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          x[k] = (*feVel_->cubPts())(i,j,k);
        }
        if ( insideDomain(x) ) {
          inside = true;
          break;
        }
      }
      for (int j = 0; j < p; ++j) {
        (*ctrlWeight_)(i,j) = (inside ? one : zero);
      }
    }
  }

  bool insideDomain(const std::vector<Real> &x) const {
    bool inside(false);
    if(x[0] > initiate_control_)
      {
	inside = true;
      }
    else
      {
	inside = false;
      }
    return inside;
  }

public:
  PDE_Stokes(Teuchos::ParameterList &parlist, Real advection_coeff) {
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
    viscosity_ = parlist.sublist("Problem").get("Viscosity", 5e-3);
    gravity_ = parlist.sublist("Problem").get("Gravity", 9.81);
    advection_coeff_ = advection_coeff;
    channelW_ = parlist.sublist("Geometry").get(    "Channel width", 8.0);
    stepW_    = parlist.sublist("Geometry").get(       "Step width", 1.0);
    initiate_control_ = parlist.sublist("Geometry").get("Initiate Control", -1.0);

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
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    fieldHelper_->splitFieldCoeff(Z, z_coeff);

    // Evaluate/interpolate finite element fields on cells.
    ROL::Ptr<Intrepid::FieldContainer<Real> > nu =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > nuGradVelX_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    ROL::Ptr<Intrepid::FieldContainer<Real> > nuGradVelY_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valVel_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valVelX_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valVelY_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valPres_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradVelX_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradVelY_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    ROL::Ptr<Intrepid::FieldContainer<Real> > divVel_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valVelDotgradVelX_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valVelDotgradVelY_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    feVel_->evaluateValue(valVelX_eval, U[0]);
    feVel_->evaluateValue(valVelY_eval, U[1]);
    fePrs_->evaluateValue(valPres_eval, U[2]);
    feVel_->evaluateGradient(gradVelX_eval, U[0]);
    feVel_->evaluateGradient(gradVelY_eval, U[1]);
    computeViscosity(nu);

    // Assemble the velocity vector and its divergence.
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*valVel_eval)(i,j,0) = (*valVelX_eval)(i,j);
        (*valVel_eval)(i,j,1) = (*valVelY_eval)(i,j);
        (*divVel_eval)(i,j)   = (*gradVelX_eval)(i,j,0) + (*gradVelY_eval)(i,j,1);
      }
    }
    // Negative pressure
    Intrepid::RealSpaceTools<Real>::scale(*valPres_eval,static_cast<Real>(-1));

    // Multiply velocity gradients with viscosity.
    Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(*nuGradVelX_eval,
                                                               *nu,
                                                               *gradVelX_eval);
    Intrepid::FunctionSpaceTools::tensorMultiplyDataData<Real>(*nuGradVelY_eval,
                                                               *nu,
                                                               *gradVelY_eval);

    // Compute nonlinear terms in the Navier-Stokes equations.
    Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(*valVelDotgradVelX_eval, *valVel_eval, *gradVelX_eval);
    Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(*valVelDotgradVelY_eval, *valVel_eval, *gradVelY_eval);

    Intrepid::RealSpaceTools<Real>::scale(*valVelDotgradVelX_eval,advection_coeff_);
    Intrepid::RealSpaceTools<Real>::scale(*valVelDotgradVelY_eval,advection_coeff_);

    /*** Evaluate weak form of the residual. ***/
    // X component of velocity equation.
    Intrepid::FunctionSpaceTools::integrate<Real>(velX_res,
                                                  *nuGradVelX_eval,         // nu gradUX
                                                  *(feVel_->gradNdetJ()),   // gradPhi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>(velX_res,
                                                  *valVelDotgradVelX_eval,  // (U . gradUX)
                                                  *(feVel_->NdetJ()),       // Phi
                                                  Intrepid::COMP_CPP,
                                                  true);
    Intrepid::FunctionSpaceTools::integrate<Real>(velX_res,
                                                  *valPres_eval,            // p
                                                  *(feVel_->DNDdetJ(0)),    // dPhi/dx
                                                  Intrepid::COMP_CPP,
                                                  true);
    
    ROL::Ptr<Intrepid::FieldContainer<Real> > valZX_eval = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    feVel_->evaluateValue(valZX_eval, Z[0]);
    Intrepid::RealSpaceTools<Real>::scale(*valZX_eval,static_cast<Real>(-1));
    ROL::Ptr<Intrepid::FieldContainer<Real> > weight_valZX_eval = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*weight_valZX_eval, *ctrlWeight_, *valZX_eval);
    Intrepid::FunctionSpaceTools::integrate<Real>(velX_res, *weight_valZX_eval, *(feVel_->NdetJ()), Intrepid::COMP_CPP, true);

    // Y component of velocity equation.
    Intrepid::FunctionSpaceTools::integrate<Real>(velY_res,
                                                  *nuGradVelY_eval,         // nu gradUY
                                                  *(feVel_->gradNdetJ()),   // gradPhi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>(velY_res,
                                                  *valVelDotgradVelY_eval,  // (U . gradUY)
                                                  *(feVel_->NdetJ()),       // Phi
                                                  Intrepid::COMP_CPP,
                                                  true);
    Intrepid::FunctionSpaceTools::integrate<Real>(velY_res,
                                                  *valPres_eval,            // p
                                                  *(feVel_->DNDdetJ(1)),    // dPhi/dy
                                                  Intrepid::COMP_CPP,
                                                  true);

    ROL::Ptr<Intrepid::FieldContainer<Real> > valGrav_eval = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    for(int i = 0; i < c; i++)
      {
	for(int j = 0; j < p; j++)
	  {
	    (*valGrav_eval)(i,j) = gravity_;
	  }
      }
    Intrepid::FunctionSpaceTools::integrate<Real>(velY_res, *valGrav_eval, *(feVel_->NdetJ()), Intrepid::COMP_CPP, true);

    ROL::Ptr<Intrepid::FieldContainer<Real> > valZY_eval = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    feVel_->evaluateValue(valZY_eval, Z[1]);
    Intrepid::RealSpaceTools<Real>::scale(*valZY_eval,static_cast<Real>(-1));
    ROL::Ptr<Intrepid::FieldContainer<Real> > weight_valZY_eval = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*weight_valZY_eval, *ctrlWeight_, *valZY_eval);
    Intrepid::FunctionSpaceTools::integrate<Real>(velY_res, *weight_valZY_eval, *(feVel_->NdetJ()), Intrepid::COMP_CPP, true);

    // Pressure equation.
    Intrepid::FunctionSpaceTools::integrate<Real>(pres_res,
                                                  *divVel_eval,             // divU
                                                  *(fePrs_->NdetJ()),       // Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::RealSpaceTools<Real>::scale(pres_res,static_cast<Real>(-1));

    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      // APPLY DIRICHLET CONDITIONS
      computeDirichlet();
      for (int i = 0; i < numSideSets; ++i) {
        // Apply Dirichlet conditions
        if (i==4) {
	    int numLocalSideIds = bdryCellLocIds_[i].size();
	    for (int j = 0; j < numLocalSideIds; ++j) {
	      int numCellsSide = bdryCellLocIds_[i][j].size();
	      int numBdryDofs = fvidx_[j].size();
	      for (int k = 0; k < numCellsSide; ++k) {
		int cidx = bdryCellLocIds_[i][j][k];
		for (int l = 0; l < numBdryDofs; ++l) {
		  (*R[0])(cidx,fvidx_[j][l]) = (*U[0])(cidx,fvidx_[j][l]) - (*bdryCellDofValues_[i][j])(k,fvidx_[j][l],0);
		}
	      }
	    }
	}
        if ((i==4) || (i==5)) {
	    int numLocalSideIds = bdryCellLocIds_[i].size();
	    for (int j = 0; j < numLocalSideIds; ++j) {
	      int numCellsSide = bdryCellLocIds_[i][j].size();
	      int numBdryDofs = fvidx_[j].size();
	      for (int k = 0; k < numCellsSide; ++k) {
		int cidx = bdryCellLocIds_[i][j][k];
		for (int l = 0; l < numBdryDofs; ++l) {
		  (*R[1])(cidx,fvidx_[j][l]) = (*U[1])(cidx,fvidx_[j][l]) - (*bdryCellDofValues_[i][j])(k,fvidx_[j][l],1);
		}
	      }
	    }
	}
        if ( (i==4) ) {
	  int numLocalSideIds = bdryCellLocIds_[i].size();
	  for (int j = 0; j < numLocalSideIds; ++j) {
	    int numCellsSide = bdryCellLocIds_[i][j].size();
	    int numBdryDofs = fpidx_[j].size();
	    for (int k = 0; k < numCellsSide; ++k) {
	      int cidx = bdryCellLocIds_[i][j][k];
	      for (int l = 0; l < numBdryDofs; ++l) {
		(*R[d])(cidx,fpidx_[j][l]) = (*U[d])(cidx,fpidx_[j][l]) - (*bdryCellDofValues_[i][j])(k,fpidx_[j][l],d);
	      }
	    }
	  }
	}
      }
    }

    // int numSideSets = bdryCellLocIds_.size();
    // if (numSideSets > 0) {
    //   // APPLY DIRICHLET CONDITIONS
    //   computeDirichlet();
    //   for (int i = 0; i < numSideSets; ++i) {
    //     // Apply Dirichlet conditions
    //     if (i!=9) {
    // 	    int numLocalSideIds = bdryCellLocIds_[i].size();
    // 	    for (int j = 0; j < numLocalSideIds; ++j) {
    // 	      int numCellsSide = bdryCellLocIds_[i][j].size();
    // 	      int numBdryDofs = fvidx_[j].size();
    // 	      for (int k = 0; k < numCellsSide; ++k) {
    // 		int cidx = bdryCellLocIds_[i][j][k];
    // 		for (int l = 0; l < numBdryDofs; ++l) {
    // 		  //std::cout << "\n  i=" << i << "   cidx=" << cidx << "   j=" << j << "  l=" << l << "  " << fvidx_[j][l] << " " << (*bdryCellDofValues_[i][j])(k,fvidx_[j][l],0);
    // 		  //std::cout << "\n  i=" << i << "   cidx=" << cidx << "   j=" << j << "  l=" << l << "  " << fvidx_[j][l] << " " << (*bdryCellDofValues_[i][j])(k,fvidx_[j][l],1);
    // 		  for (int m=0; m < d; ++m) {
    // 		    (*R[m])(cidx,fvidx_[j][l]) = (*U[m])(cidx,fvidx_[j][l]) - (*bdryCellDofValues_[i][j])(k,fvidx_[j][l],m);
    // 		  }
    // 		}
    // 	      }
    // 	    }
    // 	  }
       
    //     // Pressure pinning
    //     if (i==9) {
    //       int j = 2, l = 0;
    //       if (bdryCellLocIds_[i][0].size() > 0) { 
    //         int cidx = bdryCellLocIds_[i][0][0];
    //         (*R[d])(cidx,fpidx_[j][l]) = (*U[d])(cidx,fpidx_[j][l]);
    //       }
    //     }
    //   }
    // }

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

    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);

    // Evaluate/interpolate finite element fields on cells.
    ROL::Ptr<Intrepid::FieldContainer<Real> > nu =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > nuGradPhiX_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, p, d);
    ROL::Ptr<Intrepid::FieldContainer<Real> > nuGradPhiY_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, p, d);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valVel_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valVelX_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valVelY_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradVelX_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    ROL::Ptr<Intrepid::FieldContainer<Real> > ddxVelX_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > ddyVelX_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > ddxVelXPhiX_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > ddyVelXPhiY_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradVelY_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    ROL::Ptr<Intrepid::FieldContainer<Real> > ddxVelY_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > ddyVelY_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > ddyVelYPhiY_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > ddxVelYPhiX_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valVelDotgradPhiX_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valVelDotgradPhiY_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, p);
    feVel_->evaluateValue(valVelX_eval, U[0]);
    feVel_->evaluateValue(valVelY_eval, U[1]);
    feVel_->evaluateGradient(gradVelX_eval, U[0]);
    feVel_->evaluateGradient(gradVelY_eval, U[1]);
    computeViscosity(nu);

    // Assemble the velocity vector and its divergence.
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*valVel_eval)(i,j,0) = (*valVelX_eval)(i,j);
        (*valVel_eval)(i,j,1) = (*valVelY_eval)(i,j);
        (*ddxVelX_eval)(i,j)  = (*gradVelX_eval)(i,j,0);
        (*ddyVelX_eval)(i,j)  = (*gradVelX_eval)(i,j,1);
        (*ddxVelY_eval)(i,j)  = (*gradVelY_eval)(i,j,0);
        (*ddyVelY_eval)(i,j)  = (*gradVelY_eval)(i,j,1);
      }
    }

    // Multiply velocity gradients with viscosity.
    Intrepid::FunctionSpaceTools::tensorMultiplyDataField<Real>(*nuGradPhiX_eval,
                                                                *nu,
                                                                *(feVel_->gradN()));
    Intrepid::FunctionSpaceTools::tensorMultiplyDataField<Real>(*nuGradPhiY_eval,
                                                                *nu,
                                                                *(feVel_->gradN()));

    // Compute nonlinear terms in the Navier-Stokes equations.
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(*valVelDotgradPhiX_eval, *valVel_eval, *(feVel_->gradN()));
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(*valVelDotgradPhiY_eval, *valVel_eval, *(feVel_->gradN()));
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*ddxVelXPhiX_eval, *ddxVelX_eval, *(feVel_->N()));
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*ddyVelXPhiY_eval, *ddyVelX_eval, *(feVel_->N()));
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*ddxVelYPhiX_eval, *ddxVelY_eval, *(feVel_->N()));
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*ddyVelYPhiY_eval, *ddyVelY_eval, *(feVel_->N()));

    Intrepid::RealSpaceTools<Real>::scale(*valVelDotgradPhiX_eval,advection_coeff_);
    Intrepid::RealSpaceTools<Real>::scale(*valVelDotgradPhiY_eval,advection_coeff_);
    Intrepid::RealSpaceTools<Real>::scale(*ddxVelXPhiX_eval,advection_coeff_);
    Intrepid::RealSpaceTools<Real>::scale(*ddyVelXPhiY_eval,advection_coeff_);
    Intrepid::RealSpaceTools<Real>::scale(*ddxVelYPhiX_eval,advection_coeff_);
    Intrepid::RealSpaceTools<Real>::scale(*ddyVelYPhiY_eval,advection_coeff_);

    // Negative pressure basis.
    Intrepid::FieldContainer<Real> negPrsN(*(fePrs_->N()));
    Intrepid::RealSpaceTools<Real>::scale(negPrsN,static_cast<Real>(-1));

    /*** Evaluate weak form of the Jacobian. ***/
    // X component of velocity equation.
    Intrepid::FunctionSpaceTools::integrate<Real>(velXvelX_jac,
                                                  *nuGradPhiX_eval,         // nu gradPsi
                                                  *(feVel_->gradNdetJ()),   // gradPhi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>(velXvelX_jac,
                                                  *(feVel_->NdetJ()),       // Phi
                                                  *valVelDotgradPhiX_eval,  // (U . gradPhiX)
                                                  Intrepid::COMP_CPP,
                                                  true);
    Intrepid::FunctionSpaceTools::integrate<Real>(velXvelX_jac,
                                                  *(feVel_->NdetJ()),       // Phi
                                                  *ddxVelXPhiX_eval,        // (Phi . gradU)
                                                  Intrepid::COMP_CPP,
                                                  true);
    Intrepid::FunctionSpaceTools::integrate<Real>(velXvelY_jac,
                                                  *(feVel_->NdetJ()),       // Phi
                                                  *ddyVelXPhiY_eval,        // (Phi . gradU)
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>(velXpres_jac,
                                                  *(feVel_->DNDdetJ(0)),    // dPhi/dx
                                                  negPrsN,                  // -Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    // Y component of velocity equation.
    Intrepid::FunctionSpaceTools::integrate<Real>(velYvelY_jac,
                                                  *nuGradPhiY_eval,         // nu gradPsi
                                                  *(feVel_->gradNdetJ()),   // gradPhi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>(velYvelY_jac,
                                                  *(feVel_->NdetJ()),       // Phi
                                                  *valVelDotgradPhiY_eval,  // (U . gradPhiX)
                                                  Intrepid::COMP_CPP,
                                                  true);
    Intrepid::FunctionSpaceTools::integrate<Real>(velYvelY_jac,
                                                  *(feVel_->NdetJ()),       // Phi
                                                  *ddyVelYPhiY_eval,        // (Phi . gradU)
                                                  Intrepid::COMP_CPP,
                                                  true);
    Intrepid::FunctionSpaceTools::integrate<Real>(velYvelX_jac,
                                                  *(feVel_->NdetJ()),       // Phi
                                                  *ddxVelYPhiX_eval,        // (Phi . gradU)
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>(velYpres_jac,
                                                  *(feVel_->DNDdetJ(1)),    // dPhi/dx
                                                  negPrsN,                  // -Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    // Pressure equation.
    Intrepid::FunctionSpaceTools::integrate<Real>(presvelX_jac,
                                                  *(fePrs_->NdetJ()),       // Phi
                                                  *(feVel_->DND(0)),        // dPhi/dx
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>(presvelY_jac,
                                                  *(fePrs_->NdetJ()),       // Phi
                                                  *(feVel_->DND(1)),        // dPhi/dx
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::RealSpaceTools<Real>::scale(presvelX_jac,static_cast<Real>(-1));
    Intrepid::RealSpaceTools<Real>::scale(presvelY_jac,static_cast<Real>(-1));

    // APPLY DIRICHLET CONDITIONS
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if (i==4) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m=0; m < fv; ++m) {
		  for (int p=0; p < d; ++p) {
		    (*J[0][p])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
		  }
		  (*J[0][0])(cidx,fvidx_[j][l],fvidx_[j][l]) = static_cast<Real>(1);
		  for (int m=0; m < fp; ++m) {
		    (*J[0][2])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
		  }
		}
              }
            }
	  }
	}
        if ((i==4) || (i==5)) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m=0; m < fv; ++m) {
		  for (int p=0; p < d; ++p) {
		    (*J[1][p])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
		  }
		  (*J[1][1])(cidx,fvidx_[j][l],fvidx_[j][l]) = static_cast<Real>(1);
		  for (int m=0; m < fp; ++m) {
		    (*J[1][2])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
		  }
		}
              }
            }
	  }
	}
        if ( (i==4) ) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
	    int numBdryDofs = fpidx_[j].size();
	    for (int k = 0; k < numCellsSide; ++k) {
	      int cidx = bdryCellLocIds_[i][j][k];
	      for (int l = 0; l < numBdryDofs; ++l) {
		for (int m = 0; m < fv; ++m) {
		  for (int n = 0; n < d; ++n) {
		    (*J[d][n])(cidx,fpidx_[j][l],m) = static_cast<Real>(0);
		  }
		}
		for (int m = 0; m < fp; ++m) {
		  (*J[d][d])(cidx,fpidx_[j][l],m) = static_cast<Real>(0);
		}
		(*J[d][d])(cidx,fpidx_[j][l],fpidx_[j][l]) = static_cast<Real>(1);
	      }
	    }
	  }
	}
      }
    }
    // int numSideSets = bdryCellLocIds_.size();
    // if (numSideSets > 0) {
    //   for (int i = 0; i < numSideSets; ++i) {
    //     if (i!=9) {
    //       int numLocalSideIds = bdryCellLocIds_[i].size();
    //       for (int j = 0; j < numLocalSideIds; ++j) {
    //         int numCellsSide = bdryCellLocIds_[i][j].size();
    //         int numBdryDofs = fvidx_[j].size();
    //         for (int k = 0; k < numCellsSide; ++k) {
    //           int cidx = bdryCellLocIds_[i][j][k];
    //           for (int l = 0; l < numBdryDofs; ++l) {
    //             for (int m=0; m < fv; ++m) {
    //               for (int n=0; n < d; ++n) {
    //                 for (int p=0; p < d; ++p) {
    //                   (*J[n][p])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
    //                 }
    //                 (*J[n][n])(cidx,fvidx_[j][l],fvidx_[j][l]) = static_cast<Real>(1);
    //               }
    //             }
    //             for (int m=0; m < fp; ++m) {
    //               for (int n=0; n < d; ++n) {
    //                 (*J[n][2])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
    //               }
    //             }
    //           }
    //         }
    //       }
    //     }
    //     // Pressure pinning
    //     if (i==9) {
    //       int j = 2, l = 0;
    //       if (bdryCellLocIds_[i][0].size() > 0) { 
    //         int cidx = bdryCellLocIds_[i][0][0];
    //         for (int m = 0; m < fv; ++m) {
    //           for (int n = 0; n < d; ++n) {
    //             (*J[d][n])(cidx,fpidx_[j][l],m) = static_cast<Real>(0);
    //           }
    //         }
    //         for (int m = 0; m < fp; ++m) {
    //           (*J[d][d])(cidx,fpidx_[j][l],m) = static_cast<Real>(0);
    //         }
    //         (*J[d][d])(cidx,fpidx_[j][l],fpidx_[j][l]) = static_cast<Real>(1);
    //       }
    //     }
    //   }
    // }

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

    // ADD CONTROL TERM    
    ROL::Ptr<Intrepid::FieldContainer<Real> > weight_Val_eval = ROL::makePtr<Intrepid::FieldContainer<Real> >(c,fv, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*weight_Val_eval, *ctrlWeight_, *(feVel_->N()));
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][0], *weight_Val_eval, *(feVel_->NdetJ()), Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*J[0][0],static_cast<Real>(-1));
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[1][1], *weight_Val_eval, *(feVel_->NdetJ()), Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*J[1][1],static_cast<Real>(-1));

    // APPLY DIRICHLET CONDITIONS
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        // Apply Dirichlet conditions
        if (i==4) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m=0; m < fv; ++m) {
		  (*J[0][0])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                }
              }
            }
          }
        }
        if ((i==4) || (i==5)) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fvidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m=0; m < fv; ++m) {
		  (*J[1][1])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
                }
              }
            }
          }
        }
        if ( (i==4) ) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
	    int numBdryDofs = fpidx_[j].size();
	    for (int k = 0; k < numCellsSide; ++k) {
	      int cidx = bdryCellLocIds_[i][j][k];
	      for (int l = 0; l < numBdryDofs; ++l) {
		for (int m=0; m < fp; m++){
		  (*J[d][d])(cidx,fpidx_[j][l],m) = static_cast<Real>(0);
		}
	      }
	    }
	  }
	}
      }
    }
    // int numSideSets = bdryCellLocIds_.size();
    // if (numSideSets > 0) {
    //   for (int i = 0; i < numSideSets; ++i) {
    //     // Apply Dirichlet conditions
    //     if (i!=9) {
    //       int numLocalSideIds = bdryCellLocIds_[i].size();
    //       for (int j = 0; j < numLocalSideIds; ++j) {
    //         int numCellsSide = bdryCellLocIds_[i][j].size();
    //         int numBdryDofs = fvidx_[j].size();
    //         for (int k = 0; k < numCellsSide; ++k) {
    //           int cidx = bdryCellLocIds_[i][j][k];
    //           for (int l = 0; l < numBdryDofs; ++l) {
    //             //std::cout << "\n   j=" << j << "  l=" << l << "  " << fidx[j][l];
    //             for (int m=0; m < fv; ++m) {
    //               for (int n=0; n < d; ++n) {
    //                 (*J[n][n])(cidx,fvidx_[j][l],m) = static_cast<Real>(0);
    //               }
    //             }
    //           }
    //         }
    //       }
    //     }
    //     // Pressure pinning
    //     if (i==9) {
    //       int j = 2, l = 0;
    //       if (bdryCellLocIds_[i][0].size() > 0) { 
    //         int cidx = bdryCellLocIds_[i][0][0];
    //         for (int m = 0; m < fv; ++m) {
    //           for (int n = 0; n < d; ++n) {
    //             (*J[d][n])(cidx,fpidx_[j][l],m) = static_cast<Real>(0);
    //           }
    //         }
    //         for (int m = 0; m < fp; ++m) {
    //           (*J[d][d])(cidx,fpidx_[j][l],m) = static_cast<Real>(0);
    //         }
    //       }
    // 	}
    //   }
    // }

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);
  }

  void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    //    throw Exception::NotImplemented("");
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

    // Split u_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > L;
    fieldHelper_->splitFieldCoeff(L, l_coeff);

    // Apply Dirichlet conditions
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
	if (i!=9) {
	  int numLocalSideIds = bdryCellLocIds_[i].size();
	  for (int j = 0; j < numLocalSideIds; ++j) {
	    int numCellsSide = bdryCellLocIds_[i][j].size();
	    int numBdryDofs = fvidx_[j].size();
	    for (int k = 0; k < numCellsSide; ++k) {
	      int cidx = bdryCellLocIds_[i][j][k];
	      for (int l = 0; l < numBdryDofs; ++l) {
		//std::cout << "\n   j=" << j << "  l=" << l << "  " << fidx[j][l];
		for (int m=0; m < d; ++m) {
		  (*L[m])(cidx,fvidx_[j][l]) = static_cast<Real>(0);
		}
	      }
	    }
	    numBdryDofs = fpidx_[j].size();
	    for (int k = 0; k < numCellsSide; ++k) {
	      int cidx = bdryCellLocIds_[i][j][k];
	      for (int l = 0; l < numBdryDofs; ++l) {
		(*L[d])(cidx,fpidx_[j][l]) = static_cast<Real>(0);
	      }
	    }
	  }
        }
      }
    }
    // int numSideSets = bdryCellLocIds_.size();
    // if (numSideSets > 0) {
    //   for (int i = 0; i < numSideSets; ++i) {
    // 	if (i!=9) {
    // 	  int numLocalSideIds = bdryCellLocIds_[i].size();
    // 	  for (int j = 0; j < numLocalSideIds; ++j) {
    // 	    int numCellsSide = bdryCellLocIds_[i][j].size();
    // 	    int numBdryDofs = fvidx_[j].size();
    // 	    for (int k = 0; k < numCellsSide; ++k) {
    // 	      int cidx = bdryCellLocIds_[i][j][k];
    // 	      for (int l = 0; l < numBdryDofs; ++l) {
    // 		//std::cout << "\n   j=" << j << "  l=" << l << "  " << fidx[j][l];
    // 		for (int m=0; m < d; ++m) {
    // 		  (*L[m])(cidx,fvidx_[j][l]) = static_cast<Real>(0);
    // 		}
    // 	      }
    // 	    }
    // 	  }
    //     }
    //     // Pressure pinning
    //     if (i==9) {
    //       int j = 2, l = 0;
    //       if (bdryCellLocIds_[i][0].size() > 0) { 
    //         int cidx = bdryCellLocIds_[i][0][0];
    //         (*L[d])(cidx,fpidx_[j][l]) = static_cast<Real>(0);
    //       }
    //     }
    //   }
    // }

    // Evaluate/interpolate finite element fields on cells.
    ROL::Ptr<Intrepid::FieldContainer<Real> > valVelX_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valVelY_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valVelXPhi_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valVelYPhi_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv, p);
    feVel_->evaluateValue(valVelX_eval, L[0]);
    feVel_->evaluateValue(valVelY_eval, L[1]);

    // Compute nonlinear terms in the Navier-Stokes equations.
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*valVelXPhi_eval, *valVelX_eval, *(feVel_->N()));
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*valVelYPhi_eval, *valVelY_eval, *(feVel_->N()));

    Intrepid::RealSpaceTools<Real>::scale(*valVelXPhi_eval,advection_coeff_);
    Intrepid::RealSpaceTools<Real>::scale(*valVelYPhi_eval,advection_coeff_);

    /*** Evaluate weak form of the Hessian. ***/
    // X component of velocity equation.
    Intrepid::FunctionSpaceTools::integrate<Real>(velXvelX_jac,
                                                  *valVelXPhi_eval,        // L Phi
                                                  *(feVel_->DNDdetJ(0)),   // dPhi/dx
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>(velXvelX_jac,
                                                  *(feVel_->DNDdetJ(0)),   // dPhi/dx
                                                  *valVelXPhi_eval,        // L Phi
                                                  Intrepid::COMP_CPP,
                                                  true);
    Intrepid::FunctionSpaceTools::integrate<Real>(velXvelY_jac,
                                                  *(feVel_->DNDdetJ(1)),   // dPhi/dy
                                                  *valVelXPhi_eval,        // L Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>(velXvelY_jac,
                                                  *valVelYPhi_eval,        // L Phi
                                                  *(feVel_->DNDdetJ(0)),   // dPhi/dx
                                                  Intrepid::COMP_CPP,
                                                  true);
    // Y component of velocity equation.
    Intrepid::FunctionSpaceTools::integrate<Real>(velYvelY_jac,
                                                  *valVelYPhi_eval,        // L Phi
                                                  *(feVel_->DNDdetJ(1)),   // dPhi/dy
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>(velYvelY_jac,
                                                  *(feVel_->DNDdetJ(1)),   // dPhi/dy
                                                  *valVelYPhi_eval,        // L Phi
                                                  Intrepid::COMP_CPP,
                                                  true);
    Intrepid::FunctionSpaceTools::integrate<Real>(velYvelX_jac,
                                                  *(feVel_->DNDdetJ(0)),   // dPhi/dx
                                                  *valVelYPhi_eval,        // L Phi
                                                  Intrepid::COMP_CPP,
                                                  false);
    Intrepid::FunctionSpaceTools::integrate<Real>(velYvelX_jac,
                                                  *valVelXPhi_eval,        // L Phi
                                                  *(feVel_->DNDdetJ(1)),   // dPhi/dy
                                                  Intrepid::COMP_CPP,
                                                  true);

    // Combine the Hessians.
    fieldHelper_->combineFieldCoeff(hess, J);

  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Stokes::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Stokes::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Stokes::Hessian_22): Hessian is zero.");
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
    // Construct control boundary FE
    int sideset = 6;
    int numLocSides = bdryCellNodes[sideset].size();
    feVelBdry_.resize(numLocSides);
    for (int j = 0; j < numLocSides; ++j) {
      if (bdryCellNodes[sideset][j] != ROL::nullPtr) {
        feVelBdry_[j] = ROL::makePtr<FE<Real>>(bdryCellNodes[sideset][j],basisPtrVel_,bdryCub_,j);
      }
    }
    computeControlWeight();
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

  const std::vector<ROL::Ptr<FE<Real> > > getVelocityBdryFE(void) const {
    return feVelBdry_;
  }

  const std::vector<std::vector<int> > getBdryCellLocIds(const int sideset = -1) const {
    int side = sideset;
    if ( sideset < 0 ) {
      side = 6;
    }
    return bdryCellLocIds_[side];
  }

  const ROL::Ptr<FieldHelper<Real> > getFieldHelper(void) const {
    return fieldHelper_;
  }

}; // PDE_Stokes

#endif
