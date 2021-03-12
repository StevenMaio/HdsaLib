#ifndef PDE_BRANDES_GRIESSE_HPP
#define PDE_BRANDES_GRIESSE_HPP

#include "../../../../PDE-OPT/TOOLS/pde.hpp"
#include "../../../../PDE-OPT/TOOLS/fe.hpp"

#include "Intrepid_HGRAD_LINE_Cn_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C2_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

template <class Real>
class PDE_Brandes_Griesse : public PDE<Real> {
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
  ROL::Ptr<FE<Real> > fe_vol_;
  std::vector<std::vector<ROL::Ptr<FE<Real> > > > feBdry_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int> > fidx_;
  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)

  bool useStateRiesz_;
  bool useControlRiesz_;

  Real AlphaFunc(const std::vector<Real> & coords, const std::vector<Real> & param) {
    Real a(0);
      if(coords[0]>-0.75 && coords[0]<0.0 && coords[1]==1.0)
  	{
  	  a = param[0];
  	}
      else if(coords[0]>-0.75 && coords[0]<0.75 && coords[1]==-1.0)
  	{
  	  a = param[1];
  	}
      else
  	{
  	  a = 0.0;
  	}
    return a;
  }

  ROL::Ptr<Intrepid::FieldContainer<Real> > ctrlWeight_;

  void computeControlWeight(void) {
    int c = fe_vol_->cubPts()->dimension(0);
    int p = fe_vol_->cubPts()->dimension(1);
    int d = fe_vol_->cubPts()->dimension(2);

    ctrlWeight_ = ROL::makePtr<Intrepid::FieldContainer<Real> >(c,p);
 
    const Real zero(0), one(1);
    bool inside(false);
    std::vector<Real> x(d);
    for (int i = 0; i < c; ++i) {
      inside = false;
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          x[k] = (*fe_vol_->cubPts())(i,j,k);
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
    if(x[0]>-0.75 && x[0]<0.75 && x[1]>-0.8 && x[1]<-0.6)
      {
	inside = true;
      }
    else if(x[0]>-0.8 && x[0]<0.0 && x[1]>0.4 && x[1]<0.8)
      {
    	inside = true;
      }
    else
      {
	inside = false;
      }
    return inside;
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

  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > Generate_Boundary_Alpha_Values( )
  {
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellAlphaValues;
    std::vector<Real> param = PDE<Real>::getParameter();
    if(param.size() == 0)
      {
	param.resize(2);
	param[0] = 2.0;
	param[1] = 1.0;
      }
    
    // Compute alpha values at DOFs.
    int d = basisPtr_->getBaseCellTopology().getDimension();
    int numSidesets = bdryCellLocIds_.size();
    bdryCellAlphaValues.resize(numSidesets);
    for (int i=0; i<numSidesets; ++i) {
      int numLocSides = bdryCellLocIds_[i].size();
      bdryCellAlphaValues[i].resize(numLocSides);
      for (int j=0; j<numLocSides; ++j) {
        int c = bdryCellLocIds_[i][j].size();
        int f = basisPtr_->getCardinality();
        bdryCellAlphaValues[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);
        ROL::Ptr<Intrepid::FieldContainer<Real> > coords =
          ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f, d);
        if (c > 0) {
          fe_vol_->computeDofCoords(coords, bdryCellNodes_[i][j]);
        }
        for (int k=0; k<c; ++k) {
          for (int l=0; l<f; ++l) {
            std::vector<Real> dofpoint(d);
            for (int m=0; m<d; ++m) {
              dofpoint[m] = (*coords)(k, l, m);
            }
	    (*bdryCellAlphaValues[i][j])(k, l) = AlphaFunc(dofpoint,param);
          }
        }
      }
    }
    return bdryCellAlphaValues;
  }

public:
  PDE_Brandes_Griesse(Teuchos::ParameterList &parlist) {
    // Finite element fields.
    int basisOrder = parlist.sublist("Problem").get("Basis Order",1);
    int cubDegree  = parlist.sublist("Problem").get("Cubature Degree",4);
    if (basisOrder > 2 || basisOrder < 1) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> PDE-OPT/poisson/pde_poisson.hpp: Basis order is not 1 or 2!");
    }
  
    if (basisOrder == 1) {
      basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> > >();
    }
    else if (basisOrder == 2) {
      basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real> > >();
    }

    basisPtrs_.clear(); basisPtrs_.push_back(basisPtr_);
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();
    Intrepid::DefaultCubatureFactory<Real> cubFactory;
    cellCub_ = cubFactory.create(cellType, cubDegree);
    // Problem data.
    useStateRiesz_   = parlist.sublist("Problem").get("Use State Riesz Map", true);      // use Riesz map for state variables?
    useControlRiesz_ = parlist.sublist("Problem").get("Use Control Riesz Map", true);    // use Riesz map for control variables?

    int d = cellType.getDimension();
    shards::CellTopology bdryCellType = cellType.getCellTopologyData(d-1, 0);
    int bdryCubDegree = parlist.sublist("Problem").get("Boundary Cubature Degree",2); // set cubature degree, e.g., 2
    bdryCub_ = cubFactory.create(bdryCellType, bdryCubDegree);
  }

  void residual(ROL::Ptr<Intrepid::FieldContainer<Real> > & res,
                const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellAlphaValues = Generate_Boundary_Alpha_Values();

    // GET DIMENSIONS
    int c = u_coeff->dimension(0);
    int p = cellCub_->getNumPoints();
    int f = basisPtr_->getCardinality();
    int d = cellCub_->getDimension();
    // INITIALIZE RESIDUAL
    res = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);
    // COMPUTE STIFFNESS TERM
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradU_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p, d);
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    Intrepid::FunctionSpaceTools::integrate<Real>(*res, *gradU_eval, *(fe_vol_->gradNdetJ()), Intrepid::COMP_CPP, false);
    // ADD CONTROL TERM TO RESIDUAL
    if ( z_coeff != ROL::nullPtr ) {
      ROL::Ptr<Intrepid::FieldContainer<Real> > valZ_eval =
        ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
      fe_vol_->evaluateValue(valZ_eval, z_coeff);
      Intrepid::RealSpaceTools<Real>::scale(*valZ_eval,static_cast<Real>(-1.0));

      ROL::Ptr<Intrepid::FieldContainer<Real> > weight_valZ_eval = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, p);
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*weight_valZ_eval, *ctrlWeight_, *valZ_eval);
      Intrepid::FunctionSpaceTools::integrate<Real>(*res, *weight_valZ_eval, *(fe_vol_->NdetJ()), Intrepid::COMP_CPP, true);
    }
    // APPLY ROBIN CONDITIONS
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      const int numCubPerSide = bdryCub_->getNumPoints();
      for (int i = 0; i < numSideSets; ++i) {
        int numLocalSideIds = bdryCellLocIds_[i].size();
        for (int j = 0; j < numLocalSideIds; ++j) {
          int numCellsSide = bdryCellLocIds_[i][j].size();
          int numBdryDofs = fidx_[j].size();

    	  if (numCellsSide) {

    	  ROL::Ptr<Intrepid::FieldContainer<Real> > robinRes
    	    = ROL::makePtr<Intrepid::FieldContainer<Real> >(numCellsSide, f);

    	  // Get u coefficients on Robin boundary
    	  ROL::Ptr<Intrepid::FieldContainer<Real> > u_coeff_bdry
            = getBoundaryCoeff(*u_coeff, i, j);    	  
    	  // Evaluate u on FE basis
    	  ROL::Ptr<Intrepid::FieldContainer<Real> > valU_eval_bdry
    	    = ROL::makePtr<Intrepid::FieldContainer<Real> >(numCellsSide, numCubPerSide);
    	  feBdry_[i][j]->evaluateValue(valU_eval_bdry, u_coeff_bdry);
    	  ROL::Ptr<Intrepid::FieldContainer<Real> > valAlpha_eval_bdry
    	    = ROL::makePtr<Intrepid::FieldContainer<Real> >(numCellsSide, numCubPerSide);
    	  feBdry_[i][j]->evaluateValue(valAlpha_eval_bdry, bdryCellAlphaValues[i][j]);
    	  ROL::Ptr<Intrepid::FieldContainer<Real> > valUalpha_eval_bdry
    	    = ROL::makePtr<Intrepid::FieldContainer<Real> >(numCellsSide, numCubPerSide);
    	  Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*valUalpha_eval_bdry, *valU_eval_bdry, *valAlpha_eval_bdry);
  
    	  // Compute Neumann residual
          Intrepid::FunctionSpaceTools::integrate<Real>(*robinRes,
                                                        *valUalpha_eval_bdry,
                                                        *(feBdry_[i][j]->NdetJ()),
                                                        Intrepid::COMP_CPP, false);
	  

          for (int k = 0; k < numCellsSide; ++k) 
    	    {
    	      int cidx = bdryCellLocIds_[i][j][k];
    	      for (int l = 0; l < numBdryDofs; ++l) 
    		{
    		  (*res)(cidx,fidx_[j][l]) += (*robinRes)(k,fidx_[j][l]);
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
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellAlphaValues = Generate_Boundary_Alpha_Values();
    
    // GET DIMENSIONS
    int c = u_coeff->dimension(0);
    int f = basisPtr_->getCardinality();
    // INITILAIZE JACOBIAN
    jac = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f, f);
    // COMPUTE STIFFNESS TERM
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac, *(fe_vol_->gradN()), *(fe_vol_->gradNdetJ()), Intrepid::COMP_CPP, false);
    // APPLY ROBIN CONDITIONS
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
    		(*jac)(cidx,fidx_[j][l],m) += (*bdryCellAlphaValues[i][j])(k,fidx_[j][l])*(*feBdry_[i][j]->massMat())(k,fidx_[j][l],m);
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
    if ( z_coeff != ROL::nullPtr ) {
      // GET DIMENSIONS
      int c = u_coeff->dimension(0);
      int f = basisPtr_->getCardinality();
      // INITIALIZE JACOBIAN
      jac = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f, f);
      // ADD CONTROL TERM
      int p = fe_vol_->gradN()->dimension(2);
      ROL::Ptr<Intrepid::FieldContainer<Real> > F = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f, p);
      Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*F, *ctrlWeight_, *fe_vol_->N());
      Intrepid::FunctionSpaceTools::integrate<Real>(*jac, *F, *(fe_vol_->NdetJ()), Intrepid::COMP_CPP, false);
      Intrepid::RealSpaceTools<Real>::scale(*jac,static_cast<Real>(-1.0));
    }
    else {
      throw Exception::Zero(">>> (PDE_Brandes_Griesse::Jacobian_2): Jacobian is zero.");
    }
  }

  void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Brandes_Griesse::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Brandes_Griesse::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Brandes_Griesse::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_Brandes_Griesse::Hessian_22): Hessian is zero.");
  }

  void RieszMap_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & riesz) {
    // Optionally disable Riesz map ...
    if (!useStateRiesz_) {
      throw Exception::NotImplemented(">>> (PDE_Brandes_Griesse::RieszMap_1): Not implemented.");
    }

    // ...otherwise ...

    int c = fe_vol_->N()->dimension(0);
    int f = fe_vol_->N()->dimension(1);
    // INITILAIZE JACOBIAN
    riesz = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f, f);
    *riesz = *fe_vol_->stiffMat();
    Intrepid::RealSpaceTools<Real>::add(*riesz,*(fe_vol_->massMat()));
  }

  void RieszMap_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & riesz) {
    // Optionally disable Riesz map ...
    if (!useControlRiesz_) {
      throw Exception::NotImplemented(">>> (PDE_Brandes_Griesse::RieszMap_2): Not implemented.");
    }

    // ...otherwise ...

    int c = fe_vol_->N()->dimension(0);
    int f = fe_vol_->N()->dimension(1);
    // INITILAIZE JACOBIAN
    riesz = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f, f);
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
    fe_vol_ = ROL::makePtr<FE<Real> >(volCellNodes_,basisPtr_,cellCub_);
    // Set local boundary DOFs.
    fidx_ = fe_vol_->getBoundaryDofs();
    // Compute alpha values at DOFs.
    int d = basisPtr_->getBaseCellTopology().getDimension();
    int numSidesets = bdryCellLocIds_.size();

    feBdry_.resize(numSidesets);
    for (int i=0; i<numSidesets; i++){
      int numLocSides = bdryCellNodes_[i].size();
      feBdry_[i].resize(numLocSides);
      for (int j=0; j<numLocSides; j++){
	if (bdryCellNodes_[i][j] != ROL::nullPtr) {
	  feBdry_[i][j] = ROL::makePtr<FE<Real> >(bdryCellNodes_[i][j],basisPtr_,bdryCub_,j);
	}
      }
    }

    computeControlWeight();
  }

  const ROL::Ptr<FE<Real> > getFE(void) const {
    return fe_vol_;
  }

  const std::vector<std::vector<ROL::Ptr<FE<Real> > > >  getFEbdry(void) const {
    return feBdry_;
  }

}; // PDE_Brandes_Griesse

#endif


