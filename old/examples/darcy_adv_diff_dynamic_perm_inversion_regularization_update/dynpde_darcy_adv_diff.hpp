#ifndef PDEOPT_DYNAMIC_DARCY_ADV_DIFF_HPP
#define PDEOPT_DYNAMIC_DARCY_ADV_DIFF_HPP

#include "../../../PDE-OPT/TOOLS/dynpde.hpp"
#include "../../../PDE-OPT/TOOLS/fe.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "ROL_Ptr.hpp"

#include "pde_darcy_adv_diff.hpp"

template <class Real>
class DynamicPDE_darcy_adv_diff : public DynamicPDE<Real> {
private:
  // Cell node information
  std::vector<std::vector<std::vector<int> > > bdryCellLocIds_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellNodes_;
  // Finite element definition
  ROL::Ptr<FE<Real> > fePrs_;
  ROL::Ptr<FE<Real> > feCntm_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int> > fpidx_;
  std::vector<std::vector<int> > fcidx_;

  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellPDofValues_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellCDofValues_;

  // Field pattern, offsets, etc.
  std::vector<std::vector<int> > fieldPattern_;  // local Field/DOF pattern; set from DOF manager 
  int numFields_;                                // number of fields (equations in the PDE)
  int numDofs_;                                  // total number of degrees of freedom for all (local) fields
  std::vector<int> offset_;                      // for each field, a counting offset
  std::vector<int> numFieldDofs_;                // for each field, number of degrees of freedom

  ROL::Ptr<FieldHelper<Real> > fieldHelper_;

  // Steady PDE without Dirichlet BC
  ROL::Ptr<PDE_darcy_adv_diff<Real> > pde_;
  Real theta_;

  Real a;
  int L;
  std::vector<Real> uncertain_basis_grid;

public:
  DynamicPDE_darcy_adv_diff(Teuchos::ParameterList &parlist) {
    pde_ = ROL::makePtr<PDE_darcy_adv_diff<Real>>(parlist);
    // Time-dependent coefficients
    theta_  = parlist.sublist("Time Discretization").get("Theta",1.0);

    L = parlist.sublist("Problem").get("Number of Uncertainty Basis Function", 10);
    a = parlist.sublist("Problem").get("Noise Level", .2);
    uncertain_basis_grid.resize(L+1);
    for(int i = 0; i < L+1; i++)
      {
	uncertain_basis_grid[i] = static_cast<Real>(i)/static_cast<Real>(L);
      }
  }

  void residual(ROL::Ptr<Intrepid::FieldContainer<Real>> & res,
                const ROL::TimeStamp<Real> & ts,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const Real one(1);
    // Retrieve dimensions.
    int c = fePrs_->gradN()->dimension(0);
    int f = fePrs_->gradN()->dimension(1);
    int p = fePrs_->gradN()->dimension(2);
 
    // GET TIME STEP INFORMATION
    Real told = ts.t[0], tnew = ts.t[1], dt = tnew-told;

    updateDirichlet();

    // COMPUTE OLD RESIDUAL
    ROL::Ptr<Intrepid::FieldContainer<Real> > rold;
    pde_->setTime(told);
    pde_->residual(rold,uo_coeff,z_coeff,z_param);

    // COMPUTE NEW RESIDUAL
    ROL::Ptr<Intrepid::FieldContainer<Real> > rnew;
    pde_->setTime(tnew);
    pde_->residual(rnew,un_coeff,z_coeff,z_param);

    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > Rold;
    fieldHelper_->splitFieldCoeff(Rold, rold);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > Rnew;
    fieldHelper_->splitFieldCoeff(Rnew, rnew);
    
    Intrepid::RealSpaceTools<Real>::scale(*Rold[1], (one-theta_)*dt);
    Intrepid::RealSpaceTools<Real>::scale(*Rnew[1], theta_*dt);

    // Integrate Uold * N
    // Split uo_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > Uold;
    fieldHelper_->splitFieldCoeff(Uold, uo_coeff);  
    // Compute Cntm
    ROL::Ptr<Intrepid::FieldContainer<Real> > Cntm_old_eval;
    Cntm_old_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    feCntm_->evaluateValue(Cntm_old_eval, Uold[1]);
    ROL::Ptr<Intrepid::FieldContainer<Real> > Cold = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);
    Intrepid::FunctionSpaceTools::integrate<Real>(*Cold,
                                                  *Cntm_old_eval,
                                                  *(feCntm_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);

    // Integrate Unew * N
    // Split un_coeff into components.
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > Unew;
    fieldHelper_->splitFieldCoeff(Unew, un_coeff);  
    // Compute Cntm
    ROL::Ptr<Intrepid::FieldContainer<Real> > Cntm_new_eval;
    Cntm_new_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    feCntm_->evaluateValue(Cntm_new_eval, Unew[1]);
    ROL::Ptr<Intrepid::FieldContainer<Real> > Cnew = ROL::makePtr<Intrepid::FieldContainer<Real> >(c, f);
    Intrepid::FunctionSpaceTools::integrate<Real>(*Cnew,
                                                  *Cntm_new_eval,
                                                  *(feCntm_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);

    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > R(2);
    R[0] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    R[1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f); 

    Intrepid::RealSpaceTools<Real>::add(*R[0], *Rnew[0]);

    Intrepid::RealSpaceTools<Real>::add(*R[1], *Cnew);
    Intrepid::RealSpaceTools<Real>::subtract(*R[1], *Cold);
    Intrepid::RealSpaceTools<Real>::add(*R[1], *Rold[1]);
    Intrepid::RealSpaceTools<Real>::add(*R[1], *Rnew[1]);
    
    // APPLY DIRICHLET CONDITIONS
    int numLocalSideIds = bdryCellLocIds_[0].size();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[0][j].size();
      int numBdryDofs = fpidx_[j].size();
      for (int k = 0; k < numCellsSide; ++k) {
        int cidx = bdryCellLocIds_[0][j][k];
        for (int l = 0; l < numBdryDofs; ++l) {
          (*R[0])(cidx,fpidx_[j][l])
            = (*Unew[0])(cidx,fpidx_[j][l]) - (*bdryCellPDofValues_[0][j])(k,fpidx_[j][l]);
        }
      }
    }
    
    // APPLY DIRICHLET CONDITIONS
    numLocalSideIds = bdryCellLocIds_[2].size();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[2][j].size();
      int numBdryDofs = fcidx_[j].size();
      for (int k = 0; k < numCellsSide; ++k) {
	int cidx = bdryCellLocIds_[2][j][k];
	for (int l = 0; l < numBdryDofs; ++l) {
	  (*R[1])(cidx,fcidx_[j][l])
	    = (*Unew[1])(cidx,fcidx_[j][l]) - (*bdryCellCDofValues_[2][j])(k,fcidx_[j][l]);
        }
      }
    }

    // Combine the residuals.
    fieldHelper_->combineFieldCoeff(res, R);

  }

  void Jacobian_uo(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                   const ROL::TimeStamp<Real> & ts,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                   const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const Real one(1);
    // Retrieve dimensions.
    int c = fePrs_->gradN()->dimension(0);
    int f = fePrs_->gradN()->dimension(1);
    // GET TIME STEP INFORMATION
    Real told = ts.t[0], tnew = ts.t[1], dt = tnew-told;
    // COMPUTE OLD RESIDUAL
    ROL::Ptr<Intrepid::FieldContainer<Real> > j;
    pde_->setTime(told);
    pde_->Jacobian_1(j,uo_coeff,z_coeff,z_param);
    Intrepid::RealSpaceTools<Real>::scale(*j, (one-theta_)*dt);
    // Initialize jacobians.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > J;
    fieldHelper_->splitFieldCoeff(J, j); 
    J[0][0]->initialize();
    J[0][1]->initialize();

    ROL::Ptr<Intrepid::FieldContainer<Real> > Nold = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);

    // Integrate N * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*Nold,
                                                  *(feCntm_->N()),
                                                  *(feCntm_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);

    Intrepid::RealSpaceTools<Real>::subtract(*J[1][1], *Nold);

    // APPLY DIRICHLET CONDITIONS
    int numLocalSideIds = bdryCellLocIds_[0].size();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[0][j].size();
      int numBdryDofs = fpidx_[j].size();
      for (int k = 0; k < numCellsSide; ++k) {
        int cidx = bdryCellLocIds_[0][j][k];
        for (int l = 0; l < numBdryDofs; ++l) {
          for (int m = 0; m < f; ++m) {
            (*J[0][0])(cidx,fpidx_[j][l],m) = static_cast<Real>(0);
          }
          (*J[0][0])(cidx,fpidx_[j][l],fpidx_[j][l]) = static_cast<Real>(0); 
        }
      }
    }

    // APPLY DIRICHLET CONDITIONS
    numLocalSideIds = bdryCellLocIds_[2].size();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[2][j].size();
      int numBdryDofs = fcidx_[j].size();
      for (int k = 0; k < numCellsSide; ++k) {
        int cidx = bdryCellLocIds_[2][j][k];
        for (int l = 0; l < numBdryDofs; ++l) {
          for (int m = 0; m < f; ++m) {
            (*J[1][1])(cidx,fcidx_[j][l],m) = static_cast<Real>(0);
          }
          (*J[1][1])(cidx,fcidx_[j][l],fcidx_[j][l]) = static_cast<Real>(0); 
        }
      }
    }

    // APPLY DIRICHLET CONDITIONS
    numLocalSideIds = bdryCellLocIds_[2].size();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[2][j].size();
      int numBdryDofs = fcidx_[j].size();
      for (int k = 0; k < numCellsSide; ++k) {
        int cidx = bdryCellLocIds_[2][j][k];
        for (int l = 0; l < numBdryDofs; ++l) {
          for (int m = 0; m < f; ++m) {
            (*J[1][0])(cidx,fcidx_[j][l],m) = static_cast<Real>(0);
          }
          (*J[1][0])(cidx,fcidx_[j][l],fcidx_[j][l]) = static_cast<Real>(0);
        }
      }
    }
    
    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);

  }

  void Jacobian_un(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                   const ROL::TimeStamp<Real> & ts,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                   const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Retrieve dimensions.
    int c = fePrs_->gradN()->dimension(0);
    int f = fePrs_->gradN()->dimension(1);
    // GET TIME STEP INFORMATION
    Real told = ts.t[0], tnew = ts.t[1], dt = tnew-told;
    // COMPUTE NEW RESIDUAL
    ROL::Ptr<Intrepid::FieldContainer<Real> > j;
    pde_->setTime(tnew);
    pde_->Jacobian_1(j,un_coeff,z_coeff,z_param);

    // Initialize jacobians.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > J;
    fieldHelper_->splitFieldCoeff(J, j);
    Intrepid::RealSpaceTools<Real>::scale(*J[1][0], theta_*dt);
    Intrepid::RealSpaceTools<Real>::scale(*J[1][1], theta_*dt);

    // Integrate N * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[1][1],
                                                  *(feCntm_->N()),
                                                  *(feCntm_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);


    // APPLY DIRICHLET CONDITIONS
    int numLocalSideIds = bdryCellLocIds_[0].size();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[0][j].size();
      int numBdryDofs = fpidx_[j].size();
      for (int k = 0; k < numCellsSide; ++k) {
        int cidx = bdryCellLocIds_[0][j][k];
        for (int l = 0; l < numBdryDofs; ++l) {
          for (int m = 0; m < f; ++m) {
            (*J[0][0])(cidx,fpidx_[j][l],m) = static_cast<Real>(0);
          }
          (*J[0][0])(cidx,fpidx_[j][l],fpidx_[j][l]) = static_cast<Real>(1);
        }
      }
    }

    // APPLY DIRICHLET CONDITIONS
    numLocalSideIds = bdryCellLocIds_[2].size();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[2][j].size();
      int numBdryDofs = fcidx_[j].size();
      for (int k = 0; k < numCellsSide; ++k) {
        int cidx = bdryCellLocIds_[2][j][k];
        for (int l = 0; l < numBdryDofs; ++l) {
          for (int m = 0; m < f; ++m) {
            (*J[1][1])(cidx,fcidx_[j][l],m) = static_cast<Real>(0);
          }
          (*J[1][1])(cidx,fcidx_[j][l],fcidx_[j][l]) = static_cast<Real>(1);
        }
      }
    }

    // APPLY DIRICHLET CONDITIONS
    numLocalSideIds = bdryCellLocIds_[2].size();
    for (int j = 0; j < numLocalSideIds; ++j) {
      int numCellsSide = bdryCellLocIds_[2][j].size();
      int numBdryDofs = fcidx_[j].size();
      for (int k = 0; k < numCellsSide; ++k) {
        int cidx = bdryCellLocIds_[2][j][k];
        for (int l = 0; l < numBdryDofs; ++l) {
          for (int m = 0; m < f; ++m) {
            (*J[1][0])(cidx,fcidx_[j][l],m) = static_cast<Real>(0);
          }
          (*J[1][0])(cidx,fcidx_[j][l],fcidx_[j][l]) = static_cast<Real>(0);
        }
      }
    }

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);
  }

  void Jacobian_zf(ROL::Ptr<Intrepid::FieldContainer<Real>> & jac,
                   const ROL::TimeStamp<Real> & ts,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                   const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const Real one(1);
    // GET DIMENSIONS
    int c = fePrs_->gradN()->dimension(0);
    int f = fePrs_->gradN()->dimension(1);
    // GET TIME STEP INFORMATION
    Real told = ts.t[0], tnew = ts.t[1], dt = tnew-told;
    // INITILAIZE JACOBIAN
    ROL::Ptr<Intrepid::FieldContainer<Real> > jo;
    ROL::Ptr<Intrepid::FieldContainer<Real> > jn;
    // ASSEMBLE OLD TIME JACOBIAN
    pde_->setTime(told);
    pde_->Jacobian_2(jo,uo_coeff,z_coeff,z_param);
    // ASSEMBLE NEW TIME JACOBIAN
    pde_->setTime(tnew);
    pde_->Jacobian_2(jn,un_coeff,z_coeff,z_param);

    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > Jo;
    fieldHelper_->splitFieldCoeff(Jo, jo);
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > Jn;
    fieldHelper_->splitFieldCoeff(Jn, jn);

    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > J(2);
    J[0].resize(2);
    J[1].resize(2);
    J[0][0] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    J[0][1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    J[1][0] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    J[1][1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);

    Intrepid::RealSpaceTools<Real>::add(*J[0][0], *Jn[0][0]);
    Intrepid::RealSpaceTools<Real>::scale(*Jo[1][0],(one-theta_)*dt);
    Intrepid::RealSpaceTools<Real>::scale(*Jn[1][0],theta_*dt);
    Intrepid::RealSpaceTools<Real>::add(*J[1][0], *Jo[1][0]);
    Intrepid::RealSpaceTools<Real>::add(*J[1][0], *Jn[1][0]);

    // APPLY DIRICHLET CONDITIONS
    int numSideSets = bdryCellLocIds_[0].size();
    if (numSideSets > 0) {
      int numLocalSideIds = bdryCellLocIds_[0].size();
      for (int j = 0; j < numLocalSideIds; ++j) {
    	int numCellsSide = bdryCellLocIds_[0][j].size();
    	int numBdryDofs = fpidx_[j].size();
    	for (int k = 0; k < numCellsSide; ++k) {
    	  int cidx = bdryCellLocIds_[0][j][k];
    	  for (int l = 0; l < numBdryDofs; ++l) {
    	    //std::cout << "\n   j=" << j << "  l=" << l << "  " << fidx[j][l];
    	    for (int m = 0; m < f; ++m) {
    	      (*J[0][0])(cidx,fpidx_[j][l],m) = static_cast<Real>(0);
    	    }
    	  }
    	}
      }
    }

      // APPLY DIRICHLET CONDITIONS
    numSideSets = bdryCellLocIds_[0].size();
    if (numSideSets > 0) {
      int numLocalSideIds = bdryCellLocIds_[2].size();
      for (int j = 0; j < numLocalSideIds; ++j) {
	int numCellsSide = bdryCellLocIds_[2][j].size();
	int numBdryDofs = fcidx_[j].size();
	for (int k = 0; k < numCellsSide; ++k) {
	  int cidx = bdryCellLocIds_[2][j][k];
	  for (int l = 0; l < numBdryDofs; ++l) {
	    //std::cout << "\n   j=" << j << "  l=" << l << "  " << fidx[j][l];
	    for (int m = 0; m < f; ++m) {
	      (*J[1][0])(cidx,fcidx_[j][l],m) = static_cast<Real>(0);
	    }
	    (*J[1][0])(cidx,fcidx_[j][l],fcidx_[j][l]) = static_cast<Real>(0);
	  }
	}
      }
    }

      // Combine the jacobians.
      fieldHelper_->combineFieldCoeff(jac, J);
  }

  void Jacobian_zp(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & jac,
                   const ROL::TimeStamp<Real> & ts,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                   const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_darcy_adv_diff::Jacobian_zp): Jacobian is zero.");
  }

  void Hessian_uo_uo(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> (PDE_darcy_adv_diff::Hessian_uo_uo): Not implemented.");
  }

  void Hessian_uo_un(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_darcy_adv_diff::Hessian_uo_un): Hessian is zero.");
  }

  void Hessian_uo_zf(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> (PDE_darcy_adv_diff::Hessian_uo_zf): Not implemented.");
  }

  void Hessian_uo_zp(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_darcy_adv_diff::Hessian_uo_zp): Hessian is zero.");
  }

  void Hessian_un_uo(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_darcy_adv_diff::Hessian_un_uo): Hessian is zero.");
  }

  void Hessian_un_un(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> (PDE_darcy_adv_diff::Hessian_un_un): Not implemented.");
  }

  void Hessian_un_zf(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> (PDE_darcy_adv_diff::Hessian_un_zf): Not implemented.");
  }

  void Hessian_un_zp(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_darcy_adv_diff::Hessian_un_zp): Hessian is zero.");
  }

  void Hessian_zf_uo(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> (PDE_darcy_adv_diff::Hessian_zf_uo): Not implemented.");
  }

  void Hessian_zf_un(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> (PDE_darcy_adv_diff::Hessian_zf_un): Not implemented.");
  }

  void Hessian_zf_zf(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::NotImplemented(">>> (PDE_darcy_adv_diff::Hessian_zf_zf): Not implemented.");
  }

  void Hessian_zf_zp(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_darcy_adv_diff::Hessian_zf_zp): Hessian is zero.");
  }

  void Hessian_zp_uo(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_darcy_adv_diff::Hessian_zp_uo): Hessian is zero.");
  }

  void Hessian_zp_un(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_darcy_adv_diff::Hessian_zp_un): Hessian is zero.");
  }

  void Hessian_zp_zf(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_darcy_adv_diff::Hessian_zp_zf): Hessian is zero.");
  }

  void Hessian_zp_zp(std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> & hess,
                     const ROL::TimeStamp<Real> & ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> (PDE_darcy_adv_diff::Hessian_zp_zp): Hessian is zero.");
  }

  void RieszMap_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & riesz) {
    pde_->RieszMap_1(riesz);
  }

  void RieszMap_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & riesz) {
    pde_->RieszMap_2(riesz);
  }

  const ROL::Ptr<PDE<Real> > getPDE(void) const {
    return pde_;
  }

  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> getFields() {
    return pde_->getFields();
  }

  const ROL::Ptr<FieldHelper<Real> > getFieldHelper(void) const {
    return fieldHelper_;
  }

  void setParameter(const std::vector<Real> & param)
  {
    DynamicPDE<Real>::setParameter(param);
    pde_->setParameter(param);
  }

  void setCellNodes(const ROL::Ptr<Intrepid::FieldContainer<Real> > &volCellNodes,
                    const std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int> > > &bdryCellLocIds) {
    bdryCellNodes_ = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    pde_->setCellNodes(volCellNodes,bdryCellNodes,bdryCellLocIds);
    // Finite element definition.
    fePrs_ = pde_->getPressureFE();
    feCntm_ = pde_->getContaminantFE();
    // Get boundary degrees of freedom.
    fpidx_ = fePrs_->getBoundaryDofs();
    fcidx_ = feCntm_->getBoundaryDofs();

    computeDirichlet();
  }

  void setFieldPattern(const std::vector<std::vector<int>> &fieldPattern) {
    pde_->setFieldPattern(fieldPattern);
    fieldHelper_ = pde_->getFieldHelper();
  }

  const ROL::Ptr<FE<Real> > getPressureFE(void) const {
    return fePrs_;
  }

  const ROL::Ptr<FE<Real> > getContaminantFE(void) const {
    return feCntm_;
  }

  const std::vector<std::vector<std::vector<int> > > getBdryCellLocIds(void) const {
    return bdryCellLocIds_;
  }


private:

 void computeDirichlet(void) {
    // Compute Dirichlet values at DOFs.
    // Retrieve dimensions.
    int f = fePrs_->gradN()->dimension(1);
    int d = fePrs_->gradN()->dimension(3);
    int numSidesets = bdryCellLocIds_.size();
    bdryCellPDofValues_.resize(numSidesets);
    bdryCellCDofValues_.resize(numSidesets);
    for (int i=0; i<numSidesets; ++i) {
      int numLocSides = bdryCellLocIds_[i].size();
      bdryCellPDofValues_[i].resize(numLocSides);
      bdryCellCDofValues_[i].resize(numLocSides);
      for (int j=0; j<numLocSides; ++j) {
        int c = bdryCellLocIds_[i][j].size();
        bdryCellPDofValues_[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
        bdryCellCDofValues_[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
        ROL::Ptr<Intrepid::FieldContainer<Real> > Pcoords, Ccoords;
        Pcoords = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f,d);
        Ccoords = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f,d);
        if (c > 0) {
          fePrs_->computeDofCoords(Pcoords, bdryCellNodes_[i][j]);
          feCntm_->computeDofCoords(Ccoords, bdryCellNodes_[i][j]);
        }
        for (int k=0; k<c; ++k) {
          for (int l=0; l<f; ++l) {
            std::vector<Real> dofpoint(d);
            for (int m=0; m<d; ++m) {
              dofpoint[m] = (*Pcoords)(k, l, m);
            }
	    (*bdryCellPDofValues_[i][j])(k, l) = pressueDirichletFunc(dofpoint, i, j);
          }

          for (int l=0; l<f; ++l) {
            std::vector<Real> dofpoint(d);
            for (int m=0; m<d; ++m) {
              dofpoint[m] = (*Ccoords)(k, l, m);
            }
            (*bdryCellCDofValues_[i][j])(k, l) = contaminantDirichletFunc(dofpoint, i, j);
          }
        }
      }
    }
  }
  
  void updateDirichlet(void) {
    // Compute Dirichlet values at DOFs.
    // Retrieve dimensions.
    int f = fePrs_->gradN()->dimension(1);
    int d = fePrs_->gradN()->dimension(3);
    int numSidesets = bdryCellLocIds_.size();
    for (int i=0; i<numSidesets; ++i) {
      int numLocSides = bdryCellLocIds_[i].size();
      for (int j=0; j<numLocSides; ++j) {
        int c = bdryCellLocIds_[i][j].size();
	ROL::Ptr<Intrepid::FieldContainer<Real> > Pcoords, Ccoords;
        Pcoords = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f,d);
        Ccoords = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f,d);
        if (c > 0) {
          fePrs_->computeDofCoords(Pcoords, bdryCellNodes_[i][j]);
          feCntm_->computeDofCoords(Ccoords, bdryCellNodes_[i][j]);
        }
        for (int k=0; k<c; ++k) {
          for (int l=0; l<f; ++l) {
            std::vector<Real> dofpoint(d);
            for (int m=0; m<d; ++m) {
              dofpoint[m] = (*Pcoords)(k, l, m);
            }
	    (*bdryCellPDofValues_[i][j])(k, l) = pressueDirichletFunc(dofpoint, i, j);
          }
	  
          for (int l=0; l<f; ++l) {
            std::vector<Real> dofpoint(d);
            for (int m=0; m<d; ++m) {
              dofpoint[m] = (*Ccoords)(k, l, m);
            }
            (*bdryCellCDofValues_[i][j])(k, l) = contaminantDirichletFunc(dofpoint, i, j);
          }
        }
      }
    }
  }
  
  Real pressueDirichletFunc(const std::vector<Real> & x, int sideset, int locSideId) const {
    Real val = 0.0;
    Real pi(M_PI);

    if(x[0] == 0.0)
      {
	val = 10.0 + 2.0*std::cos(2*pi*x[1]); 

	const std::vector<Real> param = DynamicPDE<Real>::getParameter();
	if((int)param.size()>0)
	  {
	    Real val_uncertain = 0.0;
	    int offset = 16*9 + 1;
	    int count = 0;
	    Real basis_fun = 0.0;
	    for(int i = 0; i < L+1; i++)
	      {
		basis_fun = Parameter_Basis_Fun_Eval(x[1],i);
		val_uncertain += param[offset + count]*basis_fun;
		count = count+1;
	      }
	    val_uncertain = a*val_uncertain;
	    val_uncertain = 1+val_uncertain;
	    val = val*val_uncertain;
	  }

      }
    else if(x[0] == 1.0)
      {
	val = 15 + 1.0*std::cos(2*pi*x[1]) + 0.5*std::cos(4*pi*x[1]); 

	const std::vector<Real> param = DynamicPDE<Real>::getParameter();
	if((int)param.size()>0)
	  {
	    Real val_uncertain = 0.0;
	    int offset = 16*9 + 1 + (L+1);
	    int count = 0;
	    Real basis_fun = 0.0;
	    for(int i = 0; i < L+1; i++)
	      {
		basis_fun = Parameter_Basis_Fun_Eval(x[1],i);
		val_uncertain += param[offset + count]*basis_fun;
		count = count+1;
	      }
	    val_uncertain = a*val_uncertain;
	    val_uncertain = 1+val_uncertain;
	    val = val*val_uncertain;
	  }

      }

    return val;
  }

  Real contaminantDirichletFunc(const std::vector<Real> & x, int sideset, int locSideId) const {
    Real val = 0.0;
    return val;
  }

  Real Parameter_Basis_Fun_Eval(const Real & x, int i) const
  {
    // Evaluates ith 1D basis function on a grid of L+1 nodes at point x
    Real val = 0.0;
    Real dist = std::abs(x-uncertain_basis_grid[i]);
    if( dist < 1.0/static_cast<Real>(L) )
      {
	val = 1.0 - dist*static_cast<Real>(L);
      }
    return val;
  }

}; // DynamicPDE_darcy_adv_diff

#endif
