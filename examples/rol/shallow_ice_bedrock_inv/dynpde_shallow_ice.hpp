#ifndef PDEOPT_DYNAMIC_SHALLOW_ICE_HPP
#define PDEOPT_DYNAMIC_SHALLOW_ICE_HPP

#include "../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/dynpde.hpp"
#include "../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/fe.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "ROL_Ptr.hpp"

#include "pde_shallow_ice.hpp"

template <class Real>
class DynamicPDE_shallow_ice : public DynamicPDE<Real>
{
private:
  // Cell node information
  std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> bdryCellNodes_;
  // Finite element definition
  ROL::Ptr<FE<Real>> fe_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int>> fidx_;

  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> bdryCellDofValues_;

  // Field pattern, offsets, etc.
  std::vector<std::vector<int>> fieldPattern_; // local Field/DOF pattern; set from DOF manager
  ROL::Ptr<FieldHelper<Real>> fieldHelper_;

  // Steady PDE without Dirichlet BC
  ROL::Ptr<PDE_shallow_ice<Real>> pde_;
  Real theta_;

  Real width_;
  Real height_;

  int num_coeff_load_;
  std::vector<Real> surface_height_coeff_;
  int N_;
  Real h_;

public:
  DynamicPDE_shallow_ice(Teuchos::ParameterList &parlist)
  {
    pde_ = ROL::makePtr<PDE_shallow_ice<Real>>(parlist);
    // Time-dependent coefficients
    theta_ = parlist.sublist("Time Discretization").get("Theta", 1.0);

    width_ = parlist.sublist("Geometry").get("Width", 1.0);
    height_ = parlist.sublist("Geometry").get("Height", 1.0);

    num_coeff_load_ = parlist.sublist("Problem").get("Number of Coefficients in Loaded Fields", 10);
    N_ = std::sqrt(num_coeff_load_) - 1;
    h_ = 1.0 / static_cast<Real>(N_);
    surface_height_coeff_.resize(num_coeff_load_);
    // read in data
    std::ifstream in("Surface_Height.txt");
    // read the elements in the file into a vector
    // test file open
    if (in)
    {
      for (int j = 0; j < num_coeff_load_; j++)
      {
        in >> surface_height_coeff_[j];
      }
    }
    else
    {
      std::cout << "Error loading the data from Bedrock_Topography.txt" << std::endl;
    }
  }

  void residual(ROL::Ptr<Intrepid::FieldContainer<Real>> &res,
                const ROL::TimeStamp<Real> &ts,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> &uo_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> &un_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    const Real one(1);
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);

    // GET TIME STEP INFORMATION
    Real told = ts.t[0], tnew = ts.t[1], dt = tnew - told;

    // COMPUTE OLD RESIDUAL
    ROL::Ptr<Intrepid::FieldContainer<Real>> rold;
    pde_->setTime(told);
    pde_->residual(rold, uo_coeff, z_coeff, z_param);

    // COMPUTE NEW RESIDUAL
    ROL::Ptr<Intrepid::FieldContainer<Real>> rnew;
    pde_->setTime(tnew);
    pde_->residual(rnew, un_coeff, z_coeff, z_param);

    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> Rold;
    fieldHelper_->splitFieldCoeff(Rold, rold);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> Rnew;
    fieldHelper_->splitFieldCoeff(Rnew, rnew);

    Intrepid::RealSpaceTools<Real>::scale(*Rold[0], (one - theta_) * dt);
    Intrepid::RealSpaceTools<Real>::scale(*Rnew[0], theta_ * dt);

    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> Uold;
    fieldHelper_->splitFieldCoeff(Uold, uo_coeff);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> Unew;
    fieldHelper_->splitFieldCoeff(Unew, un_coeff);

    // Integrate Uold * N
    // Split uo_coeff into components.
    ROL::Ptr<Intrepid::FieldContainer<Real>> old_eval;
    old_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(old_eval, Uold[0]);
    ROL::Ptr<Intrepid::FieldContainer<Real>> U_old = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    Intrepid::FunctionSpaceTools::integrate<Real>(*U_old,
                                                  *old_eval,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);

    // Integrate Unew * N
    ROL::Ptr<Intrepid::FieldContainer<Real>> new_eval;
    new_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(new_eval, Unew[0]);
    ROL::Ptr<Intrepid::FieldContainer<Real>> U_new = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    Intrepid::FunctionSpaceTools::integrate<Real>(*U_new,
                                                  *new_eval,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);

    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> R(3);
    R[0] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    R[1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    R[2] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);

    Intrepid::RealSpaceTools<Real>::add(*R[0], *U_new);
    Intrepid::RealSpaceTools<Real>::subtract(*R[0], *U_old);
    Intrepid::RealSpaceTools<Real>::add(*R[0], *Rold[0]);
    Intrepid::RealSpaceTools<Real>::add(*R[0], *Rnew[0]);

    // APPLY DIRICHLET CONDITIONS
    int numLocalSideIds = bdryCellLocIds_[0].size();
    for (int j = 0; j < numLocalSideIds; ++j)
    {
      int numCellsSide = bdryCellLocIds_[0][j].size();
      int numBdryDofs = fidx_[j].size();
      for (int k = 0; k < numCellsSide; ++k)
      {
        int cidx = bdryCellLocIds_[0][j][k];
        for (int l = 0; l < numBdryDofs; ++l)
        {
          (*R[0])(cidx, fidx_[j][l]) = (*Unew[0])(cidx, fidx_[j][l]) - (*bdryCellDofValues_[0][j])(k, fidx_[j][l]);
        }
      }
    }

    ROL::Ptr<Intrepid::FieldContainer<Real>> v1_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(v1_eval, Unew[1]);
    Intrepid::FunctionSpaceTools::integrate<Real>(*R[1],
                                                  *v1_eval,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);

    ROL::Ptr<Intrepid::FieldContainer<Real>> v2_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(v2_eval, Unew[2]);
    Intrepid::FunctionSpaceTools::integrate<Real>(*R[2],
                                                  *v2_eval,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);

    Intrepid::RealSpaceTools<Real>::add(*R[1], *Rnew[1]);
    Intrepid::RealSpaceTools<Real>::add(*R[2], *Rnew[2]);

    Intrepid::RealSpaceTools<Real>::scale(*R[1], theta_ * dt);
    Intrepid::RealSpaceTools<Real>::scale(*R[2], theta_ * dt);

    // Combine the residuals.
    fieldHelper_->combineFieldCoeff(res, R);
  }

  void Jacobian_uo(ROL::Ptr<Intrepid::FieldContainer<Real>> &jac,
                   const ROL::TimeStamp<Real> &ts,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> &uo_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> &un_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                   const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    const Real one(1);
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    // GET TIME STEP INFORMATION
    Real told = ts.t[0], tnew = ts.t[1], dt = tnew - told;
    // COMPUTE OLD RESIDUAL
    pde_->setTime(told);
    ROL::Ptr<Intrepid::FieldContainer<Real>> j;
    pde_->Jacobian_1(j, uo_coeff, z_coeff, z_param);

    // Initialize jacobians.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> J;
    fieldHelper_->splitFieldCoeff(J, j);
    Intrepid::RealSpaceTools<Real>::scale(*J[0][0], (one - theta_) * dt);

    ROL::Ptr<Intrepid::FieldContainer<Real>> Nold = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);

    // Integrate N * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*Nold,
                                                  *(fe_->N()),
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);

    Intrepid::RealSpaceTools<Real>::subtract(*J[0][0], *Nold);
    Intrepid::RealSpaceTools<Real>::scale(*J[1][0], 0.0);
    Intrepid::RealSpaceTools<Real>::scale(*J[2][0], 0.0);

    // APPLY DIRICHLET CONDITIONS
    int numLocalSideIds = bdryCellLocIds_[0].size();
    for (int j = 0; j < numLocalSideIds; ++j)
    {
      int numCellsSide = bdryCellLocIds_[0][j].size();
      int numBdryDofs = fidx_[j].size();
      for (int k = 0; k < numCellsSide; ++k)
      {
        int cidx = bdryCellLocIds_[0][j][k];
        for (int l = 0; l < numBdryDofs; ++l)
        {
          for (int m = 0; m < f; ++m)
          {
            (*J[0][0])(cidx, fidx_[j][l], m) = static_cast<Real>(0);
          }
          (*J[0][0])(cidx, fidx_[j][l], fidx_[j][l]) = static_cast<Real>(0);
        }
      }
    }

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);
  }

  void Jacobian_un(ROL::Ptr<Intrepid::FieldContainer<Real>> &jac,
                   const ROL::TimeStamp<Real> &ts,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> &uo_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> &un_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                   const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    // Retrieve dimensions.
    // int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    // GET TIME STEP INFORMATION
    Real told = ts.t[0], tnew = ts.t[1], dt = tnew - told;
    // COMPUTE NEW RESIDUAL
    pde_->setTime(tnew);
    ROL::Ptr<Intrepid::FieldContainer<Real>> j;
    pde_->Jacobian_1(j, un_coeff, z_coeff, z_param);

    // Initialize jacobians.
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> J;
    fieldHelper_->splitFieldCoeff(J, j);
    Intrepid::RealSpaceTools<Real>::scale(*J[0][0], theta_ * dt);

    // Integrate N * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][0],
                                                  *(fe_->N()),
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);

    // APPLY DIRICHLET CONDITIONS
    int numLocalSideIds = bdryCellLocIds_[0].size();
    for (int j = 0; j < numLocalSideIds; ++j)
    {
      int numCellsSide = bdryCellLocIds_[0][j].size();
      int numBdryDofs = fidx_[j].size();
      for (int k = 0; k < numCellsSide; ++k)
      {
        int cidx = bdryCellLocIds_[0][j][k];
        for (int l = 0; l < numBdryDofs; ++l)
        {
          for (int m = 0; m < f; ++m)
          {
            (*J[0][0])(cidx, fidx_[j][l], m) = static_cast<Real>(0);
          }
          (*J[0][0])(cidx, fidx_[j][l], fidx_[j][l]) = static_cast<Real>(1);
        }
      }
    }

    Intrepid::FunctionSpaceTools::integrate<Real>(*J[1][1],
                                                  *(fe_->N()),
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);

    Intrepid::FunctionSpaceTools::integrate<Real>(*J[2][2],
                                                  *(fe_->N()),
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, true);

    Intrepid::RealSpaceTools<Real>::scale(*J[1][0], theta_ * dt);
    Intrepid::RealSpaceTools<Real>::scale(*J[1][1], theta_ * dt);
    Intrepid::RealSpaceTools<Real>::scale(*J[2][0], theta_ * dt);
    Intrepid::RealSpaceTools<Real>::scale(*J[2][2], theta_ * dt);

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);
  }

  void Jacobian_zf(ROL::Ptr<Intrepid::FieldContainer<Real>> &jac,
                   const ROL::TimeStamp<Real> &ts,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> &uo_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> &un_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                   const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    const Real one(1);
    // GET DIMENSIONS
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    // GET TIME STEP INFORMATION
    Real told = ts.t[0], tnew = ts.t[1], dt = tnew - told;
    // INITILAIZE JACOBIAN
    ROL::Ptr<Intrepid::FieldContainer<Real>> jo;
    ROL::Ptr<Intrepid::FieldContainer<Real>> jn;
    // ASSEMBLE OLD TIME JACOBIAN
    pde_->setTime(told);
    pde_->Jacobian_2(jo, uo_coeff, z_coeff, z_param);
    // ASSEMBLE NEW TIME JACOBIAN
    pde_->setTime(tnew);
    pde_->Jacobian_2(jn, un_coeff, z_coeff, z_param);

    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> Jo;
    fieldHelper_->splitFieldCoeff(Jo, jo);
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> Jn;
    fieldHelper_->splitFieldCoeff(Jn, jn);

    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> J(3);
    for (int i = 0; i < 3; i++)
    {
      J[i].resize(3);
      for (int j = 0; j < 3; j++)
      {
        J[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
      }
    }

    Intrepid::RealSpaceTools<Real>::scale(*Jo[0][0], (one - theta_) * dt);
    Intrepid::RealSpaceTools<Real>::scale(*Jn[0][0], theta_ * dt);
    Intrepid::RealSpaceTools<Real>::add(*J[0][0], *Jo[0][0]);
    Intrepid::RealSpaceTools<Real>::add(*J[0][0], *Jn[0][0]);

    // APPLY DIRICHLET CONDITIONS
    int numSidesets = bdryCellLocIds_[0].size();
    if (numSidesets > 0)
    {
      int numLocalSideIds = bdryCellLocIds_[0].size();
      for (int j = 0; j < numLocalSideIds; ++j)
      {
        int numCellsSide = bdryCellLocIds_[0][j].size();
        int numBdryDofs = fidx_[j].size();
        for (int k = 0; k < numCellsSide; ++k)
        {
          int cidx = bdryCellLocIds_[0][j][k];
          for (int l = 0; l < numBdryDofs; ++l)
          {
            for (int m = 0; m < f; ++m)
            {
              (*J[0][0])(cidx, fidx_[j][l], m) = static_cast<Real>(0);
            }
          }
        }
      }
    }

    Intrepid::RealSpaceTools<Real>::add(*J[1][0], *Jn[1][0]);
    Intrepid::RealSpaceTools<Real>::add(*J[2][0], *Jn[2][0]);
    Intrepid::RealSpaceTools<Real>::scale(*J[1][0], theta_ * dt);
    Intrepid::RealSpaceTools<Real>::scale(*J[2][0], theta_ * dt);

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);
  }

  void Jacobian_zp(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> &jac,
                   const ROL::TimeStamp<Real> &ts,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> &uo_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> &un_coeff,
                   const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                   const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> (PDE_shallow_ice::Jacobian_zp): Jacobian is zero.");
  }

  void Hessian_uo_uo(ROL::Ptr<Intrepid::FieldContainer<Real>> &hess,
                     const ROL::TimeStamp<Real> &ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> L;
    fieldHelper_->splitFieldCoeff(L, l_coeff);
    // APPLY DIRICHLET CONDITIONS
    int numLocalSideIds = bdryCellLocIds_[0].size();
    for (int j = 0; j < numLocalSideIds; ++j)
    {
      int numCellsSide = bdryCellLocIds_[0][j].size();
      int numBdryDofs = fidx_[j].size();
      for (int k = 0; k < numCellsSide; ++k)
      {
        int cidx = bdryCellLocIds_[0][j][k];
        for (int l = 0; l < numBdryDofs; ++l)
        {
          (*L[0])(cidx, fidx_[j][l]) = 0.0;
        }
      }
    }
    Intrepid::RealSpaceTools<Real>::scale(*L[1], 0.0);
    Intrepid::RealSpaceTools<Real>::scale(*L[2], 0.0);
    ROL::Ptr<Intrepid::FieldContainer<Real>> l0_coeff;
    fieldHelper_->combineFieldCoeff(l0_coeff, L);
    // GET TIME STEP INFORMATION
    Real told = ts.t[0], tnew = ts.t[1], dt = tnew - told;
    // COMPUTE OLD RESIDUAL
    pde_->setTime(told);
    ROL::Ptr<Intrepid::FieldContainer<Real>> ho;
    pde_->Hessian_11(ho, l0_coeff, uo_coeff, z_coeff, z_param);
    std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> H;
    fieldHelper_->splitFieldCoeff(H, ho);
    Intrepid::RealSpaceTools<Real>::scale(*H[0][0], (1.0 - theta_) * dt);
    // Combine the hessians
    fieldHelper_->combineFieldCoeff(hess, H);
  }

  void Hessian_uo_un(ROL::Ptr<Intrepid::FieldContainer<Real>> &hess,
                     const ROL::TimeStamp<Real> &ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> (PDE_shallow_ice::Hessian_uo_un): Hessian is zero.");
  }

  void Hessian_uo_zf(ROL::Ptr<Intrepid::FieldContainer<Real>> &hess,
                     const ROL::TimeStamp<Real> &ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> L;
    fieldHelper_->splitFieldCoeff(L, l_coeff);
    // APPLY DIRICHLET CONDITIONS
    int numLocalSideIds = bdryCellLocIds_[0].size();
    for (int j = 0; j < numLocalSideIds; ++j)
    {
      int numCellsSide = bdryCellLocIds_[0][j].size();
      int numBdryDofs = fidx_[j].size();
      for (int k = 0; k < numCellsSide; ++k)
      {
        int cidx = bdryCellLocIds_[0][j][k];
        for (int l = 0; l < numBdryDofs; ++l)
        {
          (*L[0])(cidx, fidx_[j][l]) = 0.0;
        }
      }
    }
    ROL::Ptr<Intrepid::FieldContainer<Real>> l0_coeff;
    fieldHelper_->combineFieldCoeff(l0_coeff, L);

    // GET TIME STEP INFORMATION
    Real told = ts.t[0], tnew = ts.t[1], dt = tnew - told;
    // COMPUTE OLD RESIDUAL
    pde_->setTime(told);
    pde_->Hessian_12(hess, l0_coeff, uo_coeff, z_coeff, z_param);
    Intrepid::RealSpaceTools<Real>::scale(*hess, (1.0 - theta_) * dt);
  }

  void Hessian_uo_zp(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> &hess,
                     const ROL::TimeStamp<Real> &ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> (PDE_shallow_ice::Hessian_uo_zp): Hessian is zero.");
  }

  void Hessian_un_uo(ROL::Ptr<Intrepid::FieldContainer<Real>> &hess,
                     const ROL::TimeStamp<Real> &ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> (PDE_shallow_ice::Hessian_un_uo): Hessian is zero.");
  }

  void Hessian_un_un(ROL::Ptr<Intrepid::FieldContainer<Real>> &hess,
                     const ROL::TimeStamp<Real> &ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> L;
    fieldHelper_->splitFieldCoeff(L, l_coeff);
    // APPLY DIRICHLET CONDITIONS
    int numLocalSideIds = bdryCellLocIds_[0].size();
    for (int j = 0; j < numLocalSideIds; ++j)
    {
      int numCellsSide = bdryCellLocIds_[0][j].size();
      int numBdryDofs = fidx_[j].size();
      for (int k = 0; k < numCellsSide; ++k)
      {
        int cidx = bdryCellLocIds_[0][j][k];
        for (int l = 0; l < numBdryDofs; ++l)
        {
          (*L[0])(cidx, fidx_[j][l]) = 0.0;
        }
      }
    }
    ROL::Ptr<Intrepid::FieldContainer<Real>> l0_coeff;
    fieldHelper_->combineFieldCoeff(l0_coeff, L);

    // GET TIME STEP INFORMATION
    Real told = ts.t[0], tnew = ts.t[1], dt = tnew - told;
    // COMPUTE NEW RESIDUAL
    pde_->setTime(tnew);
    pde_->Hessian_11(hess, l0_coeff, un_coeff, z_coeff, z_param);
    Intrepid::RealSpaceTools<Real>::scale(*hess, theta_ * dt);
  }

  void Hessian_un_zf(ROL::Ptr<Intrepid::FieldContainer<Real>> &hess,
                     const ROL::TimeStamp<Real> &ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> L;
    fieldHelper_->splitFieldCoeff(L, l_coeff);
    // APPLY DIRICHLET CONDITIONS
    int numLocalSideIds = bdryCellLocIds_[0].size();
    for (int j = 0; j < numLocalSideIds; ++j)
    {
      int numCellsSide = bdryCellLocIds_[0][j].size();
      int numBdryDofs = fidx_[j].size();
      for (int k = 0; k < numCellsSide; ++k)
      {
        int cidx = bdryCellLocIds_[0][j][k];
        for (int l = 0; l < numBdryDofs; ++l)
        {
          (*L[0])(cidx, fidx_[j][l]) = 0.0;
        }
      }
    }
    ROL::Ptr<Intrepid::FieldContainer<Real>> l0_coeff;
    fieldHelper_->combineFieldCoeff(l0_coeff, L);

    // GET TIME STEP INFORMATION
    Real told = ts.t[0], tnew = ts.t[1], dt = tnew - told;
    // COMPUTE OLD RESIDUAL
    pde_->setTime(tnew);
    pde_->Hessian_12(hess, l0_coeff, un_coeff, z_coeff, z_param);
    Intrepid::RealSpaceTools<Real>::scale(*hess, theta_ * dt);
  }

  void Hessian_un_zp(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> &hess,
                     const ROL::TimeStamp<Real> &ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> (PDE_shallow_ice::Hessian_un_zp): Hessian is zero.");
  }

  void Hessian_zf_uo(ROL::Ptr<Intrepid::FieldContainer<Real>> &hess,
                     const ROL::TimeStamp<Real> &ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> L;
    fieldHelper_->splitFieldCoeff(L, l_coeff);
    // APPLY DIRICHLET CONDITIONS
    int numLocalSideIds = bdryCellLocIds_[0].size();
    for (int j = 0; j < numLocalSideIds; ++j)
    {
      int numCellsSide = bdryCellLocIds_[0][j].size();
      int numBdryDofs = fidx_[j].size();
      for (int k = 0; k < numCellsSide; ++k)
      {
        int cidx = bdryCellLocIds_[0][j][k];
        for (int l = 0; l < numBdryDofs; ++l)
        {
          (*L[0])(cidx, fidx_[j][l]) = 0.0;
        }
      }
    }
    ROL::Ptr<Intrepid::FieldContainer<Real>> l0_coeff;
    fieldHelper_->combineFieldCoeff(l0_coeff, L);

    // GET TIME STEP INFORMATION
    Real told = ts.t[0], tnew = ts.t[1], dt = tnew - told;
    // COMPUTE OLD RESIDUAL
    pde_->setTime(told);
    pde_->Hessian_21(hess, l0_coeff, uo_coeff, z_coeff, z_param);
    Intrepid::RealSpaceTools<Real>::scale(*hess, (1.0 - theta_) * dt);
  }

  void Hessian_zf_un(ROL::Ptr<Intrepid::FieldContainer<Real>> &hess,
                     const ROL::TimeStamp<Real> &ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> L;
    fieldHelper_->splitFieldCoeff(L, l_coeff);
    // APPLY DIRICHLET CONDITIONS
    int numLocalSideIds = bdryCellLocIds_[0].size();
    for (int j = 0; j < numLocalSideIds; ++j)
    {
      int numCellsSide = bdryCellLocIds_[0][j].size();
      int numBdryDofs = fidx_[j].size();
      for (int k = 0; k < numCellsSide; ++k)
      {
        int cidx = bdryCellLocIds_[0][j][k];
        for (int l = 0; l < numBdryDofs; ++l)
        {
          (*L[0])(cidx, fidx_[j][l]) = 0.0;
        }
      }
    }
    ROL::Ptr<Intrepid::FieldContainer<Real>> l0_coeff;
    fieldHelper_->combineFieldCoeff(l0_coeff, L);

    // GET TIME STEP INFORMATION
    Real told = ts.t[0], tnew = ts.t[1], dt = tnew - told;
    // COMPUTE OLD RESIDUAL
    pde_->setTime(tnew);
    pde_->Hessian_21(hess, l0_coeff, un_coeff, z_coeff, z_param);
    Intrepid::RealSpaceTools<Real>::scale(*hess, theta_ * dt);
  }

  void Hessian_zf_zf(ROL::Ptr<Intrepid::FieldContainer<Real>> &hess,
                     const ROL::TimeStamp<Real> &ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> L;
    fieldHelper_->splitFieldCoeff(L, l_coeff);
    // APPLY DIRICHLET CONDITIONS
    int numLocalSideIds = bdryCellLocIds_[0].size();
    for (int j = 0; j < numLocalSideIds; ++j)
    {
      int numCellsSide = bdryCellLocIds_[0][j].size();
      int numBdryDofs = fidx_[j].size();
      for (int k = 0; k < numCellsSide; ++k)
      {
        int cidx = bdryCellLocIds_[0][j][k];
        for (int l = 0; l < numBdryDofs; ++l)
        {
          (*L[0])(cidx, fidx_[j][l]) = 0.0;
        }
      }
    }
    ROL::Ptr<Intrepid::FieldContainer<Real>> l0_coeff;
    fieldHelper_->combineFieldCoeff(l0_coeff, L);

    // GET TIME STEP INFORMATION
    Real told = ts.t[0], tnew = ts.t[1], dt = tnew - told;
    // COMPUTE OLD RESIDUAL
    pde_->setTime(told);
    pde_->Hessian_22(hess, l_coeff, uo_coeff, z_coeff, z_param);
    Intrepid::RealSpaceTools<Real>::scale(*hess, (1.0 - theta_) * dt);
    // COMPUTE NEW RESIDUAL
    ROL::Ptr<Intrepid::FieldContainer<Real>> hn;
    pde_->setTime(tnew);
    pde_->Hessian_22(hn, l0_coeff, un_coeff, z_coeff, z_param);
    Intrepid::RealSpaceTools<Real>::scale(*hn, theta_ * dt);
    Intrepid::RealSpaceTools<Real>::add(*hess, *hn);
  }

  void Hessian_zf_zp(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> &hess,
                     const ROL::TimeStamp<Real> &ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> (PDE_shallow_ice::Hessian_zf_zp): Hessian is zero.");
  }

  void Hessian_zp_uo(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> &hess,
                     const ROL::TimeStamp<Real> &ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> (PDE_shallow_ice::Hessian_zp_uo): Hessian is zero.");
  }

  void Hessian_zp_un(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> &hess,
                     const ROL::TimeStamp<Real> &ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> (PDE_shallow_ice::Hessian_zp_un): Hessian is zero.");
  }

  void Hessian_zp_zf(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> &hess,
                     const ROL::TimeStamp<Real> &ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> (PDE_shallow_ice::Hessian_zp_zf): Hessian is zero.");
  }

  void Hessian_zp_zp(std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> &hess,
                     const ROL::TimeStamp<Real> &ts,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &uo_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &un_coeff,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                     const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> (PDE_shallow_ice::Hessian_zp_zp): Hessian is zero.");
  }

  void RieszMap_1(ROL::Ptr<Intrepid::FieldContainer<Real>> &riesz)
  {
    pde_->RieszMap_1(riesz);
  }

  void RieszMap_2(ROL::Ptr<Intrepid::FieldContainer<Real>> &riesz)
  {
    pde_->RieszMap_2(riesz);
  }

  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> getFields()
  {
    return pde_->getFields();
  }

  const ROL::Ptr<PDE<Real>> getPDE(void) const
  {
    return pde_;
  }

  void setParameter(const std::vector<Real> &param)
  {
    DynamicPDE<Real>::setParameter(param);
    pde_->setParameter(param);
  }

  void setFieldPattern(const std::vector<std::vector<int>> &fieldPattern)
  {
    pde_->setFieldPattern(fieldPattern);
    fieldHelper_ = pde_->getFieldHelper();
  }

  void setCellNodes(const ROL::Ptr<Intrepid::FieldContainer<Real>> &volCellNodes,
                    const std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int>>> &bdryCellLocIds)
  {
    bdryCellNodes_ = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    pde_->setCellNodes(volCellNodes, bdryCellNodes, bdryCellLocIds);
    // Finite element definition.
    fe_ = pde_->getFE();
    // Get boundary degrees of freedom.
    fidx_ = fe_->getBoundaryDofs();

    computeDirichlet();
  }

  const ROL::Ptr<FE<Real>> getFE(void) const
  {
    return fe_;
  }

  const ROL::Ptr<FieldHelper<Real>> getFieldHelper(void) const
  {
    return fieldHelper_;
  }

  const std::vector<std::vector<std::vector<int>>> getBdryCellLocIds(void) const
  {
    return bdryCellLocIds_;
  }

private:
  void computeDirichlet(void)
  {
    // Compute Dirichlet values at DOFs.
    // Retrieve dimensions.
    int f = fe_->gradN()->dimension(1);
    int d = fe_->gradN()->dimension(3);
    int numSidesets = bdryCellLocIds_.size();
    bdryCellDofValues_.resize(numSidesets);
    for (int i = 0; i < numSidesets; ++i)
    {
      int numLocSides = bdryCellLocIds_[i].size();
      bdryCellDofValues_[i].resize(numLocSides);
      for (int j = 0; j < numLocSides; ++j)
      {
        int c = bdryCellLocIds_[i][j].size();
        bdryCellDofValues_[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
        ROL::Ptr<Intrepid::FieldContainer<Real>> coords;
        coords = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, d);
        if (c > 0)
        {
          fe_->computeDofCoords(coords, bdryCellNodes_[i][j]);
        }
        for (int k = 0; k < c; ++k)
        {
          for (int l = 0; l < f; ++l)
          {
            std::vector<Real> dofpoint(d);
            for (int m = 0; m < d; ++m)
            {
              dofpoint[m] = (*coords)(k, l, m);
            }
            (*bdryCellDofValues_[i][j])(k, l) = DirichletFunc(dofpoint, i, j);
          }
        }
      }
    }
  }

  Real DirichletFunc(const std::vector<Real> &coords, int sideset, int locSideId) const
  {
    std::vector<int> I = std::vector<int>(N_ + 1);
    Real v = 0;
    if (locSideId == 0) // y=0 boundary
    {
      for (int i = 0; i < N_ + 1; i++)
      {
        I[i] = i;
      }
      v = coords[0] / width_;
    }
    else if (locSideId == 1) // x=1 boundary
    {
      for (int i = 0; i < N_ + 1; i++)
      {
        I[i] = (N_ + 1) * (i + 1) - 1;
      }
      v = coords[1] / height_;
    }
    else if (locSideId == 2) // y=1 boundary
    {
      for (int i = 0; i < N_ + 1; i++)
      {
        I[i] = (N_ + 1) * N_ + i;
      }
      v = coords[0] / width_;
    }
    else if (locSideId == 3) // x=0 boundary
    {
      for (int i = 0; i < N_ + 1; i++)
      {
        I[i] = (N_ + 1) * i;
      }
      v = coords[1] / height_;
    }

    Real val = 0.0;
    for (int i = 0; i < N_ + 1; i++)
    {
      val += surface_height_coeff_[I[i]] * Linear_FE_Basis_Fun_Eval(v, i);
    }
    return val;
  }

  Real Linear_FE_Basis_Fun_Eval(const Real &x, int i) const
  {
    // Evaluates ith 1D FE basis function at point x
    Real val = 0.0;
    Real xi = static_cast<Real>(i) * h_;
    Real dist = std::abs(x - xi);
    if (dist < h_)
    {
      val = 1.0 - dist / h_;
    }
    return val;
  }

}; // DynamicPDE_shallow_ice

#endif
