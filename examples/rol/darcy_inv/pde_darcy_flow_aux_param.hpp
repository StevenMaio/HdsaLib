#ifndef PDE_DARCY_FLOW_AUX_PARAM_HPP
#define PDE_DARCY_FLOW_AUX_PARAM_HPP

#include "../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/pde.hpp"
#include "../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/fe.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "ROL_Ptr.hpp"
#include <cmath>

// Defines the porous media flow PDE
// -\nabla \cDot (\kappa * \nabla u) = z
// where u is the state, z is the source, and \kappa is the permeability

template <class Real>
class PDE_darcy_flow_aux_param : public PDE<Real>
{
private:
  // Finite element basis information
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>> basisPtr_;
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> basisPtrs_;
  // Cell cubature information
  ROL::Ptr<Intrepid::Cubature<Real>> cellCub_;
  // Cell node information
  ROL::Ptr<Intrepid::FieldContainer<Real>> volCellNodes_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> bdryCellNodes_;
  std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;
  // Finite element definition
  ROL::Ptr<FE<Real>> fe_vol_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int>> fidx_;
  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> bdryCellDofValues_;

public:
  PDE_darcy_flow_aux_param(Teuchos::ParameterList &parlist)
  {
    // Finite element fields.
    int basisOrder = parlist.sublist("PDE Poisson").get("Basis Order", 1);
    if (basisOrder == 1)
    {
      basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real>>>();
    }
    else if (basisOrder == 2)
    {
      basisPtr_ = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real>>>();
    }
    basisPtrs_.clear();
    basisPtrs_.push_back(basisPtr_);
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();         // get the cell type from any basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                        // create cubature factory
    int cubDegree = parlist.sublist("PDE Poisson").get("Cubature Degree", 2); // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create(cellType, cubDegree);                        // create default cubature
  }

  void residual(ROL::Ptr<Intrepid::FieldContainer<Real>> &res,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    // GET DIMENSIONS
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    int p = fe_vol_->gradN()->dimension(2);
    // int d = fe_vol_->gradN()->dimension(3);
    //  INITIALIZE RESIDUAL
    res = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // COMPUTE PDE COEFFICIENTS
    ROL::Ptr<Intrepid::FieldContainer<Real>> rhs = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);

    // Add source term to residual
    computeSource(rhs, *z_param);
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *rhs,
                                                  *(fe_vol_->NdetJ()),
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
          (*res)(cidx, fidx_[j][l]) = (*u_coeff)(cidx, fidx_[j][l]) - (*bdryCellDofValues_[0][j])(k, fidx_[j][l]);
        }
      }
    }
  }

  void Jacobian_1(ROL::Ptr<Intrepid::FieldContainer<Real>> &jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    // GET DIMENSIONS
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    // int p = fe_vol_->gradN()->dimension(2);
    // int d = fe_vol_->gradN()->dimension(3);
    //  INITILAIZE JACOBIAN
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);

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
          // std::cout << "\n   j=" << j << "  l=" << l << "  " << fidx[j][l];
          for (int m = 0; m < f; ++m)
          {
            (*jac)(cidx, fidx_[j][l], m) = static_cast<Real>(0);
          }
          (*jac)(cidx, fidx_[j][l], fidx_[j][l]) = static_cast<Real>(1);
        }
      }
    }
  }

  void Jacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real>> &jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> (PDE_darcy_flow::Jacobian_2): Jacobian is zero.");
  }

  void Jacobian_3(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> &jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    // GET DIMENSIONS
    int c = fe_vol_->gradN()->dimension(0);
    int f = fe_vol_->gradN()->dimension(1);
    int p = fe_vol_->gradN()->dimension(2);
    // int d = fe_vol_->gradN()->dimension(3);
    int size = z_param->size();
    for (int i = 0; i < size; i++)
    {
      jac[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);

      ROL::Ptr<Intrepid::FieldContainer<Real>> rhs = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);

      // Add source term to residual
      computeSource(rhs, *z_param, true, i);
      Intrepid::FunctionSpaceTools::integrate<Real>(*jac[i],
                                                    *rhs,
                                                    *(fe_vol_->NdetJ()),
                                                    Intrepid::COMP_CPP, false);

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
            (*jac[i])(cidx, fidx_[j][l]) = 0.0;
          }
        }
      }
    }
  }

  void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real>> &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> (PDE_darcy_flow::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real>> &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> (PDE_darcy_flow::Hessian_12): Hessian is zero.");
  }

  void Hessian_13(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> (PDE_darcy_flow::Hessian_13): Hessian is zero.");
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real>> &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> (PDE_darcy_flow::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real>> &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> (PDE_darcy_flow::Hessian_22): Hessian is zero.");
  }

  void Hessian_23(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> (PDE_darcy_flow::Hessian_23): Hessian is zero.");
  }

  void Hessian_31(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> (PDE_darcy_flow::Hessian_31): Hessian is zero.");
  }

  void Hessian_32(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> (PDE_darcy_flow::Hessian_32): Hessian is zero.");
  }

  void Hessian_33(std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> (PDE_darcy_flow::Hessian_33): Hessian is zero.");
  }

  void RieszMap_1(ROL::Ptr<Intrepid::FieldContainer<Real>> &riesz)
  {
    // GET DIMENSIONS
    int c = fe_vol_->N()->dimension(0);
    int f = fe_vol_->N()->dimension(1);
    // INITIALIZE RIESZ
    riesz = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    *riesz = *fe_vol_->stiffMat();
    Intrepid::RealSpaceTools<Real>::add(*riesz, *(fe_vol_->massMat()));
  }

  void RieszMap_2(ROL::Ptr<Intrepid::FieldContainer<Real>> &riesz)
  {
    // GET DIMENSIONS
    int c = fe_vol_->N()->dimension(0);
    int f = fe_vol_->N()->dimension(1);
    // INITIALIZE RIESZ
    riesz = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    *riesz = *fe_vol_->massMat();
  }

  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> getFields()
  {
    return basisPtrs_;
  }

  void setCellNodes(const ROL::Ptr<Intrepid::FieldContainer<Real>> &volCellNodes,
                    const std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int>>> &bdryCellLocIds)
  {
    volCellNodes_ = volCellNodes;
    bdryCellNodes_ = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    // Finite element definition.
    fe_vol_ = ROL::makePtr<FE<Real>>(volCellNodes_, basisPtr_, cellCub_);
    // set local boundary DOFs.
    fidx_ = fe_vol_->getBoundaryDofs();
    // Compute Dirichlet values at DOFs.
    int d = basisPtr_->getBaseCellTopology().getDimension();
    int numSidesets = bdryCellLocIds_.size();
    bdryCellDofValues_.resize(numSidesets);
    for (int i = 0; i < numSidesets; ++i)
    {
      int numLocSides = bdryCellLocIds_[i].size();
      bdryCellDofValues_[i].resize(numLocSides);
      for (int j = 0; j < numLocSides; ++j)
      {
        int c = bdryCellLocIds_[i][j].size();
        int f = basisPtr_->getCardinality();
        bdryCellDofValues_[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
        ROL::Ptr<Intrepid::FieldContainer<Real>> coords =
            ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, d);
        if (c > 0)
        {
          fe_vol_->computeDofCoords(coords, bdryCellNodes_[i][j]);
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
            (*bdryCellDofValues_[i][j])(k, l) = evaluateDirichlet(dofpoint, i, j);
          }
        }
      }
    }
  }

  const ROL::Ptr<FE<Real>> getFE(void) const
  {
    return fe_vol_;
  }

private:
private:
  void computeSource(ROL::Ptr<Intrepid::FieldContainer<Real>> &rhs, const std::vector<Real> &param, const bool is_derv = false, const int &diff_coord = 0) const
  {
    int c = fe_vol_->gradN()->dimension(0);
    int p = fe_vol_->gradN()->dimension(2);
    int d = fe_vol_->gradN()->dimension(3);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i)
    {
      for (int j = 0; j < p; ++j)
      {
        for (int k = 0; k < d; ++k)
        {
          pt[k] = (*fe_vol_->cubPts())(i, j, k);
        }
        // Compute forcing term f
        (*rhs)(i, j) = -evaluateRHS(pt, param, is_derv, diff_coord);
      }
    }
  }

  Real evaluateRHS(const std::vector<Real> &x, const std::vector<Real> &param, const bool is_derv, const int &diff_coord) const
  {

    Real val = 0.0;
    int n = std::sqrt(param.size());
    if (is_derv)
    {
      int j = diff_coord % n;
      int i = (diff_coord - j) / n;
      val = 100.0 * sin(static_cast<Real>(2 * (i + 1)) * 3.14159 * x[0]) * sin(static_cast<Real>(2 * (j + 1)) * 3.14159 * x[1]) / static_cast<Real>((i + 1) * (j + 1));
    }
    else
    {
      for (int i = 0; i < n; i++)
      {
        for (int j = 0; j < n; j++)
        {
          int k = i * n + j;
          val += param[k] * 100.0 * sin(static_cast<Real>(2 * (i + 1)) * 3.14159 * x[0]) * sin(static_cast<Real>(2 * (j + 1)) * 3.14159 * x[1]) / static_cast<Real>((i + 1) * (j + 1));
        }
      }
    }
    return val;
  }

  Real evaluateDirichlet(const std::vector<Real> &coords, int sideset, int locSideId) const
  {
    Real val = 0.0;
    if (coords[1] == 1.0)
    {
      Real x = 2 * 3.14159 * coords[0];
      val = std::pow(sin(x), 2.0);
    }
    else
    {
      Real x = 2 * 3.14159 * coords[0];
      val = std::pow(cos(2 * x), 2.0);
    }
    return val;
  }

}; // PDE_darcy_flow

#endif
