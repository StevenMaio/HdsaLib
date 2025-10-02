#ifndef ELLIPTIC_PRIOR_REGULARIZATION_OBJECTIVE_HPP
#define ELLIPTIC_PRIOR_REGULARIZATION_OBJECTIVE_HPP

#include "ROL_Objective_SimOpt.hpp"
#include "../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/linearpdeconstraint.hpp"

template <class Real>
class PDE_Reg_Op : public PDE<Real>
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

  Real gamma_, alpha_;
  bool dirichlet_;

public:
  PDE_Reg_Op(Teuchos::ParameterList &parlist)
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

    gamma_ = parlist.sublist("Problem").get("gamma", 1.0);
    alpha_ = parlist.sublist("Problem").get("alpha", 1.0);
    dirichlet_ = parlist.sublist("Problem").get("Dirichlet Prior", true);
  }

  void residual(ROL::Ptr<Intrepid::FieldContainer<Real>> &res,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    // GET DIMENSIONS
    int c = u_coeff->dimension(0);
    int p = cellCub_->getNumPoints();
    int f = basisPtr_->getCardinality();
    int d = cellCub_->getDimension();

    // Initialize residuals.
    res = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);

    // COMPUTE STIFFNESS TERM
    ROL::Ptr<Intrepid::FieldContainer<Real>> U_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_vol_->evaluateValue(U_eval, u_coeff);
    ROL::Ptr<Intrepid::FieldContainer<Real>> gradU_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);

    ROL::Ptr<Intrepid::FieldContainer<Real>> laplace_term = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    Intrepid::FunctionSpaceTools::integrate<Real>(*laplace_term, *gradU_eval, *(fe_vol_->gradNdetJ()), Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*laplace_term, gamma_);

    // ADD REACTION TERM TO RESIDUAL
    ROL::Ptr<Intrepid::FieldContainer<Real>> eye_term = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // Integrate U * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*eye_term,
                                                  *U_eval,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*eye_term, alpha_);

    Intrepid::RealSpaceTools<Real>::scale(*res, 0.0);
    Intrepid::RealSpaceTools<Real>::add(*res, *laplace_term);
    Intrepid::RealSpaceTools<Real>::add(*res, *eye_term);
  }

  void Jacobian_1(ROL::Ptr<Intrepid::FieldContainer<Real>> &jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    // GET DIMENSIONS
    int c = u_coeff->dimension(0);
    // int p = cellCub_->getNumPoints();
    int f = basisPtr_->getCardinality();
    // int d = cellCub_->getDimension();

    // INITILAIZE JACOBIAN
    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);

    // COMPUTE STIFFNESS TERM
    ROL::Ptr<Intrepid::FieldContainer<Real>> laplace_term = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    Intrepid::FunctionSpaceTools::integrate<Real>(*laplace_term, *(fe_vol_->gradN()), *(fe_vol_->gradNdetJ()), Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*laplace_term, gamma_);

    // ADD REACTION TERM TO JACOBIAN
    ROL::Ptr<Intrepid::FieldContainer<Real>> eye_term = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    // Integrate N * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*eye_term,
                                                  *(fe_vol_->N()),
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*eye_term, alpha_);

    Intrepid::RealSpaceTools<Real>::scale(*jac, 0.0);
    Intrepid::RealSpaceTools<Real>::add(*jac, *laplace_term);
    Intrepid::RealSpaceTools<Real>::add(*jac, *eye_term);
  }

  void Jacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real>> &jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    // GET DIMENSIONS
    int c = u_coeff->dimension(0);
    // int p = cellCub_->getNumPoints();
    int f = basisPtr_->getCardinality();
    // int d = cellCub_->getDimension();

    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
  }

  void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real>> &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> (PDE_Reg_Op::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real>> &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> (PDE_Reg_Op::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real>> &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> (PDE_Reg_Op::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real>> &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> (PDE_Reg_Op::Hessian_22): Hessian is zero.");
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
            (*bdryCellDofValues_[i][j])(k, l) = 0.0;
          }
        }
      }
    }
  }

  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> getFields()
  {
    return basisPtrs_;
  }

  const ROL::Ptr<FE<Real>> getFE(void) const
  {
    return fe_vol_;
  }

}; // PDE_Reg_Op

template <class Real>
class PDE_Mass_Mat : public PDE<Real>
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
  PDE_Mass_Mat(Teuchos::ParameterList &parlist)
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
    int c = u_coeff->dimension(0);
    int p = cellCub_->getNumPoints();
    int f = basisPtr_->getCardinality();
    // int d = cellCub_->getDimension();

    // Initialize residuals.
    res = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);

    // COMPUTE STIFFNESS TERM
    ROL::Ptr<Intrepid::FieldContainer<Real>> U0_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_vol_->evaluateValue(U0_eval, u_coeff);
    // Integrate U * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*res,
                                                  *U0_eval,
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

  void Jacobian_1(ROL::Ptr<Intrepid::FieldContainer<Real>> &jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    // GET DIMENSIONS
    int c = u_coeff->dimension(0);
    // int p = cellCub_->getNumPoints();
    int f = basisPtr_->getCardinality();
    // int d = cellCub_->getDimension();

    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);

    // Integrate N * N
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac,
                                                  *(fe_vol_->N()),
                                                  *(fe_vol_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

  void Jacobian_2(ROL::Ptr<Intrepid::FieldContainer<Real>> &jac,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    // GET DIMENSIONS
    int c = u_coeff->dimension(0);
    // int p = cellCub_->getNumPoints();
    int f = basisPtr_->getCardinality();
    // int d = cellCub_->getDimension();

    jac = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
  }

  void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real>> &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> (PDE_Reg_Op::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(ROL::Ptr<Intrepid::FieldContainer<Real>> &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> (PDE_Reg_Op::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(ROL::Ptr<Intrepid::FieldContainer<Real>> &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> (PDE_Reg_Op::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(ROL::Ptr<Intrepid::FieldContainer<Real>> &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &l_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr)
  {
    throw Exception::Zero(">>> (PDE_Reg_Op::Hessian_22): Hessian is zero.");
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
            (*bdryCellDofValues_[i][j])(k, l) = 0.0;
          }
        }
      }
    }
  }

  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> getFields()
  {
    return basisPtrs_;
  }

  const ROL::Ptr<FE<Real>> getFE(void) const
  {
    return fe_vol_;
  }

}; // PDE_Mass_Mat

template <class Real>
class Elliptic_Prior_Regularization_Objective : public ROL::Objective<Real>
{

private:
  ROL::Ptr<ROL::Constraint_SimOpt<Real>> con_;
  ROL::Ptr<ROL::Constraint_SimOpt<Real>> mass_mat_con_;
  ROL::Ptr<ROL::Vector<Real>> prior_mean_;

public:
  Elliptic_Prior_Regularization_Objective(const ROL::Ptr<const Teuchos::Comm<int>> &comm, const ROL::Ptr<Teuchos::ParameterList> &parlist, ROL::Ptr<std::ostream> &outStream, bool construct_operators = false)
      : ROL::Objective<Real>()
  {
    ROL::Ptr<MeshManager<Real>> meshMgr = ROL::makePtr<MeshManager_Rectangle<Real>>(*parlist);

    ROL::Ptr<PDE<Real>> elliptic_pde = ROL::makePtr<PDE_Reg_Op<Real>>(*parlist);
    con_ = ROL::makePtr<Linear_PDE_Constraint<Real>>(elliptic_pde, meshMgr, comm, *parlist, *outStream, true);

    ROL::Ptr<PDE<Real>> mass_mat_pde = ROL::makePtr<PDE_Mass_Mat<Real>>(*parlist);
    mass_mat_con_ = ROL::makePtr<Linear_PDE_Constraint<Real>>(mass_mat_pde, meshMgr, comm, *parlist, *outStream, true);
  }

  void Elliptic_Solve(ROL::Vector<Real> &z_out, const ROL::Vector<Real> &z_in, Real &tol)
  {
    con_->applyInverseJacobian_1(z_out, z_in, z_in, z_in, tol);
  }

  void Adjoint_Elliptic_Solve(ROL::Vector<Real> &z_out, const ROL::Vector<Real> &z_in, Real &tol)
  {
    con_->applyInverseAdjointJacobian_1(z_out, z_in, z_in, z_in, tol);
  }

  void Apply_Mass_Matrix(ROL::Vector<Real> &z_out, const ROL::Vector<Real> &z_in, Real &tol)
  {
    mass_mat_con_->applyJacobian_1(z_out, z_in, z_in, z_in, tol);
  }

  void Apply_Prior_Covariance(ROL::Vector<Real> &z_out, const ROL::Vector<Real> &z_in, Real &tol)
  {
    ROL::Ptr<ROL::Vector<Real>> z_tmp1 = z_out.clone();
    ROL::Ptr<ROL::Vector<Real>> z_tmp2 = z_out.clone();
    con_->applyInverseAdjointJacobian_1(*z_tmp1, z_in, z_in, z_in, tol);
    mass_mat_con_->applyJacobian_1(*z_tmp2, *z_tmp1, z_in, z_in, tol);
    con_->applyInverseJacobian_1(z_out, *z_tmp2, z_in, z_in, tol);
  }

  Real value(const ROL::Vector<Real> &z, Real &tol) override
  {
    ROL::Ptr<ROL::Vector<Real>> Az = z.clone();
    con_->applyJacobian_1(*Az, z, z, z, tol);
    ROL::Ptr<ROL::Vector<Real>> MinvAz = z.clone();
    mass_mat_con_->applyInverseJacobian_1(*MinvAz, *Az, z, z, tol);
    ROL::Ptr<ROL::Vector<Real>> AtMinvAz = z.clone();
    con_->applyAdjointJacobian_1(*AtMinvAz, *MinvAz, z, z, tol);
    Real val = 0.5 * AtMinvAz->dot(z);
    return val;
  }

  void gradient(ROL::Vector<Real> &g, const ROL::Vector<Real> &z, Real &tol) override
  {
    ROL::Ptr<ROL::Vector<Real>> Az = z.clone();
    con_->applyJacobian_1(*Az, z, z, z, tol);
    ROL::Ptr<ROL::Vector<Real>> MinvAz = z.clone();
    mass_mat_con_->applyInverseJacobian_1(*MinvAz, *Az, z, z, tol);
    con_->applyAdjointJacobian_1(g, *MinvAz, z, z, tol);
  }

  void hessVec(ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &z, Real &tol) override
  {
    ROL::Ptr<ROL::Vector<Real>> Av = z.clone();
    con_->applyJacobian_1(*Av, v, z, z, tol);
    ROL::Ptr<ROL::Vector<Real>> MinvAv = z.clone();
    mass_mat_con_->applyInverseJacobian_1(*MinvAv, *Av, z, z, tol);
    con_->applyAdjointJacobian_1(hv, *MinvAv, z, z, tol);
  }

  void precond(ROL::Vector<Real> &Pv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &x, Real &tol) override
  {
    Apply_Prior_Covariance(Pv, v, tol);
  }

}; // Elliptic_Prior_Regularization_Objective

#endif
