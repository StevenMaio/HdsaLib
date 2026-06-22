/***********************************************************************
 HdsaLib - A library for Hyper-differential Sensitivity Analysis

 Questions? Contact Joseph Hart (joshart@sandia.gov)
************************************************************************/

#ifndef HDSA_BF_SOL_OP_INTERFACE_MRHYDE_HPP
#define HDSA_BF_SOL_OP_INTERFACE_MRHYDE_HPP

#include "userInterface.hpp"
#include "HDSA_BF_Sol_Op_Interface.hpp"
#include "HDSA_Tpetra_Vector.hpp"
#include "HDSA_Std_Vector.hpp"
#include "HDSA_Solver_Interface_MrHyDE.hpp"

template <class RealT>
class BF_Sol_Op_Interface_MrHyDE : public HDSA::BF_Sol_Op_Interface<RealT>
{

private:
  HDSA::Ptr<MrHyDE::SolverManager<SolverNode>> solver_;
  HDSA::Ptr<MrHyDE::PostprocessManager<SolverNode>> postproc_;
  HDSA::Ptr<MrHyDE::ParameterManager<SolverNode>> params_;
  HDSA::Ptr<Solver_Interface_MrHyDE<RealT>> solver_interface_;

public:
  BF_Sol_Op_Interface_MrHyDE(Teuchos::RCP<Teuchos::MpiComm<int>> &comm)
  {
    std::string input_file_name = "input_hifi.yaml";
    Teuchos::RCP<MrHyDE::userInterface> UI = Teuchos::rcp(new MrHyDE::userInterface());
    Teuchos::RCP<Teuchos::ParameterList> Settings = UI->UserInterface(input_file_name);
    Instantiate_HiFi_Model(comm, Settings);

    solver_interface_ = HDSA::makePtr<Solver_Interface_MrHyDE<RealT>>(solver_, params_);

    postproc_->hdsa_solop_data.resize(solver_->setnames.size());
    for (int set = 0; set < solver_->setnames.size(); set++)
    {
      postproc_->hdsa_solop_data[set] = HDSA::makePtr<MrHyDE::SolutionStorage<SolverNode>>(solver_->settings);
    }
  }

  virtual ~BF_Sol_Op_Interface_MrHyDE()
  {
  }

  void Instantiate_HiFi_Model(Teuchos::RCP<Teuchos::MpiComm<int>> comm, Teuchos::RCP<Teuchos::ParameterList> &Settings)
  {

    Teuchos::RCP<MrHyDE::MeshInterface> mesh = Teuchos::rcp(new MrHyDE::MeshInterface(Settings, comm));

    Teuchos::RCP<MrHyDE::PhysicsInterface> physics = Teuchos::rcp(new MrHyDE::PhysicsInterface(Settings, comm,
                                                                                               mesh->getBlockNames(),
                                                                                               mesh->getSideNames(),
                                                                                               mesh->getDimension()));

    mesh->finalize(physics->getVarList(), physics->getVarTypes(), physics->getDerivedList());

    Teuchos::RCP<MrHyDE::DiscretizationInterface> disc = Teuchos::rcp(new MrHyDE::DiscretizationInterface(Settings, comm,
                                                                                                          mesh, physics));
    params_ = Teuchos::rcp(new MrHyDE::ParameterManager<SolverNode>(comm, Settings, mesh, physics, disc));

    Teuchos::RCP<MrHyDE::AssemblyManager<SolverNode>> assembler = Teuchos::rcp(new MrHyDE::AssemblyManager<SolverNode>(comm, Settings, mesh, disc, physics, params_));

    assembler->setMeshData();

    Teuchos::RCP<MrHyDE::MultiscaleManager> multiscale_manager = Teuchos::rcp(new MrHyDE::MultiscaleManager(comm, mesh, Settings,
                                                                                                            assembler->groups,
#ifndef MrHyDE_NO_AD
                                                                                                            assembler->function_managers_AD));
#else
                                                                                                            assembler->function_managers));
#endif

    postproc_ = Teuchos::rcp(new MrHyDE::PostprocessManager<SolverNode>(comm, Settings, mesh,
                                                                        disc, physics,
                                                                        multiscale_manager,
                                                                        assembler, params_));

    solver_ = Teuchos::rcp(new MrHyDE::SolverManager<SolverNode>(comm, Settings, mesh, disc, physics, assembler, params_));

    solver_->multiscale_manager = multiscale_manager;
    assembler->multiscale_manager = multiscale_manager;
    solver_->postproc = postproc_;

    mesh->allocateMeshDataStructures();
    assembler->allocateGroupStorage();
    solver_->completeSetup();
    postproc_->linalg = solver_->linalg;
    solver_->setupExplicitMass();
    assembler->finalizeFunctions();
    solver_->finalizeMultiscale();

    Kokkos::fence();
    comm->barrier();
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Implementation of base class pure virtual functions:
  // State_Solve
  // Apply_Solution_Operator_z_Jacobian_Transpose
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void State_Solve(HDSA::Vector<RealT> &u, const HDSA::Vector<RealT> &z) const
  {
    solver_interface_->State_Solve(u, z);
  }

  void Apply_Solution_Operator_z_Jacobian_Transpose(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &u_in, const HDSA::Vector<RealT> &z) const
  {
    Write_Data_Solution_Operator(u_in);
    Do_Solution_Operator(true);
    RS_Gradient(z_out, z);
    z_out.Scale(-1.0);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Functions that aid in the implementation of the base class pure virtual function
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void Write_Data_Solution_Operator(const HDSA::Vector<RealT> &u) const
  {
    if (solver_->isTransient)
    {
      const HDSA::Transient_Vector<RealT> &eu = dynamic_cast<const HDSA::Transient_Vector<RealT> &>(u);
      int n_t = solver_->settings->sublist("Solver").get<int>("number of steps", 0) + 1;
      for (int i = 0; i < n_t; i++)
      {
        const HDSA::Tpetra_Vector<RealT> &eu_i = dynamic_cast<const HDSA::Tpetra_Vector<RealT> &>(*eu[i]);
        RealT currenttime = solver_->initial_time + (double)i * solver_->deltat;
        postproc_->hdsa_solop_data[0]->store(eu_i.getVector(), currenttime, 0);
      }
    }
    else
    {
      const HDSA::Tpetra_Vector<RealT> &eu = dynamic_cast<const HDSA::Tpetra_Vector<RealT> &>(u);
      postproc_->hdsa_solop_data[0]->store(eu.getVector(), 0.0, 0);
    }
  }

  void Do_Solution_Operator(bool solop_flag) const
  {
    postproc_->hdsa_solop = solop_flag;
  }

  void RS_Gradient(HDSA::Vector<RealT> &grad_z, const HDSA::Vector<RealT> &z) const
  {
    bool new_z = Check_New_Params(z);
    if (new_z)
    {
      HDSA::Ptr<MrHyDE_OptVector> z_rol = solver_interface_->Map_HDSA_Vector_to_MrHyDE_OptVector(z);
      MrHyDE_OptVector curr_z = params_->getCurrentVector();
      ROL::Ptr<ROL::Vector<RealT>> z_tmp = curr_z.clone();
      MrHyDE_OptVector ez_tmp = Teuchos::dyn_cast<MrHyDE_OptVector>(dynamic_cast<ROL::Vector<RealT> &>(*z_tmp));
      ez_tmp.set(*z_rol);

      params_->updateParams(ez_tmp);
      ScalarT val = 0.0;
      solver_->forwardModel(val);
    }

    HDSA::Ptr<MrHyDE_OptVector> grad_z_rol = solver_interface_->Map_HDSA_Vector_to_MrHyDE_OptVector(grad_z);
    grad_z_rol->zero();

    solver_->adjointModel(*grad_z_rol);
  }

  bool Check_New_Params(const HDSA::Vector<RealT> &z) const
  {
    HDSA::Ptr<MrHyDE_OptVector> z_rol = solver_interface_->Map_HDSA_Vector_to_MrHyDE_OptVector(z);
    MrHyDE_OptVector curr_z = params_->getCurrentVector();

    ROL::Ptr<ROL::Vector<RealT>> diff = curr_z.clone();
    diff->zero();
    diff->set(curr_z);
    diff->axpy(-1.0, *z_rol);

    ScalarT dnorm = diff->norm();
    ScalarT refnorm = curr_z.norm();
    dnorm = dnorm / refnorm;
    ScalarT reltol = 1.0e-12;

    bool new_z = false;
    if (dnorm > reltol)
    {
      new_z = true;
    }
    return new_z;
  }
};
#endif
