#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_Stream.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_Objective_SimOpt.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "ROL_ReducedDynamicObjective.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_DynamicConstraintCheck.hpp"
#include "ROL_DynamicObjectiveCheck.hpp"

#include "../../../PDE-OPT/TOOLS/solver.cpp"
#include "../../../PDE-OPT/TOOLS/solver_def.hpp"
#include "../../../PDE-OPT/TOOLS/assembler.cpp"
#include "../../../PDE-OPT/TOOLS/assembler_def.hpp"
#include "../../../PDE-OPT/TOOLS/meshmanager.hpp"
#include "../../../PDE-OPT/TOOLS/dynconstraint.hpp"
#include "../../../PDE-OPT/TOOLS/ltiobjective.hpp"
#include "../../../PDE-OPT/TOOLS/pdevector.hpp"
#include "../../../PDE-OPT/TOOLS/pdeobjective.hpp"

#include "../../adapters/rol/modified_rol_source_code/ROL_ReducedDynamicObjective_Stationary_Control.hpp"
#include "../../adapters/rol/modified_rol_source_code/Objective_SimOpt_TS.hpp"
#include "../../adapters/rol/modified_rol_source_code/Misfit_Regularization_TS_Objective_SimOpt.hpp"
#include "../../adapters/rol/modified_rol_source_code/ltiobjective_TS.hpp"
#include "../../adapters/rol/modified_rol_source_code/Reduced_Objective_Regularization.hpp"

#include "mesh_shallow_ice.hpp"
#include "pde_shallow_ice.hpp"
#include "dynpde_shallow_ice.hpp"
#include "obj_shallow_ice.hpp"

#include "../../src/source_file.hpp"
#include "HDSA_ROL_B_Transpose_shallow_ice.hpp"
#include "Opt_Problem_Objects_shallow_ice.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {

  HDSA::nullstream bhs;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  HDSA::Ptr<const HDSA::Comm<int> > comm = HDSA::makePtr<HDSA::Comm<int> >();
  std::string filename = "input.xml";
  HDSA::Ptr<HDSA::ParameterList> parlist = HDSA::makePtr<HDSA::ParameterList>();
  HDSA::updateParametersFromXmlFile( filename, *parlist );

  int L = parlist->sublist("Problem").get("Number of Uncertain Basis Functions", 10);
  int dim = (L+1)*(L+1);
  HDSA::Ptr<HDSA::Vector<RealT> > theta = HDSA::makePtr<Std_Vector<RealT> >(dim);

  HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects = HDSA::makePtr<Opt_Problem_Objects_shallow_ice<RealT> >(parlist,theta,comm);
  OP_Objects->Solve_Optimization_Problem();
  OP_Objects->Write_Optimal_Solution();

  return 0;
}
