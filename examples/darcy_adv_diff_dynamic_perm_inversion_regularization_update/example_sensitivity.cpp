#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_Stream.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_Objective_SimOpt.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "ROL_ReducedDynamicObjective.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_DynamicConstraintCheck.hpp"
#include "ROL_DynamicObjectiveCheck.hpp"

//#include <iostream>
//#include <random>
//#include <unordered_set>
//#include <set>

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

#include "mesh_darcy_adv_diff.hpp"
#include "pde_darcy_adv_diff.hpp"
#include "dynpde_darcy_adv_diff.hpp"
#include "obj_darcy_adv_diff.hpp"

#include "../../src/source_file.hpp"
#include "Weight_Matrices_darcy_adv_diff.hpp"
#include "Opt_Problem_Objects_darcy_adv_diff.hpp"
#include "Parameter_Sampler_darcy_adv_diff.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {

  HDSA::nullstream bhs;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  HDSA::Ptr<const HDSA::Comm<int> > comm = HDSA::makePtr<HDSA::Comm<int> >();
  /*** Read in XML input ***/
  std::string filename = "input.xml";
  HDSA::Ptr<HDSA::ParameterList> parlist = HDSA::makePtr<HDSA::ParameterList>();
  HDSA::updateParametersFromXmlFile( filename, *parlist );
  std::string filenameSensitivity = "Sensitivity_input.xml";
  HDSA::Ptr<HDSA::ParameterList> parlist_sensitivity = HDSA::makePtr<HDSA::ParameterList>();
  HDSA::updateParametersFromXmlFile( filenameSensitivity, *parlist_sensitivity );
  
  int L = parlist->sublist("Problem").get("Number of Uncertainty Basis Function", 10);
  int dim = 9*16 + 1 + (L+1) + (L+1);

  HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects_Factory = HDSA::makePtr<Opt_Problem_Objects_darcy_adv_diff<RealT> >(parlist, comm);
  HDSA::Ptr<HDSA::Weight_Matrices<RealT> > weight_matrices_factory = HDSA::makePtr<Weight_Matrices_darcy_adv_diff<RealT> >(parlist, parlist_sensitivity);
  HDSA::Ptr<HDSA::Parameter_Sampler<RealT> > sampler = HDSA::makePtr<Parameter_Sampler_darcy_adv_diff<RealT> >(dim);
  HDSA::Sample_Local_Sensitivities<RealT>(comm,parlist_sensitivity,OP_Objects_Factory,weight_matrices_factory,sampler);  

  return 0;
}
