#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_OptimizationSolver.hpp"

#include "../../../../PDE-OPT/TOOLS/solver.cpp"
#include "../../../../PDE-OPT/TOOLS/solver_def.hpp"
#include "../../../../PDE-OPT/TOOLS/assembler.cpp"
#include "../../../../PDE-OPT/TOOLS/assembler_def.hpp"
#include "../../../../PDE-OPT/TOOLS/meshmanager.hpp"
#include "../../../../PDE-OPT/TOOLS/linearpdeconstraint.hpp"
#include "../../../../PDE-OPT/TOOLS/pdeobjective.hpp"
#include "../../../../PDE-OPT/TOOLS/pdevector.hpp"
#include "../pde_darcy_flow.hpp"
#include "../obj_darcy_flow.hpp"
#include "../mesh_darcy_flow.hpp"

#include "../../../src/source_file.hpp"
#include "Opt_Problem_Objects_darcy_flow_deterministic.hpp"

#include "Weight_Matrices_darcy_flow_deterministic.hpp"
#include "Parameter_Sampler_darcy_flow_deterministic.hpp"

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
  
  HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects_Factory = HDSA::makePtr<Opt_Problem_Objects_darcy_flow_deterministic<RealT> >(parlist, comm);
  HDSA::Ptr<HDSA::Weight_Matrices<RealT> > weight_matrices_factory = HDSA::makePtr<Weight_Matrices_darcy_flow_deterministic<RealT> >(parlist_sensitivity);
  HDSA::Ptr<HDSA::Parameter_Sampler<RealT> > sampler = HDSA::makePtr<Parameter_Sampler_darcy_flow_deterministic<RealT> >(parlist);
  HDSA::Sample_Local_Sensitivities<RealT>(comm,parlist_sensitivity,OP_Objects_Factory,weight_matrices_factory,sampler);  

  return 0;
}
