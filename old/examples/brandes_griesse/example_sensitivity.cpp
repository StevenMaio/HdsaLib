#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_OptimizationSolver.hpp"

#include "../../../PDE-OPT/TOOLS/solver.cpp"
#include "../../../PDE-OPT/TOOLS/solver_def.hpp"
#include "../../../PDE-OPT/TOOLS/assembler.cpp"
#include "../../../PDE-OPT/TOOLS/assembler_def.hpp"
#include "../../../PDE-OPT/TOOLS/meshmanager.hpp"
#include "../../../PDE-OPT/TOOLS/linearpdeconstraint.hpp"
#include "../../../PDE-OPT/TOOLS/pdevector.hpp"
#include "../../../PDE-OPT/TOOLS/pdeobjective.hpp"
#include "pde_brandes_griesse.hpp"
#include "obj_brandes_griesse.hpp"

#include "../../src/source_file.hpp"
#include "Weight_Matrices_brandes_griesse.hpp"
#include "Opt_Problem_Objects_brandes_griesse.hpp"
#include "Parameter_Sampler_brandes_griesse.hpp"

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
  
  HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects_Factory = HDSA::makePtr<Opt_Problem_Objects_brandes_griesse<RealT> >(parlist, comm);
  HDSA::Ptr<HDSA::Weight_Matrices<RealT> > weight_matrices_factory = HDSA::makePtr<Weight_Matrices_brandes_griesse<RealT> >(parlist, parlist_sensitivity);
  HDSA::Ptr<HDSA::Parameter_Sampler<RealT> > sampler = HDSA::makePtr<Parameter_Sampler_brandes_griesse<RealT> >();
  HDSA::Sample_Local_Sensitivities<RealT>(comm,parlist_sensitivity,OP_Objects_Factory,weight_matrices_factory,sampler);  

  return 0;
}
