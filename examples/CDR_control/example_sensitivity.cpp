#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_OptimizationSolver.hpp"

#include "../../../PDE-OPT/TOOLS/solver.cpp"
#include "../../../PDE-OPT/TOOLS/solver_def.hpp"
#include "../../../PDE-OPT/TOOLS/assembler.cpp"
#include "../../../PDE-OPT/TOOLS/assembler_def.hpp"
#include "../../../PDE-OPT/TOOLS/meshmanager.hpp"
#include "../../../PDE-OPT/TOOLS/linearpdeconstraint.hpp"
#include "../../../PDE-OPT/TOOLS/pdeconstraint.hpp"
#include "../../../PDE-OPT/TOOLS/pdevector.hpp"
#include "../../../PDE-OPT/TOOLS/pdeobjective.hpp"
#include "../../src/source_file.hpp"
#include "mesh_cdr.hpp"
#include "pde_cdr.hpp"
#include "obj_cdr.hpp"
#include "Opt_Problem_Objects_cdr.hpp"

#include "Weight_Matrices_cdr.hpp"
#include "Parameter_Sampler_cdr.hpp"
#include "Model_Error_Objects_cdr.hpp"

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
  
  HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects_Factory = HDSA::makePtr<Opt_Problem_Objects_CDR<RealT> >(parlist, comm);
  HDSA::Ptr<HDSA::Weight_Matrices<RealT> > weight_matrices_factory = HDSA::makePtr<Weight_Matrices_CDR<RealT> >(parlist, parlist_sensitivity);

  HDSA::Ptr<HDSA::Model_Error_Objects<RealT> > model_error_obj = HDSA::makePtr<Model_Error_Objects_CDR<RealT> >(parlist_sensitivity,OP_Objects_Factory,weight_matrices_factory,parlist);
  HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects_Factory_model_error = HDSA::makePtr<HDSA::Opt_Problem_Objects_Model_Error<RealT> >(model_error_obj);

  HDSA::Ptr<HDSA::Parameter_Sampler<RealT> > sampler = HDSA::makePtr<Parameter_Sampler_CDR<RealT> >();
  HDSA::Sample_Local_Sensitivities<RealT>(comm,parlist_sensitivity,OP_Objects_Factory_model_error,weight_matrices_factory,sampler);  

  return 0;
}
