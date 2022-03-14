#include "Teuchos_GlobalMPISession.hpp"

#include "../../../PDE-OPT/TOOLS/solver.cpp"
#include "../../../PDE-OPT/TOOLS/solver_def.hpp"
#include "../../../PDE-OPT/TOOLS/assembler.cpp"
#include "../../../PDE-OPT/TOOLS/assembler_def.hpp"
#include "../../../PDE-OPT/TOOLS/meshmanager.hpp"
#include "../../../PDE-OPT/TOOLS/pdeconstraint.hpp"
#include "../../../PDE-OPT/TOOLS/pdevector.hpp"
#include "../../../PDE-OPT/TOOLS/pdeobjective.hpp"
#include "pde_stokes.hpp"
#include "obj_stokes.hpp"

#include "../../src/source_file.hpp"
#include "Opt_Problem_Objects_stokes.hpp"
#include "Weight_Matrices_stokes.hpp"
#include "Parameter_Sampler_stokes.hpp"
#include "Model_Error_Objects_stokes.hpp"

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

  HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects_Factory = HDSA::makePtr<Opt_Problem_Objects_stokes<RealT> >(parlist,comm);
  HDSA::Ptr<HDSA::Weight_Matrices<RealT> > weight_matrices_factory = HDSA::makePtr<Weight_Matrices_stokes<RealT> >(parlist_sensitivity);
  HDSA::Ptr<HDSA::Parameter_Sampler<RealT> > sampler = HDSA::makePtr<Parameter_Sampler_stokes<RealT> >();

  HDSA::Ptr<HDSA::Bayes_Model_Error_Objects<RealT> > bayes_model_error_obj 
    = HDSA::makePtr<Bayes_Model_Error_Objects_stokes<RealT> >(parlist_sensitivity,OP_Objects_Factory,weight_matrices_factory);
  HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects_Factory_model_error = HDSA::makePtr<HDSA::Opt_Problem_Objects_Bayes_Model_Error<RealT> >(bayes_model_error_obj);

  HDSA::Sample_Local_Sensitivities<RealT>(comm,parlist_sensitivity,OP_Objects_Factory_model_error,weight_matrices_factory,sampler);  

  return 0;
}
