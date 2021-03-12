#include "Teuchos_GlobalMPISession.hpp"

#include "../../src/source_file.hpp"

#include "full_objective_CDR_model_error_test.hpp"
#include "reduced_objective_CDR_model_error_test.hpp"
#include "constraint_CDR_model_error_test.hpp"

#include "Weight_Matrices_CDR_model_error_test.hpp"
#include "Opt_Problem_Objects_CDR_model_error_test.hpp"
#include "Parameter_Sampler_CDR_model_error_test.hpp"
#include "Model_Error_Objects_CDR_model_error_test.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {

  HDSA::nullstream bhs;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  HDSA::Ptr<const HDSA::Comm<int> > comm = HDSA::makePtr<HDSA::Comm<int> >();
  /*** Read in XML input ***/
  std::string filenameSensitivity = "Sensitivity_input.xml";
  HDSA::Ptr<HDSA::ParameterList> parlist_sensitivity = HDSA::makePtr<HDSA::ParameterList>();
  HDSA::updateParametersFromXmlFile( filenameSensitivity, *parlist_sensitivity );

  HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects_Factory = HDSA::makePtr<Opt_Problem_Objects_CDR_model_error_test<RealT> >();
  HDSA::Ptr<HDSA::Weight_Matrices<RealT> > weight_matrices_factory = HDSA::makePtr<Weight_Matrices_CDR_model_error_test<RealT> >(parlist_sensitivity);
  HDSA::Ptr<HDSA::Parameter_Sampler<RealT> > sampler = HDSA::makePtr<Parameter_Sampler_CDR_model_error_test<RealT> >();

  HDSA::Ptr<HDSA::Model_Error_Objects<RealT> > model_error_obj = HDSA::makePtr<Model_Error_Objects_CDR_model_error_test<RealT> >(parlist_sensitivity,OP_Objects_Factory,weight_matrices_factory);
  HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects_Factory_model_error = HDSA::makePtr<HDSA::Opt_Problem_Objects_Model_Error<RealT> >(model_error_obj);
  HDSA::Ptr<HDSA::Weight_Matrices<RealT> > weight_matrices_factory_model_error = HDSA::makePtr<HDSA::Weight_Matrices_Model_Error<RealT> >(model_error_obj);

  HDSA::Sample_Local_Sensitivities<RealT>(comm,parlist_sensitivity,OP_Objects_Factory_model_error,weight_matrices_factory_model_error,sampler);  

  return 0;
}
