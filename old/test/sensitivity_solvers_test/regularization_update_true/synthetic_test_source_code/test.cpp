#include "Teuchos_GlobalMPISession.hpp"

#include "../../../../src/source_file.hpp"
#include "reduced_objective_synthetic_test.hpp"
#include "Weight_Matrices_synthetic_test.hpp"
#include "Opt_Problem_Objects_synthetic_test.hpp"
#include "Parameter_Sampler_synthetic_test.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {

  HDSA::nullstream bhs;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  HDSA::Ptr<const HDSA::Comm<int> > comm = HDSA::makePtr<HDSA::Comm<int> >();
  /*** Read in XML input ***/
  std::string filenameSensitivity = "Sensitivity_input.xml";
  HDSA::Ptr<HDSA::ParameterList> parlist_sensitivity = HDSA::makePtr<HDSA::ParameterList>();
  HDSA::updateParametersFromXmlFile( filenameSensitivity, *parlist_sensitivity );

  int a_dim = 10;
  std::vector<RealT> a = std::vector<RealT>(a_dim,0.0);
  for(int k = 0; k < a_dim; k++)
    {
      a[k] = static_cast<RealT>(k+1);
    }

  int theta_dim = 20;  
  int z_dim = 25;

  HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects_Factory = HDSA::makePtr<Opt_Problem_Objects_synthetic_test<RealT> >(a,z_dim);
  HDSA::Ptr<HDSA::Weight_Matrices<RealT> > weight_matrices_factory = HDSA::makePtr<Weight_Matrices_synthetic_test<RealT> >(parlist_sensitivity);
  HDSA::Ptr<HDSA::Parameter_Sampler<RealT> > sampler = HDSA::makePtr<Parameter_Sampler_synthetic_test<RealT> >(theta_dim);
  HDSA::Sample_Local_Sensitivities<RealT>(comm,parlist_sensitivity,OP_Objects_Factory,weight_matrices_factory,sampler);  

  return 0;
}
