#define Is_hessVec_z_z_Implemented false
#define Is_hessVec_z_theta_Implemented false
#define Is_hessVec_theta_z_Implemented false
#define Is_Misfit_hessVec_z_z_Implemented false
#define Is_Regularization_hessVec_z_z_Implemented false

#include "../../src/source_file_matlab_adapter.hpp"
#include "Parameter_Sampler_ode_control.hpp"
#include "Weight_Matrices_ode_control.hpp"
#include "Opt_Problem_Objects_ode_control.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {
        
  HDSA::nullstream bhs;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  HDSA::Ptr<const HDSA::Comm<int> > comm = HDSA::makePtr<HDSA::Comm<int> >();
  /*** Read in XML input ***/
  std::string filenameSensitivity = "Sensitivity_input.xml";
  HDSA::Ptr<HDSA::ParameterList> parlist_sensitivity = HDSA::makePtr<HDSA::ParameterList>();
  HDSA::updateParametersFromXmlFile( filenameSensitivity, *parlist_sensitivity );
  
  long unsigned int z_dim = 10;
  long unsigned int theta_dim = 8;

  HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects_Factory = HDSA::makePtr<Opt_Problem_Objects_ode_control<RealT> >(z_dim);
  HDSA::Ptr<HDSA::Parameter_Sampler<RealT> > sampler = HDSA::makePtr<Parameter_Sampler_ode_control<RealT> >(theta_dim);
  HDSA::Ptr<HDSA::Weight_Matrices<RealT> > weight_matrices = HDSA::makePtr<Weight_Matrices_ode_control<RealT> >(parlist_sensitivity);
  HDSA::Sample_Local_Sensitivities<RealT>(comm,parlist_sensitivity,OP_Objects_Factory,weight_matrices,sampler); 
}
