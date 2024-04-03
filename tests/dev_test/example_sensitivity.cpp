#include "Teuchos_GlobalMPISession.hpp"
#include <fstream>

#include "../../src/source_file.hpp"
#include "MD_Data_Interface_model_discrepancy_synthetic_test.hpp"
#include "MD_Opt_Prob_Interface_model_discrepancy_synthetic_test.hpp"
#include "MD_u_Prior_Interface_model_discrepancy_synthetic_test.hpp"
#include "MD_z_Prior_Interface_model_discrepancy_synthetic_test.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {

  HDSA::nullstream bhs;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  HDSA::Ptr<const HDSA::Comm<int> > comm = HDSA::makePtr<HDSA::Comm<int> >();
 
  HDSA::Ptr<HDSA::MD_Data_Interface<RealT> > md_data_interface = HDSA::makePtr<MD_Data_Interface_model_discrepancy_synthetic_test<RealT> >();
  md_data_interface->Load_Data();

  // There is a segfault happening on the line below
  HDSA::Ptr<HDSA::MD_Opt_Prob_Interface<RealT> > md_opt_prob_interface = HDSA::makePtr<MD_Opt_Prob_Interface_model_discrepancy_synthetic_test<RealT> >();

  HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT> > md_u_prior_interface = HDSA::makePtr<MD_u_Prior_Interface_model_discrepancy_synthetic_test<RealT> >();
  HDSA::Ptr<HDSA::MD_z_Prior_Interface<RealT> > md_z_prior_interface = HDSA::makePtr<MD_z_Prior_Interface_model_discrepancy_synthetic_test<RealT> >();

  return 0;
}
