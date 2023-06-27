#include "Teuchos_GlobalMPISession.hpp"
#include <fstream>

#include "../../src/source_file.hpp"
#include "MD_Interface_Poisson_1D_Test.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {

  HDSA::nullstream bhs;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  HDSA::Ptr<const HDSA::Comm<int> > comm = HDSA::makePtr<HDSA::Comm<int> >();

  HDSA::Ptr<HDSA::Model_Discrepancy_Interface<RealT> > md_interface = HDSA::makePtr<Model_Discrepancy_Interface_Poisson_1D_Test<RealT> >();
  HDSA::Ptr<HDSA::Model_Discrepancy_Update<RealT> > md_update = HDSA::makePtr<HDSA::Model_Discrepancy_Update<RealT> >(md_interface);
  
  RealT alpha = 1.e-3;
  md_update->Compute_Posterior_Data(alpha);
  HDSA::Ptr<HDSA::Vector<RealT> > z_update = md_update->Posterior_Update_Mean();

  const Model_Discrepancy_Interface_Poisson_1D_Test<RealT> md_interface_example = dynamic_cast<const Model_Discrepancy_Interface_Poisson_1D_Test<RealT>&>(*md_interface);
  HDSA::Ptr<HDSA::Vector<RealT> > matlab_z_update  = md_interface_example.Load_Matlab_z_Update();
  matlab_z_update->axpy(-1.0,*z_update);
  RealT error = matlab_z_update->norm();
  std::cout << "Difference between computed solution and stored solution from Matlab = " << error << std::endl;
  // We shouldn't expect this difference to be identically zero, but it should be small, for instance, O(10^-5)

  return 0;
}
