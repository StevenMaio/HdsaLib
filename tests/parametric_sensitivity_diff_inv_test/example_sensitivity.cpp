#include "Teuchos_GlobalMPISession.hpp"
#include <fstream>

#include "../../src/source_file.hpp"
#include "Adv_Diff_Constraint.hpp"


typedef double RealT;

int main(int argc, char *argv[]) {

  HDSA::nullstream bhs;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  HDSA::Ptr<const HDSA::Comm<int> > comm = HDSA::makePtr<HDSA::Comm<int> >();

  HDSA::Ptr<Adv_Diff_Constraint<RealT> > con = HDSA::makePtr<Adv_Diff_Constraint<RealT> >();
  int m = 100;
  HDSA::Ptr<HDSA::Vector<RealT> > u = HDSA::makePtr<Std_Vector<RealT> >(m);
  HDSA::Ptr<HDSA::Vector<RealT> > z = HDSA::makePtr<Std_Vector<RealT> >(m);
  HDSA::Ptr<HDSA::Vector<RealT> > theta = HDSA::makePtr<Std_Vector<RealT> >(m);
  z->setScalar(1.0);
  theta->setScalar(5.0);
  con->State_Solve(*u,*z,*theta);

  std::string name = "u.txt";
  u->Write_to_File(name);
  
  return 0;
}
