#include "hdsa_test_problem.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_StatusTest.hpp"
#include "ROL_CompositeStep.hpp"
#include "ROL_ConstraintStatusTest.hpp"
#include "ROL_ParameterList.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iostream>
#include <fstream>
#include <math.h>
#include "/gpfs/joshart/Trilinos/hdsalib/refactor/src/source_file.hpp"
#include "HDSA_ROL_MD_Interface_hdsa_test_problem.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);
  
  std::string filename = "input.xml";
  auto parlist = ROL::getParametersFromXmlFile( filename );
  
  int m = 51;
  Constraint_HdsaTestProb<RealT> con(m);
  Objective_HdsaTestProb<RealT> obj(m);
  ROL::Ptr<ROL::Objective_SimOpt<RealT> > pobj  = ROL::makePtrFromRef(obj);
  ROL::Ptr<ROL::Constraint_SimOpt<RealT> > pcon = ROL::makePtrFromRef(con);

  ROL::Ptr<std::vector<RealT> > z_ptr    = ROL::makePtr<std::vector<RealT> >(m, 0.0);
  ROL::StdVector<RealT> z(z_ptr);
  ROL::Ptr<ROL::Vector<RealT> > zp  = ROL::makePtrFromRef(z);
  ROL::Ptr<std::vector<RealT> > u_ptr  = ROL::makePtr<std::vector<RealT> >(m, 0.0);
  ROL::StdVector<RealT> u(u_ptr);
  ROL::Ptr<ROL::Vector<RealT> > up  = ROL::makePtrFromRef(u);

  HDSA::Ptr<HDSA::Model_Discrepancy_Interface<RealT> > md_interface = HDSA::makePtr<ROL_Model_Discrepancy_Interface_hdsa_test_problem<RealT> >(pobj, pcon, up, zp, m);
  HDSA::Ptr<HDSA::Model_Discrepancy_Update<RealT> > md_update = HDSA::makePtr<HDSA::Model_Discrepancy_Update<RealT> >(md_interface);
  
  RealT alpha = 1.e-3;
  md_update->Compute_Posterior_Data(alpha);
  HDSA::Ptr<HDSA::Vector<RealT> > z_update = md_update->Posterior_Update_Mean();

  HDSA::ROL_Vector<RealT>& z_update_rol = dynamic_cast<HDSA::ROL_Vector<RealT>&>(*z_update);
  ROL::Ptr<std::vector<RealT> > z_update_std = dynamic_cast<ROL::StdVector<RealT>&>(*z_update_rol.rol_vec).getVector();

  std::string name = "z_update.txt";
  std::ofstream fout;
  fout.open(name);
  for(int k = 0; k < z_update->dimension(); k++)
    {
      fout << std::setprecision(16) << (*z_update_std)[k] << std::endl;
    }
  fout.close();

  return 0;
}
