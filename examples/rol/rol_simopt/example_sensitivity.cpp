#include "rol_simopt_test_problem.hpp"
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
#include "../../../src/source_file.hpp"
#include "Elliptic_u_Prior_Interface_rol_simopt_test_problem.hpp"
#include "Elliptic_z_Prior_Interface_rol_simopt_test_problem.hpp"
#include "Data_Interface_rol_simopt_test_problem.hpp"

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
  Constraint_SimOptTestProb<RealT> con(m);
  Objective_SimOptTestProb<RealT> obj(m);
  ROL::Ptr<ROL::Objective_SimOpt<RealT> > pobj  = ROL::makePtrFromRef(obj);
  ROL::Ptr<ROL::Constraint_SimOpt<RealT> > pcon = ROL::makePtrFromRef(con);

  ROL::Ptr<std::vector<RealT> > z_ptr    = ROL::makePtr<std::vector<RealT> >(m, 0.0);
  ROL::StdVector<RealT> z(z_ptr);
  ROL::Ptr<ROL::Vector<RealT> > zp  = ROL::makePtrFromRef(z);
  ROL::Ptr<std::vector<RealT> > u_ptr  = ROL::makePtr<std::vector<RealT> >(m, 0.0);
  ROL::StdVector<RealT> u(u_ptr);
  ROL::Ptr<ROL::Vector<RealT> > up  = ROL::makePtrFromRef(u);

  HDSA::Ptr<HDSA::MD_Opt_Prob_Interface<RealT> > opt_prob_interface = HDSA::makePtr<HDSA::MD_ROL_Opt_Prob_Interface<RealT> >(pobj, pcon, up, zp);

  // Need to check hyper-parameter values
  HDSA::Ptr<HDSA::Random_Number_Generator<RealT> > random_number_generator = HDSA::makePtr<HDSA::Random_Number_Generator<RealT> >();
  HDSA::Ptr<HDSA::MD_Data_Interface<RealT> > data_interface = HDSA::makePtr<Data_Interface_SimOptTestProb<RealT> >(m);
  RealT alpha_u = 1.0/std::pow(2.0,2.0);
  RealT alpha_z = 1.0/std::pow(100.00,2.0);
  HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT> > u_prior_interface = HDSA::makePtr<Elliptic_u_Prior_Interface_SimOptTestProb<RealT> >(alpha_u,m);
  HDSA::Ptr<HDSA::MD_z_Prior_Interface<RealT> > z_prior_interface = HDSA::makePtr<Elliptic_z_Prior_Interface_SimOptTestProb<RealT> >(alpha_z,m,random_number_generator);



  
  // HDSA::Ptr<HDSA::Model_Discrepancy_Update<RealT> > md_update = HDSA::makePtr<HDSA::Model_Discrepancy_Update<RealT> >(md_interface);
  
  // RealT alpha = 1.e-3;
  // md_update->Compute_Posterior_Data(alpha);
  // HDSA::Ptr<HDSA::Vector<RealT> > z_update = md_update->Posterior_Update_Mean();

  // HDSA::ROL_Vector<RealT>& z_update_rol = dynamic_cast<HDSA::ROL_Vector<RealT>&>(*z_update);
  // ROL::Ptr<std::vector<RealT> > z_update_std = dynamic_cast<ROL::StdVector<RealT>&>(*z_update_rol.rol_vec).getVector();

  // std::string name = "z_update.txt";
  // std::ofstream fout;
  // fout.open(name);
  // for(int k = 0; k < z_update->dimension(); k++)
  //   {
  //     fout << std::setprecision(16) << (*z_update_std)[k] << std::endl;
  //   }
  // fout.close();

  return 0;
}
