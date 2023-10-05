#include "rol_simopt_test_problem.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_StatusTest.hpp"
#include "ROL_ConstraintStatusTest.hpp"
#include "ROL_ParameterList.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iostream>
#include <fstream>
#include <math.h>

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

  int errorFlag  = 0;
  
  std::string filename = "input.xml";
  auto parlist = ROL::getParametersFromXmlFile( filename );
  
  int m = 51;
  Constraint_SimOptTestProb<RealT> con(m);
  Objective_SimOptTestProb<RealT> obj(m);
    
  // Initialize iteration vectors.
  ROL::Ptr<std::vector<RealT> > z_ptr    = ROL::makePtr<std::vector<RealT> >(m, 0.0);
  ROL::Ptr<std::vector<RealT> > yz_ptr   = ROL::makePtr<std::vector<RealT> >(m, 0.0);
  ROL::Ptr<std::vector<RealT> > soln_ptr = ROL::makePtr<std::vector<RealT> >(m, 0.0);
  ROL::StdVector<RealT> z(z_ptr);
  ROL::StdVector<RealT> yz(yz_ptr);
  ROL::StdVector<RealT> soln(soln_ptr);
  ROL::Ptr<ROL::Vector<RealT> > zp  = ROL::makePtrFromRef(z);
  ROL::Ptr<ROL::Vector<RealT> > yzp = ROL::makePtrFromRef(yz);

  ROL::Ptr<std::vector<RealT> > u_ptr  = ROL::makePtr<std::vector<RealT> >(m, 0.0);
  ROL::Ptr<std::vector<RealT> > yu_ptr = ROL::makePtr<std::vector<RealT> >(m, 0.0);
  ROL::StdVector<RealT> u(u_ptr);
  ROL::StdVector<RealT> yu(yu_ptr);
  ROL::Ptr<ROL::Vector<RealT> > up  = ROL::makePtrFromRef(u);
  ROL::Ptr<ROL::Vector<RealT> > yup = ROL::makePtrFromRef(yu);

  ROL::Ptr<std::vector<RealT> > jv_ptr  = ROL::makePtr<std::vector<RealT>>(m, 1.0);
  ROL::StdVector<RealT> jv(jv_ptr);
  ROL::Ptr<ROL::Vector<RealT> > jvp = ROL::makePtrFromRef(jv);

  ROL::Vector_SimOpt<RealT> x(up,zp);
  ROL::Vector_SimOpt<RealT> y(yup,yzp);

  // Check derivatives
  x.randomize();
  y.randomize();
  obj.checkGradient(x,x,y,true,*outStream);
  obj.checkHessVec(x,x,y,true,*outStream);

  con.checkApplyJacobian(x,y,jv,true,*outStream);
  con.checkApplyAdjointJacobian(x,yu,jv,x,true,*outStream);
  con.checkApplyAdjointHessian(x,yu,y,x,true,*outStream);
  // Check consistency of Jacobians and adjoint Jacobians.
  con.checkAdjointConsistencyJacobian_1(jv,yu,u,z,true,*outStream);
  con.checkAdjointConsistencyJacobian_2(jv,yz,u,z,true,*outStream);
  // Check consistency of solves.
  con.checkSolve(u,z,jv,true,*outStream);
  con.checkInverseJacobian_1(jv,yu,u,z,true,*outStream);
  con.checkInverseAdjointJacobian_1(yu,jv,u,z,true,*outStream);

  // Initialize reduced objective function.
  ROL::Ptr<std::vector<RealT> > p_ptr  = ROL::makePtr<std::vector<RealT> >(m, 0.0);
  ROL::StdVector<RealT> p(p_ptr);
  ROL::Ptr<ROL::Vector<RealT> > pp  = ROL::makePtrFromRef(p);
  ROL::Ptr<ROL::Objective_SimOpt<RealT> > pobj  = ROL::makePtrFromRef(obj);
  ROL::Ptr<ROL::Constraint_SimOpt<RealT> > pcon = ROL::makePtrFromRef(con);
  ROL::Reduced_Objective_SimOpt<RealT> robj(pobj,pcon,up,zp,pp);
  // Check derivatives.
  *outStream << "Derivatives of reduced objective" << std::endl;
  robj.checkGradient(z,z,yz,true,*outStream);
  robj.checkHessVec(z,z,yz,true,*outStream);
     
  // Trust Region
  ROL::Ptr<ROL::Step<RealT> > step = ROL::makePtr<ROL::TrustRegionStep<RealT> >(*parlist);
  ROL::Ptr<ROL::StatusTest<RealT> > status = ROL::makePtr<ROL::StatusTest<RealT> >(*parlist);
  ROL::Algorithm<RealT> algo_tr(step,status,false);
  std::clock_t timer_tr = std::clock();
  algo_tr.run(z,robj,true,*outStream);
  *outStream << "Trust-Region required " << (std::clock()-timer_tr)/(RealT)CLOCKS_PER_SEC
	     << " seconds.\n";

  ROL::Ptr<const std::vector<RealT> > z_std =
    dynamic_cast<const ROL::StdVector<RealT>&>(z).getVector();
  std::cout << "Optimal solution is:" << std::endl;
  for (int i = 0; i < m; i++)
    {
      std::cout << (*z_std)[i] << std::endl;
    }

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}
