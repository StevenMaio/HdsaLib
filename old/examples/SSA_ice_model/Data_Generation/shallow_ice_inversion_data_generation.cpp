#include "Teuchos_Comm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include "ROL_Stream.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_Objective_SimOpt.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "ROL_ReducedDynamicObjective.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_DynamicConstraintCheck.hpp"
#include "ROL_DynamicObjectiveCheck.hpp"

#include <iostream>
#include <random>
#include <unordered_set>
#include <set>

#include "../../../src/source_file.hpp"

#include "../../../../PDE-OPT/TOOLS/solver.cpp"
#include "../../../../PDE-OPT/TOOLS/solver_def.hpp"
#include "../../../../PDE-OPT/TOOLS/assembler.cpp"
#include "../../../../PDE-OPT/TOOLS/assembler_def.hpp"
#include "../../../../PDE-OPT/TOOLS/meshmanager.hpp"
#include "../../../../PDE-OPT/TOOLS/pdeconstraint.hpp"
#include "../../../../PDE-OPT/TOOLS/pdevector.hpp"
#include "../../../../PDE-OPT/TOOLS/pdeobjective.hpp"
#include "../mesh_shallow_ice.hpp"
#include "../pde_shallow_ice.hpp"


template<class RealT>
void Set_beta(HDSA::Ptr<Tpetra::MultiVector<> > & z_ptr, Teuchos::RCP<Teuchos::ParameterList> & parlist)
{
  int num_coeff_load = parlist->sublist("Problem").get("Number of Coefficients in Loaded Fields", 10); 

  std::vector<RealT> beta_coeff = std::vector<RealT>(num_coeff_load);
  // read in data
  std::ifstream in("Log_Basal_Sliding.txt");          
  // read the elements in the file into a vector  
  // test file open   
  if (in) 
    {   
      for(int j = 0; j < num_coeff_load; j++)
	{
	  in >> beta_coeff[j];
	}
    }
  else
    {
      std::cout << "Error loading the data from Log_Basal_Sliding.txt" << std::endl;
    }  

  for(int k = 0; k < num_coeff_load; k++)
    {
      z_ptr->replaceGlobalValue(2*k,0,beta_coeff[k]);
    }
}

int main(int argc, char *argv[]) {
  using RealT = double;

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  HDSA::Ptr<const Teuchos::Comm<int> > comm
    = Tpetra::getDefaultComm();

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  const int myRank = comm->getRank();
  HDSA::nullstream bhs;
  HDSA::Ptr<std::ostream> outStream;
  if(myRank == 0)
    {
      outStream = HDSA::makePtrFromRef(std::cout);
    }
  else
    {	
      outStream =  HDSA::makePtrFromRef(bhs);
    }
  
  /*** Read in XML input ***/
  HDSA::Ptr<Teuchos::ParameterList> parlist = Teuchos::getParametersFromXmlFile("input_data_generation.xml");
  
  /*************************************************************************/
  /***************** BUILD GOVERNING PDE ***********************************/
  /*************************************************************************/
  /*** Initialize main data structure. ***/
  HDSA::Ptr<MeshManager<RealT> > meshMgr
    = HDSA::makePtr<MeshManager_shallow_ice<RealT> >(*parlist);
  // Initialize PDE
  HDSA::Ptr<PDE_shallow_ice<RealT> > pde
    = HDSA::makePtr<PDE_shallow_ice<RealT> >(*parlist);
  
  /*************************************************************************/
  /***************** BUILD CONSTRAINT **************************************/
  /*************************************************************************/
  HDSA::Ptr<PDE_Constraint<RealT> > con
    = HDSA::makePtr<PDE_Constraint<RealT> >(pde,meshMgr,comm,*parlist,*outStream);
  
  /*************************************************************************/
  /***************** BUILD STATE VECTORS ***********************************/
  /*************************************************************************/
  HDSA::Ptr<Tpetra::MultiVector<> > u_ptr, z_ptr, c_ptr;
  u_ptr = con->getAssembler()->createStateVector();
  z_ptr = con->getAssembler()->createControlVector();
  c_ptr = con->getAssembler()->createResidualVector();
  HDSA::Ptr<ROL::Vector<RealT> > u, c, z;
  u = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(u_ptr,pde,con->getAssembler());
  c = HDSA::makePtr<PDE_DualSimVector<RealT> >(c_ptr,pde,con->getAssembler());
  z = HDSA::makePtr<PDE_PrimalOptVector<RealT> >(z_ptr,pde,con->getAssembler());
 
  con->setSolveParameters(*parlist);
  
  // Run derivative checks
  bool checkDeriv = parlist->sublist("Problem").get("Check Derivatives",false);
  if ( checkDeriv ) {
      
    // Create state vector and set to zeroes
    HDSA::Ptr<Tpetra::MultiVector<> > du_ptr = con->getAssembler()->createStateVector();     du_ptr->randomize();
    HDSA::Ptr<ROL::Vector<RealT> > dup = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(du_ptr,pde,con->getAssembler());
    HDSA::Ptr<Tpetra::MultiVector<> > p_ptr = con->getAssembler()->createStateVector();     p_ptr->randomize();
    HDSA::Ptr<ROL::Vector<RealT> > pp = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(p_ptr,pde,con->getAssembler());
      
    // Create control vectors
    HDSA::Ptr<Tpetra::MultiVector<> > dz_ptr = con->getAssembler()->createControlVector();     du_ptr->randomize();
    HDSA::Ptr<ROL::Vector<RealT> > dzp = HDSA::makePtr<PDE_PrimalOptVector<RealT> >(du_ptr,pde,con->getAssembler());
      
    // Create ROL SimOpt vectors
    ROL::Vector_SimOpt<RealT> x(u,z);
    ROL::Vector_SimOpt<RealT> d(dup,dzp);
      
    u->randomize();
    z->randomize();
    c->randomize();
    dup->randomize();
    dzp->randomize();
    pp->randomize();

    *outStream << std::endl << "Check Jacobian of Constraint" << std::endl;
    con->checkApplyJacobian(x,d,*u,true,*outStream);
    *outStream << std::endl << "Check Jacobian_1 of Constraint" << std::endl;
    con->checkApplyJacobian_1(*u,*z,*dup,*c,true,*outStream);
    *outStream << std::endl << "Check Jacobian_2 of Constraint" << std::endl;
    con->checkApplyJacobian_2(*u,*z,*dzp,*c,true,*outStream);
      
    con->outputTpetraData();
    con->getAssembler()->printMeshData(*outStream);

    *outStream << std::endl << "Check Adjoint Jacobian of Constraint" << std::endl;
    con->checkAdjointConsistencyJacobian(*dup,d,x,true,*outStream);
    *outStream << std::endl << "Check Adjoint Jacobian_1 of Constraint" << std::endl;
    con->checkAdjointConsistencyJacobian_1(*c,*dup,*u,*z,true,*outStream);
    *outStream << std::endl << "Check Adjoint Jacobian_2 of Constraint" << std::endl;
    con->checkAdjointConsistencyJacobian_2(*c,*dzp,*u,*z,true,*outStream);
      
    *outStream << std::endl << "Check Constraint Solve" << std::endl;
    con->checkSolve(*u,*z,*c,true,*outStream);
    *outStream << std::endl << "Check Inverse Jacobian_1 of Constraint" << std::endl;
    con->checkInverseJacobian_1(*c,*dup,*u,*z,true,*outStream);
    *outStream << std::endl << "Check Inverse Adjoint Jacobian_1 of Constraint" << std::endl;
    con->checkInverseAdjointJacobian_1(*c,*pp,*u,*z,true,*outStream);
    }

  /*************************************************************************/
  /***************** GENERATE DATA *****************************************/
  /*************************************************************************/
  
  // 1. Generate data with a FE mesh twice the resolution desired for the inverse problem
  // 2. Run Sparcify_Data.m matlab code to remove unwanted data and write to data.txt
  std::vector<RealT> clean_data = std::vector<RealT>(u->dimension(),0.0);
  Set_beta<RealT>(z_ptr,parlist);
  RealT tol = 1.e-8;
  con->solve(*c,*u,*z,tol);
  int dim = u->dimension();
  for(int i = 0; i < dim; i++)
    {
      clean_data[i] = u->dot(*c->basis(i));
    }

  std::string name;
  std::ofstream fout;
  name="clean_true_state.txt";
  fout.open(name);
  for(int i = 0; i < dim; i++)
    {
      fout << std::setprecision(16) << clean_data[i] << std::endl;
    }
  fout.close();

  con->outputTpetraVector(z_ptr,"true_beta.txt");
  con->outputTpetraVector(c_ptr,"residual.txt");

  return 0;
}
