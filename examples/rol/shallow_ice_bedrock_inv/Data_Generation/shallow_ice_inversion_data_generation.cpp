#include "Teuchos_Comm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include "ROL_Stream.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_Solver.hpp"
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

#include "HDSA_Stream.hpp"
#include "HDSA_Comm.hpp"
#include "HDSA_Ptr.hpp"
#include "HDSA_ParameterList.hpp"

#include "../../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/meshmanager.hpp"
#include "../../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/dynconstraint.hpp"
#include "../../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/pdeconstraint.hpp"
#include "../../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/pdeobjective.hpp"
#include "../../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/ltiobjective.hpp"
#include "../../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/pdevector.hpp"
#include "../../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/solver.cpp"
#include "../../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/solver_def.hpp"
#include "../../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/assembler.cpp"
#include "../../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/assembler_def.hpp"

#include "../mesh_shallow_ice.hpp"
#include "../pde_shallow_ice.hpp"
#include "../dynpde_shallow_ice.hpp"


template<class RealT>
void Set_bed(HDSA::Ptr<Tpetra::MultiVector<> > & z_ptr, Teuchos::RCP<Teuchos::ParameterList> & parlist)
{
  int num_coeff_load = parlist->sublist("Problem").get("Number of Coefficients in Loaded Fields", 10); 

  

  std::vector<RealT> bed_coeff = std::vector<RealT>(num_coeff_load);
  // read in data
  std::ifstream in("Bedrock_Topography.txt");          
  // read the elements in the file into a vector  
  // test file open   
  if (in) 
    {   
      for(int j = 0; j < num_coeff_load; j++)
	{
	  in >> bed_coeff[j];
	}
    }
  else
    {
      std::cout << "Error loading the data from Bedrock_Topography.txt" << std::endl;
    }  

  for(int k = 0; k < num_coeff_load; k++)
    {
      z_ptr->replaceGlobalValue(3*k,0,bed_coeff[k]);
    }
}

template<class RealT>
void Set_Initial_Condition(HDSA::Ptr<Tpetra::MultiVector<> > & u_ptr, const HDSA::Ptr<Teuchos::ParameterList> & parlist)
{
  int num_coeff_load = parlist->sublist("Problem").get("Number of Coefficients in Loaded Fields", 10); 

  std::vector<RealT> initial_iter_coeff = std::vector<RealT>(num_coeff_load);
  // read in data
  std::ifstream in("Surface_Height.txt");          
  // read the elements in the file into a vector  
  // test file open   
  if (in) 
    {   
      for(int j = 0; j < num_coeff_load; j++)
	{
	  in >> initial_iter_coeff[j];
	}
    }
  else
    {
      std::cout << "Error loading the data from Surface_Height.txt" << std::endl;
    }  
  
  int nx = parlist->sublist("Geometry").get("NX",0);
  int ny = parlist->sublist("Geometry").get("NY",0);
  RealT nx_float = static_cast<RealT>(nx);
  RealT ny_float = static_cast<RealT>(ny);
  std::vector<RealT> coords = std::vector<RealT>(2,0.0);

  for(int k = 0; k < num_coeff_load; k++)
    {
      u_ptr->replaceGlobalValue(3*k,0,initial_iter_coeff[k]);
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
  int nt         = parlist->sublist("Time Discretization").get("Number of Time Steps", 100);
  RealT T        = parlist->sublist("Time Discretization").get("End Time",             1.0);
  RealT dt       = T/static_cast<RealT>(nt);
  
  /*************************************************************************/
  /***************** BUILD GOVERNING PDE ***********************************/
  /*************************************************************************/
  /*** Initialize main data structure. ***/
  HDSA::Ptr<MeshManager<RealT> > meshMgr
    = HDSA::makePtr<MeshManager_shallow_ice<RealT> >(*parlist);
  // Initialize PDE
  HDSA::Ptr<DynamicPDE_shallow_ice<RealT> > pde
    = HDSA::makePtr<DynamicPDE_shallow_ice<RealT> >(*parlist);
  
  /*************************************************************************/
  /***************** BUILD CONSTRAINT **************************************/
  /*************************************************************************/
  HDSA::Ptr<DynConstraint<RealT> > dyn_con
    = HDSA::makePtr<DynConstraint<RealT> >(pde,meshMgr,comm,*parlist,*outStream);
  dyn_con->getAssembler()->printMeshData(*outStream);
  
  /*************************************************************************/
  /***************** BUILD STATE VECTORS ***********************************/
  /*************************************************************************/
  HDSA::Ptr<Tpetra::MultiVector<> > u0_ptr, uo_ptr, un_ptr,z_ptr, ck_ptr;
  u0_ptr = dyn_con->getAssembler()->createStateVector();
  uo_ptr = dyn_con->getAssembler()->createStateVector();
  un_ptr = dyn_con->getAssembler()->createStateVector();
  z_ptr = dyn_con->getAssembler()->createControlVector();
  ck_ptr = dyn_con->getAssembler()->createResidualVector();
  HDSA::Ptr<ROL::Vector<RealT> > u0, uo, un, ck, z;
  u0 = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(u0_ptr,pde,*dyn_con->getAssembler());
  uo = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(uo_ptr,pde,*dyn_con->getAssembler());
  un = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(un_ptr,pde,*dyn_con->getAssembler());
  ck = HDSA::makePtr<PDE_DualSimVector<RealT> >(ck_ptr,pde,*dyn_con->getAssembler());
  z = HDSA::makePtr<PDE_PrimalOptVector<RealT> >(z_ptr,pde,*dyn_con->getAssembler());
  u0->zero();
  Set_Initial_Condition<RealT>(u0_ptr,parlist);

  dyn_con->setSolveParameters(*parlist);
  
  /*************************************************************************/
  /***************** BUILD REDUCED COST FUNCTIONAL *************************/
  /*************************************************************************/
  std::vector<ROL::TimeStamp<RealT> > timeStamp(nt);
  for( int k=0; k<nt; ++k ) {
    timeStamp.at(k).t.resize(2);
    timeStamp.at(k).t.at(0) = k*dt;
    timeStamp.at(k).t.at(1) = (k+1)*dt;
  }

  // Run derivative checks                                                                                                                                                                                
  bool checkDeriv = parlist->sublist("Problem").get("Check Derivatives",true);
  if ( checkDeriv )
    {
      uo->randomize(2900.0,3200.0); un->randomize(2900.0,3200.0); z->randomize(-18.0,-15.0);
      ROL::ValidateFunction<RealT> validate(1,13,20,11,true,*outStream);
      ROL::DynamicConstraintCheck<RealT>::check(*dyn_con,validate,*uo,*un,*z);
    }

  /*************************************************************************/
  /***************** GENERATE DATA *****************************************/
  /*************************************************************************/
  
  // 1. Generate data with a FE mesh twice the resolution desired for the inverse problem
  // 2. Run Sparcify_Data.m matlab code to remove unwanted data and write to data.txt
  std::vector<std::vector<RealT> > clean_data;
  clean_data.resize(nt);
  Set_bed<RealT>(z_ptr,parlist);
  int dim = u0->dimension();
  uo->set(*u0); un->zero();
  clean_data[0].resize(dim);
  for(int i = 0; i < dim; i++)
    {
      clean_data[0][i] = uo->dot(*un->basis(i));
    }
  for(int k = 1; k < nt; k++)
    {
      // Advance time stepper
      dyn_con->solve(*ck, *uo, *un, *z, timeStamp[k-1]);
      uo->set(*un);
      clean_data[k].resize(dim);
      for(int i = 0; i < dim; i++)
  	{
	  clean_data[k][i] = uo->dot(*un->basis(i));
  	}
    }

  std::string name;
  std::ofstream fout;
  name="clean_true_state.txt";
  fout.open(name);
  for(int i=0; i < nt; i++)
    {
      for(int j=0; j < dim; j++)
  	{
  	  fout << std::setprecision(16) << clean_data[i][j] << "  " << std::setw(20);
  	}
      fout << " " << std::endl;
    }
  fout.close();

  dyn_con->outputTpetraVector(z_ptr,"true_bed.txt");
  dyn_con->outputTpetraData();
  dyn_con->getAssembler()->printMeshData(*outStream);

  return 0;
}
