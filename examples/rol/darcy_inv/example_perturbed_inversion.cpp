#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_TpetraMultiVector.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "ROL_LinearCombinationObjective.hpp"

#include "../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/meshmanager.hpp"
#include "../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/pdeconstraint.hpp"
#include "../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/pdeobjective.hpp"
#include "../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/pdevector.hpp"
#include "pde_darcy_flow.hpp"
#include "obj_darcy_flow.hpp"
#include "mesh_darcy_flow.hpp"
#include "elliptic_prior_reg_obj.hpp"
#include "../../../src/source_file.hpp"
#include "../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/solver.cpp"
#include "../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/solver_def.hpp"
#include "../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/assembler.cpp"
#include "../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/assembler_def.hpp" 

typedef double RealT;

template<class RealT>
void Load_Nominal_Solution(HDSA::Ptr<Tpetra::MultiVector<> > & z_ptr) 
{
  int num_coeff_load = z_ptr->getGlobalLength();
  
  // read in data
  std::ifstream in("z_bar.txt");          
  // read the elements in the file into a vector  
  // test file open   
  RealT val = 0.0;
  if (in) 
    {   
      for(int j = 0; j < num_coeff_load; j++)
      {
        in >> val;
        z_ptr->replaceGlobalValue(j,0,val);
      }
    }
  else
    {
      std::cout << "Error loading the data from z_bar.txt" << std::endl;
    }  
}

template<class RealT>
void Load_Perturbed_Parameters(std::vector<RealT> & param) 
{
  int theta_dim = param.size();
  
  // read in data
  std::ifstream in("theta_star.txt");          
  // read the elements in the file into a vector  
  // test file open   
  RealT val = 0.0;
  if (in) 
    {   
      for(int j = 0; j < theta_dim; j++)
      {
        in >> val;
        param[j] = val;
      }
    }
  else
    {
      std::cout << "Error loading the data from theta_star.txt" << std::endl;
    }  
}

int main(int argc, char *argv[]) {

  HDSA::nullstream bhs;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  HDSA::Ptr<const HDSA::Comm<int> > comm = HDSA::makePtr<HDSA::Comm<int> >();
  HDSA::Ptr<std::ostream> outStream;
  int myRank = comm->getRank();
  if(myRank == 0)
    {
	    outStream = HDSA::makePtrFromRef(std::cout);
    }
  else
    {	
	    outStream =  HDSA::makePtrFromRef(bhs);
    }
  std::string filename = "input.xml";
  HDSA::Ptr<HDSA::ParameterList> parlist = HDSA::makePtr<HDSA::ParameterList>();
  HDSA::updateParametersFromXmlFile( filename, *parlist );

  /*************************************************************************/
  /***************** BUILD GOVERNING PDE ***********************************/
  /*************************************************************************/
  /*** Initialize main data structure. ***/
  HDSA::Ptr<MeshManager<RealT> > meshMgr = HDSA::makePtr<MeshManager_darcy_flow<RealT> >(*parlist);
  HDSA::Ptr<PDE_darcy_flow<RealT> > pde = HDSA::makePtr<PDE_darcy_flow<RealT> >(*parlist);
  HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > con = HDSA::makePtr<PDE_Constraint<RealT> >(pde,meshMgr,comm->Get_Teuchos_Communicator(),*parlist,*outStream);
  HDSA::Ptr<PDE_Constraint<RealT> >pdecon = HDSA::dynamicPtrCast<PDE_Constraint<RealT> >(con);
  HDSA::Ptr<Assembler<RealT> > assembler = pdecon->getAssembler();
  
  int theta_modes = parlist->sublist("Problem").get("Uncertain Modes per Dimension", 1);
  std::vector<RealT> param = std::vector<RealT>(std::pow(theta_modes,2),0.0);
  Load_Perturbed_Parameters<RealT>(param);
  pdecon->setParameter(param);
  con->setSolveParameters(*parlist);
  
  /*************************************************************************/
  /***************** BUILD VECTORS *****************************************/
  /*************************************************************************/
  // Create state vectors
  HDSA::Ptr<Tpetra::MultiVector<> > u_ptr, p_ptr, r_ptr, z_ptr;
  HDSA::Ptr<ROL::Vector<RealT> > up, zp, pp, rp;
  u_ptr  = assembler->createStateVector();
  up  = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(u_ptr,pde,assembler,*parlist);
  p_ptr  = assembler->createStateVector();
  pp  = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(p_ptr,pde,assembler,*parlist);
  // Create residual vector
  r_ptr  = assembler->createStateVector();
  rp  = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(r_ptr,pde,assembler,*parlist);
  // Create control vectors
  z_ptr  = assembler->createControlVector();
  zp = HDSA::makePtr<PDE_PrimalOptVector<RealT> >(z_ptr,pde,assembler,*parlist);

   /*************************************************************************/
  /***************** BUILD COST FUNCTIONAL *********************************/
  /*************************************************************************/ 

  int num_obs = parlist->sublist("Problem").get("Number of Observations", 1);
  std::vector<RealT> target_data = std::vector<RealT>(num_obs,0.0);
  std::vector<int> target_data_ids = std::vector<int>(num_obs,0);

  // read in data
  std::ifstream in_data("Data_Generation/obs_data.txt");          
  // read the elements in the file into a vector  
  // test file open   
  if (in_data) 
    {   
      for(int j = 0; j < num_obs; j++)
	    {
	      in_data >> target_data[j];
	    }
    }
  else
    {
      std::cout << "Error loading the data from obs_data.txt" << std::endl;
    }  

  // read in data
  std::ifstream in_loc("Data_Generation/obs_locations.txt");          
  // read the elements in the file into a vector  
  // test file open   
  if (in_loc) 
    {   
      for(int j = 0; j < num_obs; j++)
	    {
	      in_loc >> target_data_ids[j];
        target_data_ids[j] = target_data_ids[j] - 1;
	    }
    }
  else
    {
      std::cout << "Error loading the data from obs_locations.txt" << std::endl;
    }  

  HDSA::Ptr<ROL::Objective_SimOpt<RealT> > obj_misfit = HDSA::makePtr<State_Cost_darcy_flow<RealT> >(target_data, target_data_ids, up);

  HDSA::Ptr<ROL::VectorController<RealT> > stateStore = HDSA::makePtr<ROL::VectorController<RealT> >();
  HDSA::Ptr<ROL::Objective<RealT> > robj_misfit = HDSA::makePtr<ROL::Reduced_Objective_SimOpt<RealT> >(obj_misfit, con, stateStore, up, zp, pp, true, false);
  HDSA::Ptr<ROL::Objective<RealT> > robj_reg = HDSA::makePtr<Elliptic_Prior_Regularization_Objective<RealT> >(comm->Get_Teuchos_Communicator(), parlist, outStream);

  std::vector<RealT> weights = std::vector<RealT>(2);
  weights[0] = parlist->sublist("Problem").get("State Cost",1.0);
  weights[1] = 1.0;
  std::vector<HDSA::Ptr<ROL::Objective<RealT> > > obj_vec;
  obj_vec.resize(2);
  obj_vec[0] = robj_misfit;
  obj_vec[1] = robj_reg;
  HDSA::Ptr<ROL::Objective<RealT> > robj = HDSA::makePtr<ROL::LinearCombinationObjective<RealT> >(weights,obj_vec);
  
  // Run derivative checks
  bool checkDeriv = parlist->sublist("Problem").get("Check Derivatives",false);
  if ( checkDeriv ) { 
    // Create state vector and set to zeroes
    HDSA::Ptr<Tpetra::MultiVector<> > du_ptr = assembler->createStateVector();     du_ptr->randomize();
    HDSA::Ptr<ROL::Vector<RealT> > dup = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(du_ptr,pde,assembler);
      
    // Create control vectors
    HDSA::Ptr<Tpetra::MultiVector<> > dz_ptr  = assembler->createControlVector();
    HDSA::Ptr<ROL::Vector<RealT> > dzp = HDSA::makePtr<PDE_PrimalOptVector<RealT> >(dz_ptr,pde,assembler,*parlist);
      
    // Create ROL SimOpt vectors
    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> d(dup,dzp);
      
    up->randomize();
    zp->randomize();
    pp->randomize();
    dup->randomize();
    dzp->randomize();
    rp->randomize();
    pp->randomize();

    *outStream << "Check Gradient of Full Misfit Objective Function" << std::endl;
    obj_misfit->checkGradient(x,d,true,*outStream);
    *outStream << std::endl << "Check Hessian of Full Misfit Objective Function" << std::endl;
    obj_misfit->checkHessVec(x,d,true,*outStream);
    *outStream << std::endl << "Check Jacobian of Constraint" << std::endl;
    con->checkApplyJacobian(x,d,*up,true,*outStream);
    *outStream << std::endl << "Check Jacobian_1 of Constraint" << std::endl;
    con->checkApplyJacobian_1(*up,*zp,*dup,*rp,true,*outStream);
    *outStream << std::endl << "Check Jacobian_2 of Constraint" << std::endl;
    con->checkApplyJacobian_2(*up,*zp,*dzp,*rp,true,*outStream);
    *outStream << std::endl << "Check Hessian of Constraint" << std::endl;
    con->checkApplyAdjointHessian(x,*dup,d,x,true,*outStream);
    *outStream << std::endl << "Check Hessian_11 of Constraint" << std::endl;
    con->checkApplyAdjointHessian_11(*up,*zp,*pp,*dup,*rp,true,*outStream);
    *outStream << std::endl << "Check Hessian_12 of Constraint" << std::endl;
    con->checkApplyAdjointHessian_12(*up,*zp,*pp,*dup,*dzp,true,*outStream);
    *outStream << std::endl << "Check Hessian_21 of Constraint" << std::endl;
    con->checkApplyAdjointHessian_21(*up,*zp,*pp,*dzp,*rp,true,*outStream);
    *outStream << std::endl << "Check Hessian_22 of Constraint" << std::endl;
    con->checkApplyAdjointHessian_22(*up,*zp,*pp,*dzp,*dzp,true,*outStream);
      
    *outStream << std::endl << "Check Adjoint Jacobian of Constraint" << std::endl;
    con->checkAdjointConsistencyJacobian(*dup,d,x,true,*outStream);
    *outStream << std::endl << "Check Adjoint Jacobian_1 of Constraint" << std::endl;
    con->checkAdjointConsistencyJacobian_1(*pp,*dup,*up,*zp,true,*outStream);
    *outStream << std::endl << "Check Adjoint Jacobian_2 of Constraint" << std::endl;
    con->checkAdjointConsistencyJacobian_2(*pp,*dzp,*up,*zp,true,*outStream);
      
    *outStream << std::endl << "Check Constraint Solve" << std::endl;
    con->checkSolve(*up,*zp,*rp,true,*outStream);
    *outStream << std::endl << "Check Inverse Jacobian_1 of Constraint" << std::endl;
    con->checkInverseJacobian_1(*rp,*dup,*up,*zp,true,*outStream);
    *outStream << std::endl << "Check Inverse Adjoint Jacobian_1 of Constraint" << std::endl;
    con->checkInverseAdjointJacobian_1(*rp,*pp,*up,*zp,true,*outStream);
      
    *outStream << std::endl << "Check Gradient of Reduced Objective Function" << std::endl;
    robj->checkGradient(*zp,*dzp,true,*outStream);
    *outStream << std::endl << "Check Hessian of Reduced Objective Function" << std::endl;
    robj->checkHessVec(*zp,*dzp,true,*outStream); 
  }
	
  // Set initial vector
  Load_Nominal_Solution<RealT>(z_ptr);  

  // Build optimization problem and check derivatives
  ROL::OptimizationProblem<RealT> optProb(robj,zp);
  // Build optimization solver and solve
  ROL::OptimizationSolver<RealT> optSolver(optProb,*parlist);
  std::clock_t timer = std::clock();
  optSolver.solve(*outStream);
  *outStream << "Trust Region Time: "
	      << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC)
	      << " seconds." << std::endl << std::endl;
  
  RealT tol = 1.e-8;
  con->solve(*pp,*up,*zp,tol);
  pdecon->outputTpetraVector(u_ptr,"optimal_u_perturbed.txt");
  pdecon->outputTpetraVector(z_ptr,"optimal_z_perturbed.txt");
  pdecon->outputTpetraData();
  assembler->printMeshData(*outStream);

  return 0;
}
