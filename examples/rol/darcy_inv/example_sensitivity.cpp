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
#include "pde_darcy_flow_aux_param.hpp"
#include "obj_darcy_flow.hpp"
#include "mesh_darcy_flow.hpp"
#include "elliptic_prior_reg_obj.hpp"
#include "HDSA_Stream.hpp"
#include "HDSA_Comm.hpp"
#include "HDSA_Ptr.hpp"
#include "HDSA_ParameterList.hpp"
#include "HDSA_PC_Quasi_Newton_Preconditioner_LIS.hpp"
#include "HDSA_PC_Pseudo_Time_Continuation.hpp"
#include "../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/solver.cpp"
#include "../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/solver_def.hpp"
#include "../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/assembler.cpp"
#include "../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/assembler_def.hpp" 

#include "PC_Sensitivity_Operator_Interface_darcy_flow.hpp"
#include "PC_LIS_Interface_darcy_flow.hpp"

typedef double RealT;

template<class RealT>
void Load_Nominal_Solution(HDSA::Ptr<Tpetra::MultiVector<> > & z_ptr) 
{
  int num_coeff_load = z_ptr->getGlobalLength();
  
  // read in data
  std::ifstream in("optimal_z.txt");         
  // read the elements in the file into a vector  
  // test file open   
  RealT val = 0.0;
  if (in) 
    {   
      // Skip the first two lines
      std::string line;
      std::getline(in, line); // Skip first line
      std::getline(in, line); // Skip second line
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
void Set_Parameters(HDSA::Ptr<HDSA::Vector<RealT> > & theta, std::string & filename) 
{
  int theta_dim = theta->dimension();
  Std_Vector<RealT>& theta_std = dynamic_cast<Std_Vector<RealT>&>(*theta); 
  
  // read in data
  std::ifstream in("Data_Generation/"+filename);          
  // read the elements in the file into a vector  
  // test file open   
  RealT val = 0.0;
  if (in) 
    {   
      for(int j = 0; j < theta_dim; j++)
      {
        in >> val;
        theta_std.Replace_Element(j,val);
      }
    }
  else
    {
      std::cout << "Error loading the data from " << filename << std::endl;
    }  
}

template<class RealT>
void Set_Parameters(std::vector<RealT> & param, std::string & filename) 
{
  int theta_dim = param.size();
  
  // read in data
  std::ifstream in("Data_Generation/"+filename);          
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
      std::cout << "Error loading the data from " << filename << std::endl;
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
  
  HDSA::Ptr<PDE_darcy_flow_aux_param<RealT> > pde_aux_param = HDSA::makePtr<PDE_darcy_flow_aux_param<RealT> >(*parlist);
  HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > con_aux_param = HDSA::makePtr<PDE_Constraint<RealT> >(pde_aux_param,meshMgr,comm->Get_Teuchos_Communicator(),*parlist,*outStream);
  int theta_modes = parlist->sublist("Problem").get("Uncertain Modes per Dimension", 1);
  std::vector<RealT> param = std::vector<RealT>(std::pow(theta_modes,2),0.0);
  std::string filename_theta_bar = "theta_bar.txt";
  Set_Parameters<RealT>(param,filename_theta_bar);
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

  int theta_dim = param.size();
  HDSA::Ptr<HDSA::Vector<RealT> > theta_bar = HDSA::makePtr<Std_Vector<RealT> >(theta_dim, comm);
  Set_Parameters<RealT>(theta_bar,filename_theta_bar);
  
  Load_Nominal_Solution<RealT>(z_ptr);  
  HDSA::Ptr<HDSA::Vector<RealT> > z_bar = HDSA::makePtr<HDSA::ROL_Vector<RealT> >(zp);

  std::string filename_sensitivity = "Sensitivity_input.xml";
  HDSA::Ptr<HDSA::ParameterList> parlist_sensitivity = HDSA::makePtr<HDSA::ParameterList>();
  HDSA::updateParametersFromXmlFile( filename_sensitivity, *parlist_sensitivity );

  RealT grad_tol = parlist_sensitivity->sublist("Problem").get("Gradient Tolerance", 1.e-7);
  bool use_qn_prec = parlist_sensitivity->sublist("Problem").get("use_qn_prec", true);
  bool print_cg_output = parlist_sensitivity->sublist("Problem").get("print_cg_output", true);
  bool print_cg_iter = parlist_sensitivity->sublist("Problem").get("print_cg_iter", true);
  int rank = parlist_sensitivity->sublist("Problem").get("rank", 0);
  int oversampling = parlist_sensitivity->sublist("Problem").get("oversampling", 0);
  int N_fe = parlist_sensitivity->sublist("Problem").get("N_fe", 1);
  int N_me = parlist_sensitivity->sublist("Problem").get("N_me", 1);
  RealT cg_tol = parlist_sensitivity->sublist("Problem").get("CG Tolerance", 1.e-5);
  int max_cg_iter = parlist_sensitivity->sublist("Problem").get("Maximum CG Iterations", 100);
  bool fd_check = parlist_sensitivity->sublist("Problem").get("Finite Difference Check", false);
  bool use_block_update = parlist_sensitivity->sublist("Problem").get("Use Block Update", true);
  RealT tau = parlist_sensitivity->sublist("Problem").get("Block Update tau", 1.e-5);
  int max_storage = parlist_sensitivity->sublist("Problem").get("Maximum Block Update Storage", 10);

  HDSA::Ptr<HDSA::PC_Sensitivity_Operator_Interface<RealT> > sen_op_interface = HDSA::makePtr<PC_Sensitivity_Operator_Interface_darcy_flow<RealT> >(obj_misfit,robj_reg,con,con_aux_param,up,z_bar,theta_bar,fd_check);
  HDSA::Ptr<HDSA::PC_LIS_Interface<RealT> > lis_interface = HDSA::makePtr<PC_LIS_Interface_darcy_flow<RealT> >(robj_misfit,robj_reg,z_bar,theta_bar);    
  HDSA::Ptr<HDSA::PC_Quasi_Newton_Preconditioner_LIS<RealT> > qn_prec = HDSA::makePtr<HDSA::PC_Quasi_Newton_Preconditioner_LIS<RealT> >(z_bar,theta_bar,lis_interface,use_block_update,tau,max_storage);

  HDSA::Ptr<HDSA::PC_Pseudo_Time_Continuation<RealT> > sen =
    HDSA::makePtr<HDSA::PC_Pseudo_Time_Continuation<RealT> >(z_bar,theta_bar,sen_op_interface,qn_prec, grad_tol,use_qn_prec,print_cg_output,print_cg_iter,cg_tol,max_cg_iter);

  if(rank > 0)
  {
    qn_prec->Compute_Hessian_GEVP(*z_bar, *theta_bar, rank, oversampling);
  }

  HDSA::Ptr<Tpetra::MultiVector<> > z_star_ptr  = assembler->createControlVector();
  HDSA::Ptr<ROL::Vector<RealT> > z_star_rol = HDSA::makePtr<PDE_PrimalOptVector<RealT> >(z_star_ptr,pde,assembler,*parlist);
  HDSA::Ptr<HDSA::Vector<RealT> > z_star = HDSA::makePtr<HDSA::ROL_Vector<RealT> >(z_star_rol);
  
  HDSA::Ptr<HDSA::Vector<RealT> > grad_star = z_bar->clone();
  HDSA::Ptr<HDSA::Vector<RealT> > theta_star = theta_bar->clone();
  std::string filename_theta_star = "theta_star.txt";
  Set_Parameters<RealT>(theta_star,filename_theta_star);
  
  if(N_fe>0)
  {
    sen->Pseudo_Time_Continuation_Forward_Euler(*z_star,*grad_star,*theta_star,N_fe); 
    pdecon->outputTpetraVector(z_star_ptr,"z_star_fe.txt");
  }
  
  if(N_me>0)
  {
    sen->Pseudo_Time_Continuation_Modified_Euler(*z_star,*grad_star,*theta_star,N_me); 
    pdecon->outputTpetraVector(z_star_ptr,"z_star_me.txt");
  }

  return 0;
}
