#include "Teuchos_GlobalMPISession.hpp"

#include "./modified_rol_source_code/ROL_ReducedDynamicObjective_Stationary_Control.hpp"

#include "ROL_Stream.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_Objective_SimOpt.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "ROL_ReducedDynamicObjective.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_DynamicConstraintCheck.hpp"
#include "ROL_DynamicObjectiveCheck.hpp"

#include "../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/dynconstraint.hpp"
#include "../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/ltiobjective.hpp"
#include "../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/pdeobjective.hpp"
#include "../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/pdevector.hpp"

#include "./modified_rol_source_code/Objective_SimOpt_TS.hpp"
#include "./modified_rol_source_code/Misfit_Regularization_TS_Objective_SimOpt.hpp"
#include "./modified_rol_source_code/ltiobjective_TS.hpp"
#include "./modified_rol_source_code/Reduced_Objective_Regularization.hpp"

#include "../../../src/source_file.hpp"

#include "mesh_shallow_ice.hpp"
#include "elliptic_prior_reg_obj.hpp"
#include "pde_shallow_ice.hpp"
#include "dynpde_shallow_ice.hpp"
#include "obj_shallow_ice.hpp"

#include "../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/solver.cpp"
#include "../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/solver_def.hpp"
#include "../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/assembler.cpp"
#include "../../../../Trilinos/packages/rol/example/PDE-OPT/TOOLS/assembler_def.hpp"


typedef double RealT;

template<class RealT>
void Set_Initial_Iterate(HDSA::Ptr<Tpetra::MultiVector<> > & z_ptr, const HDSA::Ptr<HDSA::ParameterList> & parlist) 
{
  int num_coeff_load = 3*parlist->sublist("Problem").get("Number of Coefficients in Loaded Fields", 10); 
  std::vector<RealT> initial_iter_coeff = std::vector<RealT>(num_coeff_load);
  // read in data
  std::ifstream in("z_bar.txt");          
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
      std::cout << "Error loading the data from z_bar.txt" << std::endl;
    }  
  
  for(int k = 0; k < num_coeff_load; k++)
    {
      z_ptr->replaceGlobalValue(k,0,initial_iter_coeff[k]);
    }

}

template<class RealT>
void Set_Initial_Condition(HDSA::Ptr<Tpetra::MultiVector<> > & u_ptr, const HDSA::Ptr<HDSA::ParameterList> & parlist) 
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
  
  for(int k = 0; k < num_coeff_load; k++)
    {
      u_ptr->replaceGlobalValue(3*k,0,initial_iter_coeff[k]);
    }
}

int main(int argc, char *argv[]) {

  HDSA::nullstream bhs;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  HDSA::Ptr<const HDSA::Comm<int> > comm = HDSA::makePtr<HDSA::Comm<int> >();
  std::string filename = "input.xml";
  HDSA::Ptr<HDSA::ParameterList> parlist = HDSA::makePtr<HDSA::ParameterList>();
  HDSA::updateParametersFromXmlFile( filename, *parlist );

  int L = parlist->sublist("Problem").get("Number of Uncertain Basis Functions", 10);
  int theta_dim = 2*(L+1)*(L+1);
  HDSA::Ptr<HDSA::Vector<RealT> > theta = HDSA::makePtr<Std_Vector<RealT> >(theta_dim);
  Std_Vector<RealT> &theta_std = dynamic_cast<Std_Vector<RealT> &>(*theta);
  RealT x = 0.0;
  RealT y = 0.0;
  RealT val = 0.0;
  for(int k = 0; k < (L+1)*(L+1); k++)
  {
    x = (550.0/double(L))*double((k%(L+1)));
    y = (450.0/double(L))*double(std::floor(k/(L+1)));
    val = std::sin(2.0*M_PI*x/550.0) * std::sin(2.0*M_PI*y/450.0);
    theta_std.Replace_Element(k, val);
    val = std::cos(2.0*M_PI*x/550.0) * std::cos(2.0*M_PI*y/450.0);
    theta_std.Replace_Element((L+1)*(L+1) + k, val);
  }

  int myRank = comm->getRank();
  HDSA::Ptr<std::ostream> outStream;
  if(myRank == 0)
    {
      outStream = HDSA::makePtrFromRef(std::cout);
    }
  else
    {	
      outStream =  HDSA::makePtrFromRef(bhs);
    }
  
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
  HDSA::Ptr<DynamicPDE<RealT> > pde = HDSA::makePtr<DynamicPDE_shallow_ice<RealT> >(*parlist);
  HDSA::Ptr<DynamicPDE_shallow_ice<RealT> > pde_shallow_ice = HDSA::dynamicPtrCast<DynamicPDE_shallow_ice<RealT> >(pde);  
  
  /*************************************************************************/
  /***************** BUILD CONSTRAINT **************************************/
  /*************************************************************************/
  HDSA::Ptr<DynConstraint<RealT> > dyn_con
    = HDSA::makePtr<DynConstraint<RealT> >(pde,meshMgr,comm->Get_Teuchos_Communicator(),*parlist,*outStream);
  dyn_con->getAssembler()->printMeshData(*outStream);
  
  /*************************************************************************/
  /***************** BUILD STATE VECTORS ***********************************/
  /*************************************************************************/
  HDSA::Ptr<Tpetra::MultiVector<> > u0_ptr, uo_ptr, un_ptr,z_ptr, ck_ptr;
  HDSA::Ptr<ROL::Vector<RealT> > u0, uo, un, ck, z;
  u0_ptr = dyn_con->getAssembler()->createStateVector();
  uo_ptr = dyn_con->getAssembler()->createStateVector();
  un_ptr = dyn_con->getAssembler()->createStateVector();
  z_ptr = dyn_con->getAssembler()->createControlVector();
  ck_ptr = dyn_con->getAssembler()->createResidualVector();
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
  std::vector<ROL::TimeStamp<RealT> > timeStamp;
  timeStamp.resize(nt);
  for( int k=0; k<nt; ++k ) {
    timeStamp.at(k).t.resize(2);
    timeStamp.at(k).t.at(0) = k*dt;
    timeStamp.at(k).t.at(1) = (k+1)*dt;
  }

  /*************************************************************************/
  /***************** SET SENSOR CONFIGURATION ******************************/
  /*************************************************************************/
  int nx = parlist->sublist("Geometry").get("NX",70);
  int ny = parlist->sublist("Geometry").get("NY",70);
  int nsx = parlist->sublist("Geometry").get("Sensors Per x-Dimension",10);
  int nsy = parlist->sublist("Geometry").get("Sensors Per y-Dimension",10);
  
  /*** Check sensors and mesh nodes coorespond ***/
  if (nx % (nsx-1) != 0 || ny % (nsy-1) != 0) {
    std::cout << "Error: NX and NY must be divisible by the number of sensors in their respective dimensions \n";
  }
  
  /*** Set sensor locations ***/
  std::vector<int> data_weight_id = std::vector<int>(2*nsx*nsy);
  int count = 0;
  for(int i = 0; i < nsy; i++)
    {
      for(int j = 0; j < nsx; j++)
	{
	  data_weight_id[count] = 3*( (nx/(nsx-1))*j + (ny/(nsy-1))*(nx+1)*i ) + 1;
	  data_weight_id[count+1] = 3*( (nx/(nsx-1))*j + (ny/(nsy-1))*(nx+1)*i ) + 2;
	  count += 2;
	}
    }
  
  /*************************************************************************/
  /***************** READ DATA *********************************************/
  /*************************************************************************/ 
  std::vector<std::vector<RealT> > data;
  data.resize(nt);
  int dim = u0->dimension();
  // read in data
  std::ifstream in("data.txt");          
  count = 0;
  // read the elements in the file into a vector  
  // test file open   
  if (in) {   
    for(int i=0; i < nt; i++)
      {
	data[i].resize(dim);
	for(int j=0; j < dim; j++)
	  {
	    in >> data[i][j];
	  }
      }   
  }
  else
    {
      std::cout << "Error loading the data from data.txt" << std::endl;
    }    
  
  int num_active = data_weight_id.size();
  std::vector<RealT> data_weight = std::vector<RealT>(num_active,1.0);
    
  /*************************************************************************/
  /***************** BUILD COST FUNCTIONAL *********************************/
  /*************************************************************************/
  HDSA::Ptr<Objective_SimOpt_TS<RealT> > misfit_obj = HDSA::makePtr<State_Cost_shallow_ice<RealT> >(timeStamp,data,data_weight,data_weight_id,u0);
  HDSA::Ptr<QoI<RealT> > qoiH1 = HDSA::makePtr<QoI_H1_shallow_ice<RealT> >(pde_shallow_ice->getFE(),pde_shallow_ice->getFieldHelper());
  std::vector<HDSA::Ptr<ROL::Objective_SimOpt<RealT> > > reg_obj(1);  
  bool construct_matrices = parlist->sublist("Problem").get("Construct Matrices",false);
  reg_obj[0] = HDSA::makePtr<Elliptic_Prior_Regularization_Objective<RealT> >(comm->Get_Teuchos_Communicator(), parlist, outStream, construct_matrices);
  std::vector<RealT> weights = std::vector<RealT>(2,1.0);
  weights[0] = parlist->sublist("Problem").get("State Cost",1.0);
  weights[0] = weights[0]*(static_cast<RealT>(nt)/T);
  weights[1] = weights[1]*(static_cast<RealT>(nt)/T)*(1/static_cast<RealT>(nt-1));
  HDSA::Ptr<Objective_SimOpt_TS<RealT> > obj_k = HDSA::makePtr<Misfit_Regularization_Objective_SimOpt_TS<RealT> >(timeStamp, misfit_obj,reg_obj,weights);
  HDSA::Ptr<LTI_Objective_TS<RealT> > dyn_obj = HDSA::makePtr<LTI_Objective_TS<RealT> >(*parlist,obj_k);
  /***************** BUILD REDUCED COST FUNCTIONAL *************************/
  /*************************************************************************/
  ROL::ParameterList &rpl = parlist->sublist("Reduced Dynamic Objective");
  HDSA::Ptr<ROL::Objective<RealT> > robj = HDSA::makePtr<ROL::ReducedDynamicObjective_Stationary_Control<RealT> >(dyn_obj, dyn_con, u0, z, ck, timeStamp, rpl);
  
  // Misfit and regularization objectives for LIS codes
  // Misfit
  std::vector<RealT> weights_misfit = std::vector<RealT>(2);
  weights_misfit[0] = weights[0];
  weights_misfit[1] = 0.0;
  HDSA::Ptr<Objective_SimOpt_TS<RealT> > obj_k_misfit = HDSA::makePtr<Misfit_Regularization_Objective_SimOpt_TS<RealT> >(timeStamp, misfit_obj,reg_obj,weights_misfit);
  HDSA::Ptr<LTI_Objective_TS<RealT> > dyn_obj_misfit = HDSA::makePtr<LTI_Objective_TS<RealT> >(*parlist,obj_k_misfit);
  HDSA::Ptr<ROL::Objective<RealT> > robj_misfit = HDSA::makePtr<ROL::ReducedDynamicObjective_Stationary_Control<RealT> >(dyn_obj_misfit, dyn_con, u0, z, ck, timeStamp, rpl);
  // Regularization
  std::vector<RealT> weights_reg = std::vector<RealT>(1);
  RealT time_scaling = T*static_cast<RealT>(nt-1)/static_cast<RealT>(nt);
  weights_reg[0] = 1.0; //weights[1]*time_scaling;
  HDSA::Ptr<ROL::Objective_SimOpt<RealT> > obj_reg = HDSA::makePtr<ROL::LinearCombinationObjective_SimOpt<RealT> >(weights_reg,reg_obj);
  HDSA::Ptr<ROL::Objective<RealT> > robj_reg = HDSA::makePtr<Reduced_Objective_Regularization<RealT> >(obj_reg, u0); 
  
  /*************************************************************************/
  /***************** RUN VECTOR AND DERIVATIVE CHECKS **********************/
  /*************************************************************************/
  // Run derivative checks
  bool checkDeriv = parlist->sublist("Problem").get("Check Derivatives",true);
  if ( checkDeriv ) 
    {
      // Create state vector
      HDSA::Ptr<Tpetra::MultiVector<> > du_ptr = dyn_con->getAssembler()->createStateVector();
      HDSA::Ptr<ROL::Vector<RealT> > du = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(du_ptr,pde,*dyn_con->getAssembler());
      
      // Create control vector
      HDSA::Ptr<Tpetra::MultiVector<> > dz_ptr = dyn_con->getAssembler()->createControlVector();
      HDSA::Ptr<ROL::Vector<RealT> > dz = HDSA::makePtr<PDE_PrimalOptVector<RealT> >(dz_ptr,pde,*dyn_con->getAssembler());
      
      un->randomize();
      uo->randomize();
      du->randomize();
      z->randomize();
      dz->randomize();
      obj_k->Update_current_TS(timeStamp[1].t.at(0));
      
      *outStream << std::endl << "Check State Gradient of Full Space Objective Function" << std::endl;
      obj_k->checkGradient_1(*uo,*z,*du);
      *outStream << std::endl << "Check Control Gradient of Full Space Objective Function" << std::endl;
      obj_k->checkGradient_2(*uo,*z,*dz);
      
      ROL::ValidateFunction<RealT> validate(1,13,20,11,true,*outStream);
      ROL::DynamicConstraintCheck<RealT>::check(*dyn_con,validate,*uo,*un,*z);       
      
      *outStream << std::endl << "Check Gradient of Reduced Objective Function" << std::endl;
      robj->checkGradient(*z,*dz,true,*outStream);
      *outStream << std::endl << "Check Hessian of Reduced Objective Function" << std::endl;
      robj->checkHessVec(*z,*dz,true,*outStream);
    }
 
  
  // Set initial vector
  z->zero();
  Set_Initial_Iterate<RealT>(z_ptr,parlist);
  dyn_con->outputTpetraVector(z_ptr,"initial_iterate.txt");
	
  // Build optimization problem and check derivatives
  ROL::OptimizationProblem<RealT> optProb(robj,z);
  // Build optimization solver and solve
  ROL::OptimizationSolver<RealT> optSolver(optProb,*parlist);
  std::clock_t timer = std::clock();
  optSolver.solve(*outStream);
  *outStream << "Trust Region Time: "
	      << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC)
	      << " seconds." << std::endl << std::endl;
  
  dyn_con->outputTpetraVector(z_ptr,"optimal_z.txt");
  
  return 0;
}
