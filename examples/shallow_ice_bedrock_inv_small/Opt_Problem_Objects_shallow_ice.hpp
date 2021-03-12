#ifndef OPT_PROBLEM_OBJECTS_SHALLOW_ICE_HPP
#define OPT_PROBLEM_OBJECTS_SHALLOW_ICE_HPP

// ROL instantiation of Opt_Problem_Objects

template <class RealT>
class Opt_Problem_Objects_shallow_ice : public Opt_Problem_Objects_ROL<RealT> {

private:

  const HDSA::Ptr<HDSA::ParameterList > parlist;
  HDSA::Ptr<std::ostream> outStream;
  HDSA::Ptr<DynConstraint<RealT> > dyn_con;
  HDSA::Ptr<Tpetra::MultiVector<> > u0_ptr, uo_ptr, un_ptr,z_ptr, ck_ptr;
  HDSA::Ptr<ROL::Vector<RealT> > u0, uo, un, ck, z;
  std::vector<ROL::TimeStamp<RealT> > timeStamp;
  std::vector<int> data_weight_id;

public:

  Opt_Problem_Objects_shallow_ice(const HDSA::Ptr<HDSA::ParameterList > & parlist_in, const HDSA::Ptr<const HDSA::Comm<int> > & comm): 
    Opt_Problem_Objects_ROL<RealT>(parlist_in,comm), parlist(parlist_in)
  { }

  Opt_Problem_Objects_shallow_ice(const HDSA::Ptr<HDSA::ParameterList > & parlist_in, const HDSA::Ptr<HDSA::Vector<RealT> > & theta_in, const HDSA::Ptr<const HDSA::Comm<int> > & comm):
    Opt_Problem_Objects_ROL<RealT>(parlist_in,comm), parlist(parlist_in)
  {
    HDSA::nullstream bhs;
    int myRank = comm->getRank();
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
    dyn_con = HDSA::makePtr<DynConstraint<RealT> >(pde,meshMgr,comm->Get_Teuchos_Communicator(),*parlist,*outStream);
    dyn_con->getAssembler()->printMeshData(*outStream);
  
    /*************************************************************************/
    /***************** BUILD STATE VECTORS ***********************************/
    /*************************************************************************/
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
    Set_Initial_Condition(u0_ptr,parlist);

    dyn_con->setSolveParameters(*parlist);
  
    /*************************************************************************/
    /***************** BUILD REDUCED COST FUNCTIONAL *************************/
    /*************************************************************************/
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
    data_weight_id = std::vector<int>(2*nsx*nsy);
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
    // HDSA::Ptr<Objective_SimOpt_TS<RealT> > misfit_obj = HDSA::makePtr<State_Cost_shallow_ice<RealT> >(timeStamp,data,data_weight,data_weight_id,u0);
    // HDSA::Ptr<QoI<RealT> > qoiH1 = HDSA::makePtr<QoI_H1_shallow_ice<RealT> >(pde_shallow_ice->getFE(),pde_shallow_ice->getFieldHelper());
    // std::vector<HDSA::Ptr<ROL::Objective_SimOpt<RealT> > > reg_obj(2);   
    // reg_obj[0] = HDSA::makePtr<IntegralObjective<RealT> >(qoiH1,dyn_con->getAssembler());
    // HDSA::Ptr<QoI<RealT> > qoiL2 = HDSA::makePtr<QoI_L2_shallow_ice<RealT> >(pde_shallow_ice->getFE(),pde_shallow_ice->getFieldHelper());
    // reg_obj[1] = HDSA::makePtr<IntegralObjective<RealT> >(qoiL2,dyn_con->getAssembler());
    // std::vector<RealT> weights = std::vector<RealT>(3);
    // weights[0] = parlist->sublist("Problem").get("State Cost",1.0);
    // weights[1] = parlist->sublist("Problem").get("Gradient Norm Squared Coefficient",1.e-3);
    // weights[2] = parlist->sublist("Problem").get("Function Norm Squared Coefficient",1.e-3);
    // HDSA::Ptr<Objective_SimOpt_TS<RealT> > obj_k = HDSA::makePtr<Misfit_Regularization_Objective_SimOpt_TS<RealT> >(timeStamp, misfit_obj,reg_obj,weights);
    // HDSA::Ptr<LTI_Objective_TS<RealT> > dyn_obj = HDSA::makePtr<LTI_Objective_TS<RealT> >(*parlist,obj_k);
    HDSA::Ptr<Objective_SimOpt_TS<RealT> > misfit_obj = HDSA::makePtr<State_Cost_shallow_ice<RealT> >(timeStamp,data,data_weight,data_weight_id,u0);
    HDSA::Ptr<QoI<RealT> > qoiH1 = HDSA::makePtr<QoI_H1_shallow_ice<RealT> >(pde_shallow_ice->getFE(),pde_shallow_ice->getFieldHelper());
    std::vector<HDSA::Ptr<ROL::Objective_SimOpt<RealT> > > reg_obj(1);   
    reg_obj[0] = HDSA::makePtr<Elliptic_Prior_Regularization_Objective<RealT> >(comm, parlist, outStream);
    std::vector<RealT> weights = std::vector<RealT>(2,1.0);
    weights[0] = parlist->sublist("Problem").get("State Cost",1.0);
    HDSA::Ptr<Objective_SimOpt_TS<RealT> > obj_k = HDSA::makePtr<Misfit_Regularization_Objective_SimOpt_TS<RealT> >(timeStamp, misfit_obj,reg_obj,weights);
    HDSA::Ptr<LTI_Objective_TS<RealT> > dyn_obj = HDSA::makePtr<LTI_Objective_TS<RealT> >(*parlist,obj_k);

    /*************************************************************************/
    /***************** BUILD REDUCED COST FUNCTIONAL *************************/
    /*************************************************************************/
    ROL::ParameterList &rpl = parlist->sublist("Reduced Dynamic Objective");
    HDSA::Ptr<ROL::Objective<RealT> > robj = HDSA::makePtr<ROL::ReducedDynamicObjective_Stationary_Control<RealT> >(dyn_obj, dyn_con, u0, z, ck, timeStamp, rpl);

    // // Misfit and regularization objectives for LIS codes
    // // Misfit
    // std::vector<RealT> weights_misfit = std::vector<RealT>(3);
    // weights_misfit[0] = weights[0];
    // weights_misfit[1] = 0.0;
    // weights_misfit[2] = 0.0;
    // HDSA::Ptr<Objective_SimOpt_TS<RealT> > obj_k_misfit = HDSA::makePtr<Misfit_Regularization_Objective_SimOpt_TS<RealT> >(timeStamp, misfit_obj,reg_obj,weights_misfit);
    // HDSA::Ptr<LTI_Objective_TS<RealT> > dyn_obj_misfit = HDSA::makePtr<LTI_Objective_TS<RealT> >(*parlist,obj_k_misfit);
    // HDSA::Ptr<ROL::Objective<RealT> > robj_misfit = HDSA::makePtr<ROL::ReducedDynamicObjective_Stationary_Control<RealT> >(dyn_obj_misfit, dyn_con, u0, z, ck, timeStamp, rpl);
    // // Regularization
    // std::vector<RealT> weights_reg = std::vector<RealT>(2);
    // RealT time_scaling = T*static_cast<RealT>(nt-1)/static_cast<RealT>(nt);
    // weights_reg[0] = weights[1]*time_scaling;
    // weights_reg[1] = weights[2]*time_scaling;
    // HDSA::Ptr<ROL::Objective_SimOpt<RealT> > obj_reg = HDSA::makePtr<ROL::LinearCombinationObjective_SimOpt<RealT> >(weights_reg,reg_obj);
    // HDSA::Ptr<ROL::Objective<RealT> > robj_reg = HDSA::makePtr<Reduced_Objective_Regularization<RealT> >(obj_reg, u0);  
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
    weights_reg[0] = weights[1]*time_scaling;
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
    
    HDSA::Ptr<ROL_B_Transpose<RealT> > rol_Bt = HDSA::makePtr<ROL_B_Transpose_shallow_ice<RealT> >(parlist,comm,dyn_con->getAssembler(),theta_in);
    HDSA::Opt_Problem_Objects<RealT>::rs_obj = HDSA::makePtr<ROL_RS_Objective_LIS<RealT> >(robj,robj_misfit,robj_reg,1.0,1.0,rol_Bt);
    HDSA::Opt_Problem_Objects<RealT>::z = HDSA::makePtr<ROL_PDEOPT_Tpetra_Vector<RealT> >(pde,dyn_con->getAssembler());
    HDSA::Opt_Problem_Objects<RealT>::theta = theta_in;

    // Set initial vector
    Set_Initial_Iterate(z_ptr,parlist);
    dyn_con->outputTpetraVector(z_ptr,"initial_iterate.txt");
    dynamic_cast<const ROL_Vector<RealT>&>(*HDSA::Opt_Problem_Objects<RealT>::z).get_rol_vec()->set(*z);
  }

  virtual ~Opt_Problem_Objects_shallow_ice()
  { }

  HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > Construct_Opt_Problem_Objects(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<const HDSA::Comm<int> > & comm) 
  { 
    HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects =  HDSA::makePtr<Opt_Problem_Objects_shallow_ice<RealT> >(parlist,theta,comm);
    return OP_Objects;
  }

  void Load_Optimal_Solution() 
  { 
    // read in solution and write to Opt_Problem_Objects<RealT>::z
    std::ifstream inputFile("bed_read.txt");          
    RealT value;
    int count = 0;
    // read the elements in the file into a vector  
    // test file open   
    if (inputFile) {   
      while ( inputFile >> value ) {
	HDSA::Opt_Problem_Objects<RealT>::z->Replace_Element(count,value);
	count += 1;
      }
    }
    else
      {
	std::cout << "Error loading the optimal solution" << std::endl;
      } 
  }
  
  void Write_Optimal_Solution()
  {
    int nt = parlist->sublist("Time Discretization").get("Number of Time Steps", 100);
    z->set(*dynamic_cast<const ROL_Vector<RealT>&>(*HDSA::Opt_Problem_Objects<RealT>::z).get_rol_vec());
    dyn_con->outputTpetraVector(z_ptr,"bed.txt");
    dyn_con->outputTpetraData();
    dyn_con->getAssembler()->printMeshData(*outStream);
    // Output state to file
    uo->set(*u0); un->zero();
    for (int k = 1; k < nt; ++k) {
      // Print previous state to file
      std::stringstream ufile;
      ufile << "state_" << k-1 << ".txt";
      dyn_con->outputTpetraVector(uo_ptr, ufile.str());
      // Advance time stepper
      dyn_con->solve(*ck, *uo, *un, *z, timeStamp[k-1]);
      uo->set(*un);
    }
    // Print previous state to file
    std::stringstream ufile;
    ufile << "state_" << nt-1 << ".txt";
    dyn_con->outputTpetraVector(uo_ptr, ufile.str());
    
    int num_active = data_weight_id.size();
    // Write weights to text files
    std::string name = "Active_Sensors.txt";
    std::ofstream fout;
    fout.open(name);
    for(int i = 0; i < num_active; i++)
      {
	fout << data_weight_id[i] << std::endl;
      }
    fout.close();
  }

 RealT Initial_Iterate(std::vector<RealT> & coords) const
  { 
    RealT val = -17.0;
    return val;
  }


  void Map_Index_to_Coords(int k, std::vector<RealT> & coords, int d, int nx, int ny, RealT nx_float, RealT ny_float) const
  {
    coords[0] = static_cast<RealT>( (k%(nx+1)) )*(1.0/nx_float);
    coords[1] = static_cast<RealT>( std::floor( static_cast<RealT>(k)/(nx_float+1.0) ) )*(1.0/ny_float);
  }

  void Set_Initial_Iterate(HDSA::Ptr<Tpetra::MultiVector<> > & z_ptr, const HDSA::Ptr<HDSA::ParameterList> & parlist) const
  {
    int num_coeff_load = parlist->sublist("Problem").get("Number of Coefficients in Loaded Fields", 10); 
    std::vector<RealT> initial_iter_coeff = std::vector<RealT>(num_coeff_load);
    // read in data
    std::ifstream in("matlab_initial_iterate.txt");          
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
	std::cout << "Error loading the data from matlab_initial_iterate.txt" << std::endl;
      }  

    for(int k = 0; k < num_coeff_load; k++)
      {
	z_ptr->replaceGlobalValue(3*k,0,initial_iter_coeff[k]);
      }

    // int nx = parlist->sublist("Geometry").get("NX",0);
    // int ny = parlist->sublist("Geometry").get("NY",0);
    // int d = (nx+1)*(ny+1);
    // RealT nx_float = static_cast<RealT>(nx);
    // RealT ny_float = static_cast<RealT>(ny);
    // std::vector<RealT> coords = std::vector<RealT>(2,0.0);
    // for(int k = 0; k < d; k++)
    //   {
    // 	Map_Index_to_Coords(k,coords,d,nx,ny,nx_float,ny_float);
    // 	RealT val = Initial_Iterate(coords);
    // 	z_ptr->replaceGlobalValue(3*k,0,val);
    //   }
  }

  void Set_Initial_Condition(HDSA::Ptr<Tpetra::MultiVector<> > & u_ptr, const HDSA::Ptr<HDSA::ParameterList> & parlist) const
  {
    int num_coeff_load = parlist->sublist("Problem").get("Number of Coefficients in Loaded Fields", 10); 

    std::vector<RealT> initial_iter_coeff = std::vector<RealT>(num_coeff_load);
    // read in data
    std::ifstream in("Surface_Height.txt");          
    //std::ifstream in("Synthetic_Surface_Height.txt"); 
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

};


#endif
