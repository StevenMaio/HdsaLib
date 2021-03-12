#ifndef OPT_PROBLEM_OBJECTS_DARCY_ADV_DIFF_HPP
#define OPT_PROBLEM_OBJECTS_DARCY_ADV_DIFF_HPP

// ROL instantiation of Opt_Problem_Objects

template <class RealT>
class Opt_Problem_Objects_darcy_adv_diff : public Opt_Problem_Objects_ROL<RealT> {

private:

  const HDSA::Ptr<HDSA::ParameterList > parlist;
  HDSA::Ptr<std::ostream> outStream;
  HDSA::Ptr<DynConstraint<RealT> > dyn_con;
  HDSA::Ptr<Tpetra::MultiVector<> > u0_ptr, uo_ptr, un_ptr,z_ptr, ck_ptr;
  HDSA::Ptr<ROL::Vector<RealT> > u0, uo, un, ck, z;
  std::vector<ROL::TimeStamp<RealT> > timeStamp;
  std::vector<int> data_weight_id;

public:

  Opt_Problem_Objects_darcy_adv_diff(const HDSA::Ptr<HDSA::ParameterList > & parlist_in, const HDSA::Ptr<const HDSA::Comm<int> > & comm): 
    Opt_Problem_Objects_ROL<RealT>(parlist_in,comm), parlist(parlist_in)
  { }

  Opt_Problem_Objects_darcy_adv_diff(const HDSA::Ptr<HDSA::ParameterList > & parlist_in, const HDSA::Ptr<HDSA::Vector<RealT> > & theta_in, const HDSA::Ptr<const HDSA::Comm<int> > & comm):
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
      = HDSA::makePtr<MeshManager_darcy_adv_diff<RealT> >(*parlist);
    // Initialize PDE describing advection-diffusion equation
    HDSA::Ptr<DynamicPDE<RealT> > pde = HDSA::makePtr<DynamicPDE_darcy_adv_diff<RealT> >(*parlist);
    HDSA::Ptr<DynamicPDE_darcy_adv_diff<RealT> > pde_darcy_adv_diff = HDSA::dynamicPtrCast<DynamicPDE_darcy_adv_diff<RealT> >(pde);  

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
    int nsp, nsc; 
    int nsx = parlist->sublist("Geometry").get("Sensors Per x-Dimension",10);
    int nsy = parlist->sublist("Geometry").get("Sensors Per y-Dimension",10);
    nsp = (nsx-2)*nsy;
    nsc = nsx*nsy;
    
    /*** Check sensors and mesh nodes coorespond ***/
    if (nx % (nsx-1) != 0 || ny % (nsy-1) != 0) {
      std::cout << "Error: NX and NY must be divisible by the number of sensors in their respective dimensions \n";
    }
    
    /*** Set sensor locations ***/
    std::vector<int> p_ids = std::vector<int>((nsx-2)*nsy);
    int count = 0;
    for(int i = 0; i < nsy; i++)
      {
	for(int j = 1; j < nsx-1; j++)
	  {
	    p_ids[count] = (nx/(nsx-1))*j + (ny/(nsy-1))*(nx+1)*i;
	    count += 1;
	  }
      }
    
    std::vector<int> u_ids = std::vector<int>(nsx*nsy);
    count = 0;
    for(int i = 0; i < nsy; i++)
      {
	for(int j = 0; j < nsx; j++)
	  {
	    u_ids[count] = (nx/(nsx-1))*j + (ny/(nsy-1))*(nx+1)*i;
	    count += 1;
	  }
      }
    
    int n = (nsx-2)*nsy + nsx*nsy;
    data_weight_id.resize(n);
    for(int i = 0; i < (nsx-2)*nsy; i++)
      {
	data_weight_id[i] = 2*p_ids[i];
      }
    for(int i = 0; i < nsx*nsy; i++)
      {
	data_weight_id[i+(nsx-2)*nsy] = 2*u_ids[i]+1;
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
    
    std::vector<RealT> p_data = std::vector<RealT>(nt*nsp);
    std::vector<RealT> c_data = std::vector<RealT>(nt*nsc);
    int countp = 0;
    int countc = 0;
    for(int k = 0; k < nt; k++)
      {
	for(int i = 0; i<nsp; i++)
	  {
	    int l1 = data_weight_id[i];
	    p_data[countp] = data[k][l1];
	    countp += 1;  
	  }
	for(int i = 0; i<nsc; i++)
	  {
	    int l2 = data_weight_id[nsp+i];
	    c_data[countc] = data[k][l2];
	    countc += 1;  
	  }
      }
    
    RealT av_p_data = 0;
    RealT av_c_data = 0;
    for(int k=0; k<nt*nsp; k++)
      {
	av_p_data += p_data[k];
      }
    for(int k=0; k<nt*nsc; k++)
      {
	av_c_data += c_data[k];
      }
    av_p_data = av_p_data/static_cast<RealT>(nt*nsp);
    av_c_data = av_c_data/static_cast<RealT>(nt*nsc);
    
    int num_active = data_weight_id.size();
    std::vector<RealT> data_weight = std::vector<RealT>(num_active,1.0);
    for(int k=0; k<nsp; k++)
      {
	data_weight[k] = data_weight[k]/(av_p_data*av_p_data*static_cast<RealT>(nt*nsp));
      }
    for(int k=0; k<nsc; k++)
      {
	data_weight[nsp+k] = 2.0*data_weight[nsp+k]/(av_c_data*av_c_data*static_cast<RealT>(nt*nsc));
      }
    
    /*************************************************************************/
    /***************** BUILD COST FUNCTIONAL *********************************/
    /*************************************************************************/
    HDSA::Ptr<Objective_SimOpt_TS<RealT> > misfit_obj = HDSA::makePtr<State_Cost_darcy_adv_diff<RealT> >(timeStamp,data,data_weight,data_weight_id,u0);
    HDSA::Ptr<QoI<RealT> > qoiH1 = HDSA::makePtr<QoI_H1_darcy_adv_diff<RealT> >(pde_darcy_adv_diff->getPressureFE(),pde_darcy_adv_diff->getFieldHelper());
    std::vector<HDSA::Ptr<ROL::Objective_SimOpt<RealT> > > reg_obj(3);   
    reg_obj[0] = HDSA::makePtr<IntegralObjective<RealT> >(qoiH1,dyn_con->getAssembler());
    HDSA::Ptr<QoI<RealT> > qoiL2 = HDSA::makePtr<QoI_L2_darcy_adv_diff<RealT> >(pde_darcy_adv_diff->getPressureFE(),pde_darcy_adv_diff->getFieldHelper());
    reg_obj[1] = HDSA::makePtr<IntegralObjective<RealT> >(qoiL2,dyn_con->getAssembler());
    HDSA::Ptr<QoI<RealT> > qoiInformed = HDSA::makePtr<QoI_Informed_darcy_adv_diff<RealT> >(pde_darcy_adv_diff->getPressureFE(),pde_darcy_adv_diff->getFieldHelper());
    reg_obj[2] = HDSA::makePtr<IntegralObjective<RealT> >(qoiInformed,dyn_con->getAssembler());
    std::vector<RealT> weights = std::vector<RealT>(4);
    weights[0] = parlist->sublist("Problem").get("State Cost",1.0);
    weights[1] = parlist->sublist("Problem").get("Gradient Norm Squared Coefficient",1.e-3);
    weights[2] = parlist->sublist("Problem").get("Function Norm Squared Coefficient",1.e-3);
    weights[3] = parlist->sublist("Problem").get("Informed Regularization Coefficient",1.e-3);
    HDSA::Ptr<Objective_SimOpt_TS<RealT> > obj_k = HDSA::makePtr<Misfit_Regularization_Objective_SimOpt_TS<RealT> >(timeStamp, misfit_obj,reg_obj,weights);
    HDSA::Ptr<LTI_Objective_TS<RealT> > dyn_obj = HDSA::makePtr<LTI_Objective_TS<RealT> >(*parlist,obj_k);

    /*************************************************************************/
    /***************** BUILD REDUCED COST FUNCTIONAL *************************/
    /*************************************************************************/
    ROL::ParameterList &rpl = parlist->sublist("Reduced Dynamic Objective");
    HDSA::Ptr<ROL::Objective<RealT> > robj = HDSA::makePtr<ROL::ReducedDynamicObjective_Stationary_Control<RealT> >(dyn_obj, dyn_con, u0, z, ck, timeStamp, rpl);

    // Misfit and regularization objectives for LIS codes
    // Misfit
    std::vector<RealT> weights_misfit = std::vector<RealT>(4);
    weights_misfit[0] = weights[0];
    weights_misfit[1] = 0.0;
    weights_misfit[2] = 0.0;
    weights_misfit[3] = 0.0;
    HDSA::Ptr<Objective_SimOpt_TS<RealT> > obj_k_misfit = HDSA::makePtr<Misfit_Regularization_Objective_SimOpt_TS<RealT> >(timeStamp, misfit_obj,reg_obj,weights_misfit);
    HDSA::Ptr<LTI_Objective_TS<RealT> > dyn_obj_misfit = HDSA::makePtr<LTI_Objective_TS<RealT> >(*parlist,obj_k_misfit);
    HDSA::Ptr<ROL::Objective<RealT> > robj_misfit = HDSA::makePtr<ROL::ReducedDynamicObjective_Stationary_Control<RealT> >(dyn_obj_misfit, dyn_con, u0, z, ck, timeStamp, rpl);
    // Regularization
    std::vector<RealT> weights_reg = std::vector<RealT>(3);
    RealT time_scaling = T*static_cast<RealT>(nt-1)/static_cast<RealT>(nt);
    weights_reg[0] = weights[1]*time_scaling;
    weights_reg[1] = weights[2]*time_scaling;
    weights_reg[2] = weights[3]*time_scaling;
    HDSA::Ptr<ROL::Objective_SimOpt<RealT> > obj_reg = HDSA::makePtr<ROL::LinearCombinationObjective_SimOpt<RealT> >(weights_reg,reg_obj);
    HDSA::Ptr<ROL::Objective<RealT> > robj_reg = HDSA::makePtr<Reduced_Objective_Regularization<RealT> >(obj_reg, u0);  
  
    /*************************************************************************/
    /***************** RUN VECTOR AND DERIVATIVE CHECKS **********************/
    /*************************************************************************/
    // Run derivative checks
    bool checkDeriv = parlist->sublist("Problem").get("Check Derivatives",true);
    if ( checkDeriv ) 
      {
	// Create control vector
	HDSA::Ptr<Tpetra::MultiVector<> > dz_ptr = dyn_con->getAssembler()->createControlVector();
	HDSA::Ptr<ROL::Vector<RealT> > dz = HDSA::makePtr<PDE_PrimalOptVector<RealT> >(dz_ptr,pde,*dyn_con->getAssembler());
	
	z->randomize();
	dz->randomize();
	
	for(int i = 1; i < z->dimension(); i++)
	  {
	    z_ptr->replaceGlobalValue(i,0,0.0);
	    dz_ptr->replaceGlobalValue(i,0,0.0);
	  }
	
	*outStream << std::endl << "Check Gradient of Reduced Objective Function" << std::endl;
	robj->checkGradient(*z,*dz,true,*outStream);
	*outStream << std::endl << "Check Hessian of Reduced Objective Function" << std::endl;
	robj->checkHessVec(*z,*dz,true,*outStream);
      }
    
    //HDSA::Opt_Problem_Objects<RealT>::rs_obj = HDSA::makePtr<ROL_RS_Objective<RealT> >(robj);
    HDSA::Opt_Problem_Objects<RealT>::rs_obj = HDSA::makePtr<ROL_RS_Objective_LIS<RealT> >(robj,robj_misfit,robj_reg,1.0,1.0);
    HDSA::Opt_Problem_Objects<RealT>::z = HDSA::makePtr<ROL_PDEOPT_Tpetra_Vector<RealT> >(pde,dyn_con->getAssembler());
    HDSA::Opt_Problem_Objects<RealT>::theta = theta_in;

    // Set initial vector
    RealT init_iter_noise = parlist->sublist("Problem").get("Initial Iterate Noise",0.0);
    Set_Initial_Iterate(z_ptr,parlist,init_iter_noise);
    dyn_con->outputTpetraVector(z_ptr,"initial_iterate.txt");
    dynamic_cast<const ROL_Vector<RealT>&>(*HDSA::Opt_Problem_Objects<RealT>::z).get_rol_vec()->set(*z);

  }

  virtual ~Opt_Problem_Objects_darcy_adv_diff()
  { }

  HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > Construct_Opt_Problem_Objects(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<const HDSA::Comm<int> > & comm) 
  { 
    HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects =  HDSA::makePtr<Opt_Problem_Objects_darcy_adv_diff<RealT> >(parlist,theta,comm);
    return OP_Objects;
  }

  void Load_Optimal_Solution() 
  { 
    // read in solution and write to Opt_Problem_Objects<RealT>::z
    std::ifstream inputFile("perm_read.txt");          
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
    dyn_con->outputTpetraVector(z_ptr,"perm.txt");
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
    RealT val = -8.0*coords[1]*coords[1] + 8.0*coords[1] - 1.0;
    return val;
  }

  RealT Random_Field_Noise(std::vector<RealT> & coords,std::vector<RealT> & coeff) const
  { 
    RealT pi(M_PI);
    RealT val = 0.0;
    int L = std::sqrt(coeff.size())-1;
    int count = 0;
    RealT x_fun = 0.0;
    RealT y_fun = 0.0;
    //RealT dist_x = 0.0;
    //RealT dist_y = 0.0;
    for(int j = 0; j < L+1; j++)
      {
	y_fun = 0.0;
	// dist_y = std::abs(coords[1]-static_cast<RealT>(j)/static_cast<RealT>(L));
	// if( dist_y < 1.0/static_cast<RealT>(L) )
	//   {
	//     y_fun = 1.0-static_cast<RealT>(L)*dist_y;
	//   }
	y_fun = std::cos(2.0*pi*static_cast<RealT>(j+1)*coords[1]);
	for(int i = 0; i < L+1; i++)
	  {
	    x_fun = 0.0;
	    // dist_x = std::abs(coords[0]-static_cast<RealT>(i)/static_cast<RealT>(L));
	    // if( dist_x < 1.0/static_cast<RealT>(L) )
	    //   {
	    // 	x_fun = 1.0-static_cast<RealT>(L)*dist_x;
	    //   }
	    x_fun = std::cos(2.0*pi*static_cast<RealT>(i+1)*coords[0]);
	    
	    val += coeff[count]*x_fun*y_fun/static_cast<RealT>(j+1);
	    count += 1;
	  }
      }
    return val;
  }
  
  RealT Permeability_Eval(std::vector<RealT> & coords) const
  { 
    int T = 5;
    std::vector<RealT> amp = std::vector<RealT>(T);
    std::vector<RealT> sdx = std::vector<RealT>(T);
    std::vector<RealT> sdy = std::vector<RealT>(T);
    std::vector<RealT> c = std::vector<RealT>(T);
    std::vector<RealT> h = std::vector<RealT>(T);
    std::vector<RealT> k= std::vector<RealT>(T);
    
    amp[0] = -2.0; amp[1] = -2.2; amp[2] = 2.1; amp[3] = 2.0; amp[4] = 1.9;
    sdx[0] = 1.0; sdx[1] = 2.0; sdx[2] = 17.0; sdx[3] = 7.0; sdx[4] = 29.0;
    sdy[0] = 50.0; sdy[1] = 43.0; sdy[2] = 34.0; sdy[3] = 46.0; sdy[4] = 58.0;
    c[0] = 10.0; c[1] = 6.0; c[2] = 8.0; c[3] = 6.0; c[4] = 11.0;
    h[0] = .6; h[1] = .45; h[2] = .27; h[3] = .53; h[4] = .78;
    k[0] = .75; k[1] = .15; k[2] = .43; k[3] = .55; k[4] = .62; 
    
    RealT val = 0.0;
    for(int i = 0; i < T; i++)
      {
	val += amp[i]*std::exp(-1.0*( sdx[i]*(coords[0]-h[i])*(coords[0]-h[i]) + c[i]*(coords[0]-h[i])*(coords[1]-k[i]) + sdy[i]*(coords[1]-k[i])*(coords[1]-k[i]) ) );
      }
    
    ////////////////////////
    val = val/1.2;
    //////////////////////

    return val;
  }

  void Map_Index_to_Coords(int k, std::vector<RealT> & coords, int d, int nx, int ny, RealT nx_float, RealT ny_float) const
  {
    coords[0] = static_cast<RealT>( (k%(nx+1)) )*(1.0/nx_float);
    coords[1] = static_cast<RealT>( std::floor( static_cast<RealT>(k)/(nx_float+1.0) ) )*(1.0/ny_float);
  }

  void Set_Initial_Iterate(HDSA::Ptr<Tpetra::MultiVector<> > & z_ptr, const HDSA::Ptr<HDSA::ParameterList> & parlist, RealT & init_iterate_noise) const
  {
    int L = 5;
    std::vector<RealT> coeff = std::vector<RealT>((L+1)*(L+1));
    
    unsigned seed = 13420958;
    std::default_random_engine generator;
    generator.seed(seed);
    std::normal_distribution<RealT> distribution = std::normal_distribution<RealT>(0.0,1.0);
    
    for(int k = 0; k < (L+1)*(L+1); k++)
      {
	coeff[k] = distribution(generator);
      }
    
    int nx = parlist->sublist("Geometry").get("NX",0);
    int ny = parlist->sublist("Geometry").get("NY",0);
    int d = (nx+1)*(ny+1);
    RealT nx_float = static_cast<RealT>(nx);
    RealT ny_float = static_cast<RealT>(ny);
    std::vector<RealT> coords = std::vector<RealT>(2,0.0);
    RealT m = 0.0;
    for(int k = 0; k < d; k++)
      {
	Map_Index_to_Coords(k,coords,d,nx,ny,nx_float,ny_float);
	//RealT val = Permeability_Eval(coords)+init_iterate_noise*Random_Field_Noise(coords,coeff);
	RealT val = Initial_Iterate(coords);
	z_ptr->replaceGlobalValue(2*k,0,val);
	z_ptr->replaceGlobalValue(2*k+1,0,0.0);
	m = std::max(m,std::abs(val));
      }
    //z_ptr->scale(2.0/m);
    //z_ptr->scale(1.5/m);
  }

};


#endif
