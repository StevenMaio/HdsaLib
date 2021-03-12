#ifndef OPT_PROBLEM_OBJECTS_SHALLOW_ICE_HPP
#define OPT_PROBLEM_OBJECTS_SHALLOW_ICE_HPP

// ROL instantiation of Opt_Problem_Objects

template <class RealT>
class Opt_Problem_Objects_shallow_ice : public Opt_Problem_Objects_ROL<RealT> {

private:

  const HDSA::Ptr<HDSA::ParameterList > parlist;
  HDSA::Ptr<std::ostream> outStream;
  HDSA::Ptr<PDE_Constraint<RealT> > con;
  HDSA::Ptr<Tpetra::MultiVector<> > u_ptr, z_ptr, c_ptr;
  HDSA::Ptr<ROL::Vector<RealT> > u, c, z;
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

  
    /*************************************************************************/
    /***************** BUILD GOVERNING PDE ***********************************/
    /*************************************************************************/
    /*** Initialize main data structure. ***/
    HDSA::Ptr<MeshManager<RealT> > meshMgr
      = HDSA::makePtr<MeshManager_shallow_ice<RealT> >(*parlist);
    // Initialize PDE
    HDSA::Ptr<PDE<RealT> > pde = HDSA::makePtr<PDE_shallow_ice<RealT> >(*parlist);
    HDSA::Ptr<PDE_shallow_ice<RealT> > pde_shallow_ice = HDSA::dynamicPtrCast<PDE_shallow_ice<RealT> >(pde); 

    /*************************************************************************/
    /***************** BUILD CONSTRAINT **************************************/
    /*************************************************************************/
    con = HDSA::makePtr<PDE_Constraint<RealT> >(pde,meshMgr,comm->Get_Teuchos_Communicator(),*parlist,*outStream);
  
    /*************************************************************************/
    /***************** BUILD STATE VECTORS ***********************************/
    /*************************************************************************/
    u_ptr = con->getAssembler()->createStateVector();
    z_ptr = con->getAssembler()->createControlVector();
    c_ptr = con->getAssembler()->createResidualVector();
    u = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(u_ptr,pde,con->getAssembler());
    c = HDSA::makePtr<PDE_DualSimVector<RealT> >(c_ptr,pde,con->getAssembler());
    z = HDSA::makePtr<PDE_PrimalOptVector<RealT> >(z_ptr,pde,con->getAssembler());

    con->setSolveParameters(*parlist);
  
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
    data_weight_id = std::vector<int>(nsx*nsy);
    int count = 0;
    for(int i = 0; i < nsy; i++)
      {
	for(int j = 0; j < nsx; j++)
	  {
	    data_weight_id[count] = (nx/(nsx-1))*j + (ny/(nsy-1))*(nx+1)*i;
	    count += 1;
	  }
      }
       
    /*************************************************************************/
    /***************** READ DATA *********************************************/
    /*************************************************************************/ 
    std::vector<RealT> data;
    int dim = u->dimension();
    data.resize(dim);
    // read in data
    std::ifstream in("data.txt");          
    count = 0;
    // read the elements in the file into a vector  
    // test file open   
    if (in) {   
      for(int i = 0; i < dim; i++)
	{
	  in >> data[i];
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

    // Need to construct robj_misfit
    std::vector<HDSA::Ptr<ROL::Objective_SimOpt<RealT> > > obj_vec;
    obj_vec.resize(1);
    obj_vec[0] = HDSA::makePtr<State_Cost_shallow_ice<RealT> >(data, data_weight, data_weight_id, u);
    std::vector<RealT> weights = std::vector<RealT>(1);
    weights[0] = parlist->sublist("Problem").get("State Cost",1.0);
    HDSA::Ptr<ROL::Objective_SimOpt<RealT> > obj = HDSA::makePtr<ROL::LinearCombinationObjective_SimOpt<RealT> >(weights,obj_vec);
    HDSA::Ptr<ROL::SimController<RealT> > stateStore = HDSA::makePtr<ROL::SimController<RealT> >();
    HDSA::Ptr<ROL::Objective<RealT> > robj_misfit = HDSA::makePtr<ROL::Reduced_Objective_SimOpt<RealT> >(obj, con, stateStore, u, z, c, true, false);

    std::vector<HDSA::Ptr<ROL::Objective_SimOpt<RealT> > > reg_obj(2);   
    HDSA::Ptr<QoI<RealT> > qoiH1 = HDSA::makePtr<QoI_H1_shallow_ice<RealT> >(pde_shallow_ice->getFE());
    reg_obj[0] = HDSA::makePtr<IntegralObjective<RealT> >(qoiH1,con->getAssembler());
    HDSA::Ptr<QoI<RealT> > qoiL2 = HDSA::makePtr<QoI_L2_shallow_ice<RealT> >(pde_shallow_ice->getFE());
    reg_obj[1] = HDSA::makePtr<IntegralObjective<RealT> >(qoiL2,con->getAssembler());
    std::vector<RealT> weights_reg = std::vector<RealT>(2);
    weights_reg[0] = parlist->sublist("Problem").get("Gradient Norm Squared Coefficient",1.e-3);
    weights_reg[1] = parlist->sublist("Problem").get("Function Norm Squared Coefficient",1.e-3);
    HDSA::Ptr<ROL::Objective_SimOpt<RealT> > obj_reg = HDSA::makePtr<ROL::LinearCombinationObjective_SimOpt<RealT> >(weights_reg,reg_obj);
    HDSA::Ptr<ROL::Objective<RealT> > robj_reg = HDSA::makePtr<Reduced_Objective_Regularization<RealT> >(obj_reg, u);  
      
    std::vector<HDSA::Ptr<ROL::Objective<RealT> > > total_obj_vec;
    total_obj_vec.resize(2);
    total_obj_vec[0] = robj_misfit;
    total_obj_vec[1] = robj_reg;
    HDSA::Ptr<ROL::Objective<RealT> > robj = HDSA::makePtr<ROL::LinearCombinationObjective<RealT> >(total_obj_vec);

    HDSA::Opt_Problem_Objects<RealT>::rs_obj = HDSA::makePtr<ROL_RS_Objective_LIS<RealT> >(robj,robj_misfit,robj_reg,1.0,1.0);
    HDSA::Opt_Problem_Objects<RealT>::z = HDSA::makePtr<ROL_PDEOPT_Tpetra_Vector<RealT> >(pde,con->getAssembler());
    HDSA::Opt_Problem_Objects<RealT>::theta = theta_in;

    // Set initial vector
    Set_Initial_Iterate(z_ptr,parlist);
    con->outputTpetraVector(z_ptr,"initial_iterate.txt");
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
    std::ifstream inputFile("beta_read.txt");          
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
    z->set(*dynamic_cast<const ROL_Vector<RealT>&>(*HDSA::Opt_Problem_Objects<RealT>::z).get_rol_vec());
    con->outputTpetraVector(z_ptr,"beta.txt");
    con->outputTpetraData();
    con->getAssembler()->printMeshData(*outStream);
    // Output state to file
    RealT tol = 1.e-8;
    con->solve(*c,*u,*z,tol);
    con->outputTpetraVector(u_ptr, "state.txt");
    
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
    RealT val = 1.0;
    return val;
  }


  void Map_Index_to_Coords(int k, std::vector<RealT> & coords, int d, int nx, int ny, RealT nx_float, RealT ny_float) const
  {
    coords[0] = static_cast<RealT>( (k%(nx+1)) )*(1.0/nx_float);
    coords[1] = static_cast<RealT>( std::floor( static_cast<RealT>(k)/(nx_float+1.0) ) )*(1.0/ny_float);
  }

  void Set_Initial_Iterate(HDSA::Ptr<Tpetra::MultiVector<> > & z_ptr, const HDSA::Ptr<HDSA::ParameterList> & parlist) const
  {
    int nx = parlist->sublist("Geometry").get("NX",0);
    int ny = parlist->sublist("Geometry").get("NY",0);
    int d = (nx+1)*(ny+1);
    RealT nx_float = static_cast<RealT>(nx);
    RealT ny_float = static_cast<RealT>(ny);
    std::vector<RealT> coords = std::vector<RealT>(2,0.0);
    for(int k = 0; k < d; k++)
      {
	Map_Index_to_Coords(k,coords,d,nx,ny,nx_float,ny_float);
	RealT val = Initial_Iterate(coords);
	z_ptr->replaceGlobalValue(k,0,val);
      }
  }

};


#endif
