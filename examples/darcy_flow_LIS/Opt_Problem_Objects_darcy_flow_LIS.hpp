#ifndef OPT_PROBLEM_OBJECTS_DARCY_FLOW_LIS_HPP
#define OPT_PROBLEM_OBJECTS_DARCY_FLOW_LIS_HPP

// ROL instantiation of Opt_Problem_Objects

template <class RealT>
class Opt_Problem_Objects_darcy_flow_LIS : public Opt_Problem_Objects_ROL<RealT> {

private:

  const HDSA::Ptr<HDSA::ParameterList > parlist;
  HDSA::Ptr<std::ostream> outStream;
  HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > con;
  HDSA::Ptr<Linear_PDE_Constraint<RealT> > pdecon;
  HDSA::Ptr<Assembler<RealT> > assembler;
  HDSA::Ptr<Tpetra::MultiVector<> > u_ptr, p_ptr, r_ptr;
  HDSA::Ptr<std::vector<RealT> > z_ptr;
  HDSA::Ptr<ROL::Vector<RealT> > up, zp, pp, rp;
  std::vector<RealT> target_data;
  std::vector<int> target_data_ids;
  HDSA::Ptr<ROL::Objective_SimOpt<RealT> > obj;


public:

  Opt_Problem_Objects_darcy_flow_LIS(const HDSA::Ptr<HDSA::ParameterList > & parlist_in, const HDSA::Ptr<const HDSA::Comm<int> > & comm): 
    Opt_Problem_Objects_ROL<RealT>(parlist_in,comm), parlist(parlist_in)
  { }

  Opt_Problem_Objects_darcy_flow_LIS(const HDSA::Ptr<HDSA::ParameterList > & parlist_in, const HDSA::Ptr<HDSA::Vector<RealT> > & theta_in, const HDSA::Ptr<const HDSA::Comm<int> > & comm):
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
    HDSA::Ptr<MeshManager<RealT> > meshMgr = HDSA::makePtr<MeshManager_darcy_flow<RealT>>(*parlist);
    // Initialize PDE describing advection-diffusion equation
    HDSA::Ptr<PDE<RealT> > pde = HDSA::makePtr<PDE_darcy_flow<RealT> >(*parlist);
    con = HDSA::makePtr<Linear_PDE_Constraint<RealT> >(pde,meshMgr,comm->Get_Teuchos_Communicator(),*parlist,*outStream,true);
    pdecon = HDSA::dynamicPtrCast<Linear_PDE_Constraint<RealT> >(con);
    assembler = pdecon->getAssembler();
  
    int L = parlist->sublist("Problem").get("Number of Uncertainty Basis Function", 10);
    std::vector<RealT> param = std::vector<RealT>((L+1)*(L+1),0.0);
    pdecon->setParameter(param);
    con->setSolveParameters(*parlist);
  
    /*************************************************************************/
    /***************** BUILD VECTORS *****************************************/
    /*************************************************************************/
    int controlDim = 6;
    // Create state vectors
    u_ptr  = assembler->createStateVector();
    up  = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(u_ptr,pde,assembler,*parlist);
    p_ptr  = assembler->createStateVector();
    pp  = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(p_ptr,pde,assembler,*parlist);
    // Create residual vector
    r_ptr  = assembler->createResidualVector();
    rp  = HDSA::makePtr<PDE_DualSimVector<RealT> >(r_ptr,pde,assembler,*parlist);
    // Create control vectors
    z_ptr  = HDSA::makePtr<std::vector<RealT> >(controlDim,0.0);
    zp = HDSA::makePtr<PDE_OptVector<RealT> >(HDSA::makePtr<ROL::StdVector<RealT> >(z_ptr));

    target_data = std::vector<RealT>(12,0.0);
    target_data_ids = std::vector<int>(12);
  
    int nx = parlist->sublist("Geometry").get("NX", 3);
    int ny = parlist->sublist("Geometry").get("NY", 3);
    target_data_ids[0] = Map_Coords_to_Node(0.28,0.31,nx,ny);
    target_data[0] = 8.2;
    target_data_ids[1] = Map_Coords_to_Node(0.41,0.49,nx,ny);
    target_data[1] = 8.5;
    target_data_ids[2] = Map_Coords_to_Node(0.56,0.66,nx,ny);
    target_data[2] = 8.7;
    target_data_ids[3] = Map_Coords_to_Node(0.7,0.33,nx,ny);
    target_data[3] = 7.5;
    target_data_ids[4] = Map_Coords_to_Node(0.29,0.52,nx,ny);
    target_data[4] = 7.6;
    target_data_ids[5] = Map_Coords_to_Node(0.41,0.64,nx,ny);
    target_data[5] = 8.1;
    target_data_ids[6] = Map_Coords_to_Node(0.54,0.32,nx,ny);
    target_data[6] = 7.3;
    target_data_ids[7] = Map_Coords_to_Node(0.69,0.5,nx,ny);
    target_data[7] = 7.8;
    target_data_ids[8] = Map_Coords_to_Node(0.32,0.62,nx,ny);
    target_data[8] = 8.0;
    target_data_ids[9] = Map_Coords_to_Node(0.43,0.29,nx,ny);
    target_data[9] = 7.7;
    target_data_ids[10] = Map_Coords_to_Node(0.55,0.53,nx,ny);
    target_data[10] = 7.9;
    target_data_ids[11] = Map_Coords_to_Node(0.72,0.64,nx,ny);
    target_data[11] = 8.2;
  
    /*************************************************************************/
    /***************** BUILD COST FUNCTIONAL *********************************/
    /*************************************************************************/ 
    std::vector<HDSA::Ptr<ROL::Objective_SimOpt<RealT> > > obj_vec;
    obj_vec.resize(2);
    obj_vec[0] = HDSA::makePtr<State_Cost_darcy_flow<RealT> >(target_data, target_data_ids, up);
    HDSA::Ptr<QoI<RealT> > qoi_penalty =  HDSA::makePtr<QoI_L2Penalty_darcy_flow<RealT> >(); 
    obj_vec[1] = HDSA::makePtr<IntegralObjective<RealT> >(qoi_penalty,assembler);
    std::vector<RealT> weights = std::vector<RealT>(2);
    weights[0] = parlist->sublist("Problem").get("State Cost",1.0);
    weights[1] = parlist->sublist("Problem").get("Control Cost",1.e-3);
    obj = HDSA::makePtr<ROL::LinearCombinationObjective_SimOpt<RealT> >(weights,obj_vec);

    HDSA::Ptr<ROL::VectorController<RealT> > stateStore = HDSA::makePtr<ROL::VectorController<RealT> >();
    HDSA::Ptr<ROL::Reduced_Objective_SimOpt<RealT> > robj = HDSA::makePtr<ROL::Reduced_Objective_SimOpt<RealT> >(obj, con, stateStore, up, zp, pp, true, false);

    HDSA::Ptr<ROL::Reduced_Objective_SimOpt<RealT> > robj_misfit = HDSA::makePtr<ROL::Reduced_Objective_SimOpt<RealT> >(obj_vec[0], con, stateStore, up, zp, pp, true, false);
    HDSA::Ptr<ROL::Reduced_Objective_SimOpt<RealT> > robj_reg = HDSA::makePtr<ROL::Reduced_Objective_SimOpt<RealT> >(obj_vec[1], con, stateStore, up, zp, pp, true, false);  

    // Run derivative checks
    bool checkDeriv = parlist->sublist("Problem").get("Check Derivatives",false);
    if ( checkDeriv ) {
      
      // Create state vector and set to zeroes
      HDSA::Ptr<Tpetra::MultiVector<> > du_ptr = assembler->createStateVector();     du_ptr->randomize();
      HDSA::Ptr<ROL::Vector<RealT> > dup = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(du_ptr,pde,assembler);
      
      // Create control vectors
      HDSA::Ptr<std::vector<RealT> > dz_ptr  = HDSA::makePtr<std::vector<RealT> >(controlDim,0.0);
      HDSA::Ptr<ROL::Vector<RealT> > dzp = HDSA::makePtr<PDE_OptVector<RealT> >(HDSA::makePtr<ROL::StdVector<RealT> >(dz_ptr));
      
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

      *outStream << "Check Gradient of Full Objective Function" << std::endl;
      obj->checkGradient(x,d,true,*outStream);
      *outStream << std::endl << "Check Hessian of Full Objective Function" << std::endl;
      obj->checkHessVec(x,d,true,*outStream);
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
	
    HDSA::Opt_Problem_Objects<RealT>::rs_obj = HDSA::makePtr<ROL_RS_Objective_LIS<RealT> >(robj,robj_misfit,robj_reg,weights[0],weights[1]);
    HDSA::Opt_Problem_Objects<RealT>::z = HDSA::makePtr<ROL_PDEOPT_Std_Vector<RealT> >(zp->dimension());
    HDSA::Opt_Problem_Objects<RealT>::theta = theta_in;
  }

  virtual ~Opt_Problem_Objects_darcy_flow_LIS()
  { }

  HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > Construct_Opt_Problem_Objects(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<const HDSA::Comm<int> > & comm) 
  { 
    HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects =  HDSA::makePtr<Opt_Problem_Objects_darcy_flow_LIS<RealT> >(parlist,theta,comm);
    return OP_Objects;
  }

  void Load_Optimal_Solution() 
  { 
    // read in solution and write to Opt_Problem_Objects<RealT>::z
    std::ifstream inputFile("control.txt");          
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
	std::cout << "Error loading the optimal control solution" << std::endl;
      }  
  }
  
  void Write_Optimal_Solution()
  {
    RealT tol = 1.e-8;
    up->zero();
    rp->zero();
    zp->set(*dynamic_cast<const ROL_Vector<RealT>&>(*HDSA::Opt_Problem_Objects<RealT>::z).get_rol_vec());
    con->solve(*rp,*up,*zp,tol);
    // Output.
    pdecon->outputTpetraVector(u_ptr,"state.txt");
    pdecon->outputTpetraData();
    assembler->printMeshData(*outStream);
    std::ofstream zfile;
    zfile.open("control.txt");
    for (int i = 0; i < zp->dimension(); i++) 
      {
	zfile << (*z_ptr)[i] << "\n";
      }
    zfile.close();
  
    // Write solutions to text files
    std::string name;
    std::ofstream fout;
    name = "Target_Data.txt";
    fout.open(name);
    for(unsigned i = 0; i < target_data.size(); i++)
      {
	fout << target_data_ids[i] << std::setw(20) << target_data[i] << std::endl;
      }
    fout.close();
    
    zp->zero();
    con->solve(*rp,*up,*zp,tol);
    pdecon->outputTpetraVector(u_ptr,"uncontrolled_state.txt");
  }

private:

int Map_Coords_to_Node(RealT x, RealT y, int nx, int ny)
{
  int x_cell = std::floor(x*static_cast<RealT>(nx));
  int y_cell = std::floor(y*static_cast<RealT>(ny));
  int node = (nx+1)*y_cell + x_cell;
  return node;
}

};


#endif
