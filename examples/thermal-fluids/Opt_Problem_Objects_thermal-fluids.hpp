#ifndef OPT_PROBLEM_OBJECTS_THERMAL_FLUIDS_HPP
#define OPT_PROBLEM_OBJECTS_THERMAL_FLUIDS_HPP

// ROL instantiation of Opt_Problem_Objects

template <class RealT>
class Opt_Problem_Objects_thermal_fluids : public Opt_Problem_Objects_ROL<RealT> {

private:

  const HDSA::Ptr<HDSA::ParameterList > parlist;
  const HDSA::Ptr<const HDSA::Comm<int> > comm_;
  HDSA::Ptr<std::ostream> outStream;
  HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > con;
  HDSA::Ptr<PDE_Constraint<RealT> > pdecon;
  HDSA::Ptr<Assembler<RealT> > assembler;
  HDSA::Ptr<Tpetra::MultiVector<> > u_ptr, z_ptr, p_ptr, r_ptr;
  HDSA::Ptr<ROL::Vector<RealT> > up, zp, pp, rp;

public:

  Opt_Problem_Objects_thermal_fluids(const HDSA::Ptr<HDSA::ParameterList > & parlist_in, const HDSA::Ptr<const HDSA::Comm<int> > & comm): 
    Opt_Problem_Objects_ROL<RealT>(parlist_in,comm), parlist(parlist_in), comm_(comm)
  { }

  Opt_Problem_Objects_thermal_fluids(const HDSA::Ptr<HDSA::ParameterList > & parlist_in, const HDSA::Ptr<HDSA::Vector<RealT> > & theta_in, const HDSA::Ptr<const HDSA::Comm<int> > & comm):
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

    /*** Initialize main data structure. ***/
    HDSA::Ptr<MeshManager<RealT> > meshMgr = HDSA::makePtr<MeshManager_ThermalFluids<RealT> >(*parlist);
    // Initialize PDE describing Navier-Stokes equations.
    HDSA::Ptr<PDE<RealT> > pde = HDSA::makePtr<PDE_ThermalFluids<RealT> >(*parlist);
    HDSA::Ptr<PDE_ThermalFluids<RealT> > pde_tf = HDSA::dynamicPtrCast<PDE_ThermalFluids<RealT> >(pde);
    con = HDSA::makePtr<PDE_Constraint<RealT> >(pde,meshMgr,comm->Get_Teuchos_Communicator(),*parlist,*outStream);
    // Cast the constraint and get the assembler.
    pdecon = HDSA::dynamicPtrCast<PDE_Constraint<RealT> >(con);
    assembler = pdecon->getAssembler();
    con->setSolveParameters(*parlist);
    pdecon->outputTpetraData();

    // Create state vector and set to zeroes
    u_ptr  = assembler->createStateVector();  
    p_ptr  = assembler->createStateVector();     
    r_ptr  = assembler->createResidualVector(); 
    z_ptr  = assembler->createControlVector();  
    up  = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(u_ptr,pde,assembler);
    pp  = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(p_ptr,pde,assembler);
    rp  = HDSA::makePtr<PDE_DualSimVector<RealT> >(r_ptr,pde,assembler);
    zp  = HDSA::makePtr<PDE_PrimalOptVector<RealT> >(z_ptr,pde,assembler);

    // Initialize objective function.
    std::vector<HDSA::Ptr<QoI<RealT> > > qoi_vec(2,ROL::nullPtr);
    qoi_vec[0] = HDSA::makePtr<QoI_State_ThermalFluids<RealT> >(*parlist,
								pde_tf->getVelocityFE(),
								pde_tf->getPressureFE(),
								pde_tf->getThermalFE(),
								pde_tf->getFieldHelper());
    qoi_vec[1] = HDSA::makePtr<QoI_L2Penalty_ThermalFluids<RealT> >(pde_tf->getVelocityFE(),
								    pde_tf->getPressureFE(),
								    pde_tf->getThermalFE(),
								    pde_tf->getThermalBdryFE(),
								    pde_tf->getBdryCellLocIds(),
								    pde_tf->getFieldHelper());
    HDSA::Ptr<StdObjective_ThermalFluids<RealT> > std_obj = HDSA::makePtr<StdObjective_ThermalFluids<RealT> >(*parlist);
    HDSA::Ptr<ROL::Objective_SimOpt<RealT> > obj = HDSA::makePtr<PDE_Objective<RealT> >(qoi_vec,std_obj,assembler);
    HDSA::Ptr<ROL::SimController<RealT> > stateStore = HDSA::makePtr<ROL::SimController<RealT> >();
    HDSA::Ptr<ROL::Reduced_Objective_SimOpt<RealT> > robj = HDSA::makePtr<ROL::Reduced_Objective_SimOpt<RealT> >(obj, con, stateStore, up, zp, pp, true, false);

    RealT tol(1.e-8);
    bool initSolve = parlist->sublist("Problem").get("Solve state for full space",true);
    if (initSolve) 
      {
	con->solve(*rp,*up,*zp,tol);
	pdecon->outputTpetraVector(u_ptr,"state_uncontrolled.txt");
      }

    // Initialize reduced objective function
    HDSA::Opt_Problem_Objects<RealT>::fs_obj = HDSA::makePtr<ROL_FS_Objective_Model_Error<RealT> >(obj);
    HDSA::Opt_Problem_Objects<RealT>::con = HDSA::makePtr<ROL_Constraint<RealT> >(con);
    HDSA::Opt_Problem_Objects<RealT>::u = HDSA::makePtr<ROL_PDEOPT_Tpetra_Vector<RealT> >(pde,assembler);
    HDSA::Opt_Problem_Objects<RealT>::rs_obj = HDSA::makePtr<ROL_RS_Objective_Model_Error<RealT> >(robj);
    HDSA::Opt_Problem_Objects<RealT>::z = HDSA::makePtr<ROL_PDEOPT_Tpetra_Vector<RealT> >(pde,assembler);
    HDSA::Opt_Problem_Objects<RealT>::lambda = HDSA::makePtr<ROL_PDEOPT_Tpetra_Vector<RealT> >(pde,assembler);
    HDSA::Opt_Problem_Objects<RealT>::theta = theta_in;
  }

  virtual ~Opt_Problem_Objects_thermal_fluids()
  { }

  HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > Construct_Opt_Problem_Objects(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<const HDSA::Comm<int> > & comm) 
  { 
    HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects =  HDSA::makePtr<Opt_Problem_Objects_thermal_fluids<RealT> >(parlist,theta,comm);
    return OP_Objects;
  }

  void Load_Optimal_Solution() 
  { 
    // read in solution and write to Opt_Problem_Objects<RealT>::z
    std::ifstream inputFile("control_read.txt");          
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
	std::cout << "Error loading the optimal source solution" << std::endl;
      }   
    RealT tol = 1.e-8;
    up->zero();
    rp->zero();
    zp->set(*dynamic_cast<const ROL_Vector<RealT>&>(*HDSA::Opt_Problem_Objects<RealT>::z).get_rol_vec());
    con->solve(*rp,*dynamic_cast<const ROL_Vector<RealT>&>(*HDSA::Opt_Problem_Objects<RealT>::u).get_rol_vec(),*zp,tol); 
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
    pdecon->outputTpetraVector(z_ptr,"control.txt");
    pdecon->outputTpetraData();
    assembler->printMeshData(*outStream);
  }

  const HDSA::Ptr<const Teuchos::Comm<int> > Get_comm(void) const
  {
    return comm_->Get_Teuchos_Communicator();
  }

};


#endif
