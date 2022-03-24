#ifndef OPT_PROBLEM_OBJECTS_STOKES_HPP
#define OPT_PROBLEM_OBJECTS_STOKES_HPP

// ROL instantiation of Opt_Problem_Objects

template <class RealT>
class Opt_Problem_Objects_stokes : public Opt_Problem_Objects_ROL<RealT> {

private:

  const HDSA::Ptr<HDSA::ParameterList > parlist;
  HDSA::Ptr<std::ostream> outStream;
  HDSA::Ptr<ROL::Objective_SimOpt<RealT> > obj;
  HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > con;
  HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > con_ns;
  HDSA::Ptr<Linear_PDE_Constraint<RealT> > pdecon;
  HDSA::Ptr<Assembler<RealT> > assembler;
  HDSA::Ptr<Tpetra::MultiVector<> > u_ptr, z_ptr, p_ptr, r_ptr;
  HDSA::Ptr<ROL::Vector<RealT> > up, zp, pp, rp;

public:

  Opt_Problem_Objects_stokes(const HDSA::Ptr<HDSA::ParameterList > & parlist_in, const HDSA::Ptr<const HDSA::Comm<int> > & comm): 
    Opt_Problem_Objects_ROL<RealT>(parlist_in,comm), parlist(parlist_in)
  { }

  Opt_Problem_Objects_stokes(const HDSA::Ptr<HDSA::ParameterList > & parlist_in, const HDSA::Ptr<HDSA::Vector<RealT> > & theta_in, const HDSA::Ptr<const HDSA::Comm<int> > & comm):
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
    HDSA::Ptr<MeshManager<RealT> > meshMgr = HDSA::makePtr<MeshManager_Stokes<RealT> >(*parlist);
    // Initialize PDE describing Stokes equations.
    HDSA::Ptr<PDE<RealT> > pde = HDSA::makePtr<PDE_Stokes<RealT> >(*parlist,0.0);
    HDSA::Ptr<PDE_Stokes<RealT> > pde_stokes = HDSA::dynamicPtrCast<PDE_Stokes<RealT> >(pde);
    con = HDSA::makePtr<Linear_PDE_Constraint<RealT> >(pde,meshMgr,comm->Get_Teuchos_Communicator(),*parlist,*outStream);
    // Cast the constraint and get the assembler.
    pdecon = ROL::dynamicPtrCast<Linear_PDE_Constraint<RealT> >(con);
    assembler = pdecon->getAssembler();
    con->setSolveParameters(*parlist);
    
    // Create state vector and set to zeroes
    u_ptr = assembler->createStateVector();    u_ptr->randomize();
    p_ptr = assembler->createStateVector();    p_ptr->randomize();
    z_ptr = assembler->createControlVector();  z_ptr->randomize();
    r_ptr = assembler->createResidualVector(); r_ptr->putScalar(0.0);
    up = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(u_ptr,pde,assembler,*parlist);
    pp = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(p_ptr,pde,assembler,*parlist);
    zp = HDSA::makePtr<PDE_PrimalOptVector<RealT> >(z_ptr,pde,assembler,*parlist);
    rp = HDSA::makePtr<PDE_DualSimVector<RealT> >(r_ptr,pde,assembler,*parlist);

    // Initialize quadratic objective function.
    std::vector<HDSA::Ptr<QoI<RealT> > > qoi_vec(2,ROL::nullPtr);
    qoi_vec[0] = HDSA::makePtr<QoI_Vertical_Velocity_Stokes<RealT> >(pde_stokes->getVelocityFE(),
								     pde_stokes->getPressureFE(),
								     pde_stokes->getFieldHelper(),
								     *parlist);
    qoi_vec[1] = HDSA::makePtr<QoI_L2Penalty_Stokes<RealT> >(pde_stokes->getVelocityFE(),
							     pde_stokes->getPressureFE(),
							     pde_stokes->getFieldHelper());
    HDSA::Ptr<StdObjective_Stokes<RealT> > std_obj = HDSA::makePtr<StdObjective_Stokes<RealT> >(*parlist);
    HDSA::Ptr<ROL::Objective_SimOpt<RealT> > obj = HDSA::makePtr<PDE_Objective<RealT> >(qoi_vec,std_obj,assembler);
    HDSA::Ptr<ROL::Reduced_Objective_SimOpt<RealT> > robj = HDSA::makePtr<ROL::Reduced_Objective_SimOpt<RealT> >(obj, con, up, zp, pp, true, false);
    
    // Initialize PDE describing Navier-Stokes equations.
    HDSA::Ptr<PDE_Stokes<RealT> > pde_ns = HDSA::makePtr<PDE_Stokes<RealT> >(*parlist,1.0);
    con_ns = HDSA::makePtr<PDE_Constraint<RealT> >(pde_ns,meshMgr,comm->Get_Teuchos_Communicator(),*parlist,*outStream);
    con_ns->setSolveParameters(*parlist);

    // Run derivative checks
    bool checkDeriv = parlist->sublist("Problem").get("Check derivatives",false);
    if ( checkDeriv ) {
      ROL::Vector_SimOpt<RealT> x(up,zp);
      HDSA::Ptr<ROL::Vector<RealT> > dup = up->clone(); dup->randomize(-1.0,1.0);
      HDSA::Ptr<ROL::Vector<RealT> > dzp = zp->clone(); dzp->randomize(-1.0,1.0);
      ROL::Vector_SimOpt<RealT> d(dup,dzp);
      obj->checkGradient(x,d,true,*outStream);
      obj->checkHessVec(x,d,true,*outStream);
      con->checkApplyJacobian_1(*up,*zp,*dup,*up,true,*outStream);
      con->checkApplyJacobian_2(*up,*zp,*dzp,*up,true,*outStream);
      con->checkApplyJacobian(x,d,*up,true,*outStream);
      con->checkApplyAdjointHessian(x,*dup,d,x,true,*outStream);
      con->checkAdjointConsistencyJacobian(*dup,d,x,true,*outStream);
      con->checkInverseJacobian_1(*up,*up,*up,*zp,true,*outStream);
      con->checkInverseAdjointJacobian_1(*up,*up,*up,*zp,true,*outStream);
      robj->checkGradient(*zp,*dzp,true,*outStream);
      robj->checkHessVec(*zp,*dzp,true,*outStream);
    }
    
    HDSA::Opt_Problem_Objects<RealT>::fs_obj = HDSA::makePtr<ROL_FS_Objective_Model_Error<RealT> >(obj);
    HDSA::Opt_Problem_Objects<RealT>::con = HDSA::makePtr<ROL_Constraint<RealT> >(con);
    HDSA::Opt_Problem_Objects<RealT>::u = HDSA::makePtr<ROL_PDEOPT_Tpetra_Vector<RealT> >(pde,assembler,true);
    HDSA::Opt_Problem_Objects<RealT>::rs_obj = HDSA::makePtr<ROL_RS_Objective_Model_Error<RealT> >(robj);
    HDSA::Opt_Problem_Objects<RealT>::z = HDSA::makePtr<ROL_PDEOPT_Tpetra_Vector<RealT> >(pde,assembler);
    HDSA::Opt_Problem_Objects<RealT>::theta = theta_in;
  }

  virtual ~Opt_Problem_Objects_stokes()
  { }

  HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > Construct_Opt_Problem_Objects(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<const HDSA::Comm<int> > & comm) 
  { 
    HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects =  HDSA::makePtr<Opt_Problem_Objects_stokes<RealT> >(parlist,theta,comm);
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
	std::cout << "Error loading the optimal control solution" << std::endl;
      }    

    // read in solution and write to Opt_Problem_Objects<RealT>::u
    std::ifstream inputFile2("state_read.txt");          
    count = 0;
    // read the elements in the file into a vector  
    // test file open   
    if (inputFile2) {   
      while ( inputFile2 >> value ) {
	HDSA::Opt_Problem_Objects<RealT>::u->Replace_Element(count,value);
	count += 1;
      }
    }
    else
      {
	std::cout << "Error loading the optimal state solution" << std::endl;
      }   
  }
  
  void Write_Optimal_Solution()
  {
    assembler->printMeshData(*outStream);
    pdecon->outputTpetraData();
    RealT tol(1.e-8);
    Teuchos::Array<RealT> res(1,0);

    up->zero(); zp->zero(); rp->zero();
    con->solve(*rp,*up,*zp,tol);
    pdecon->outputTpetraVector(u_ptr,"uncontrolled_state.txt");
    con->value(*rp,*up,*zp,tol);
    r_ptr->norm2(res.view(0,1));
    *outStream << "Residual Norm: " << res[0] << std::endl;

    // I'm getting a segmentation fault due in parallel when accessing the get_rol_vec() vector
    zp->set(*dynamic_cast<const ROL_Vector<RealT>&>(*HDSA::Opt_Problem_Objects<RealT>::z).get_rol_vec());
    con->solve(*rp,*up,*zp,tol);
    pdecon->outputTpetraVector(u_ptr,"state.txt");
    pdecon->outputTpetraVector(z_ptr,"control.txt");
    con->value(*rp,*up,*zp,tol);
    r_ptr->norm2(res.view(0,1));
    *outStream << "Residual Norm: " << res[0] << std::endl;

    con_ns->solve(*rp,*up,*zp,tol);
    pdecon->outputTpetraVector(u_ptr,"hifi_state.txt");
    con_ns->value(*rp,*up,*zp,tol);
    r_ptr->norm2(res.view(0,1));
    *outStream << "Residual Norm: " << res[0] << std::endl;
  }

};


#endif
