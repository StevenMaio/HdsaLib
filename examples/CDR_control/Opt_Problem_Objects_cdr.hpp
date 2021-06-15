#ifndef OPT_PROBLEM_OBJECTS_CDR_HPP
#define OPT_PROBLEM_OBJECTS_CDR_HPP

// ROL instantiation of Opt_Problem_Objects

template <class RealT>
class Opt_Problem_Objects_CDR : public Opt_Problem_Objects_ROL<RealT> {

private:

  const HDSA::Ptr<HDSA::ParameterList > parlist;
  const HDSA::Ptr<const HDSA::Comm<int> > comm_;
  HDSA::Ptr<std::ostream> outStream;
  HDSA::Ptr<ROL::Objective_SimOpt<RealT> > obj;
  HDSA::Ptr<ROL::Reduced_Objective_SimOpt<RealT> > robj;
  HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > con;
  HDSA::Ptr<PDE_Constraint<RealT> > pdecon; 
  HDSA::Ptr<Assembler<RealT> > assembler;
  HDSA::Ptr<Tpetra::MultiVector<> > u_ptr, z_ptr, p_ptr, r_ptr;
  HDSA::Ptr<ROL::Vector<RealT> > up, zp, pp, rp;

public:
  
  Opt_Problem_Objects_CDR(const HDSA::Ptr<HDSA::ParameterList > & parlist_in, const HDSA::Ptr<const HDSA::Comm<int> > & comm): 
    Opt_Problem_Objects_ROL<RealT>(parlist_in,comm), parlist(parlist_in), comm_(comm)
  { }

  Opt_Problem_Objects_CDR(const HDSA::Ptr<HDSA::ParameterList > & parlist_in, const HDSA::Ptr<HDSA::Vector<RealT> > & theta_in, const HDSA::Ptr<const HDSA::Comm<int> > & comm):
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

    // Initialize PDE
    HDSA::Ptr<MeshManager<RealT> > meshMgr = HDSA::makePtr<MeshManager_CDR<RealT> >(*parlist);
    HDSA::Ptr<PDE<RealT> > pde = HDSA::makePtr<PDE_CDR<RealT> >(*parlist);
    HDSA::Ptr<PDE_CDR<RealT> > pde_cdr = HDSA::dynamicPtrCast<PDE_CDR<RealT> >(pde);
    con = HDSA::makePtr<PDE_Constraint<RealT> >(pde,meshMgr,comm->Get_Teuchos_Communicator(),*parlist,*outStream);
    // Cast the constraint and get the assembler.
    pdecon = HDSA::dynamicPtrCast<PDE_Constraint<RealT> >(con);
    assembler = pdecon->getAssembler();

    // Create state vector and set to zeroes
    u_ptr  = assembler->createStateVector();   u_ptr->putScalar(0.0);
    z_ptr  = assembler->createControlVector(); z_ptr->putScalar(0.0);
    p_ptr  = assembler->createStateVector();   p_ptr->putScalar(0.0);
    r_ptr  = assembler->createStateVector();   r_ptr->putScalar(0.0);
    up  = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(u_ptr,pde,assembler);
    zp  = HDSA::makePtr<PDE_PrimalOptVector<RealT> >(z_ptr,pde,assembler);
    pp  = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(p_ptr,pde,assembler);
    rp  = HDSA::makePtr<PDE_DualSimVector<RealT> >(r_ptr,pde,assembler);
    
    /*************************************************************************/
    /***************** BUILD COST FUNCTIONAL *********************************/
    /*************************************************************************/

    std::vector<ROL::Ptr<QoI<RealT> > > qoi_vec(3,ROL::nullPtr);
    qoi_vec[0] = ROL::makePtr<QoI_L2Tracking_CDR<RealT> >(pde_cdr->getFE());
    qoi_vec[1] = ROL::makePtr<QoI_L2Penalty_CDR<RealT> >(pde_cdr->getFE());
    qoi_vec[2] = ROL::makePtr<QoI_H1Penalty_CDR<RealT> >(pde_cdr->getFE());
    RealT alpha1 = parlist->sublist("Problem").get("L2 penalty parameter",1e-2);
    RealT alpha2 = parlist->sublist("Problem").get("H1 penalty parameter",1e-2);
    std::vector<RealT> wt(3); wt[0] = static_cast<RealT>(1); wt[1] = alpha1; wt[2] = alpha2;
    obj = HDSA::makePtr<PDE_Objective<RealT> >(qoi_vec,wt,assembler);
    // Initialize reduced objective function
    robj = HDSA::makePtr<ROL::Reduced_Objective_SimOpt<RealT> >(obj, con, up, zp, pp);

    // Run derivative checks
    bool checkDeriv = parlist->sublist("Problem").get("Check Derivatives",false);
    if ( checkDeriv ) {
      
      // Create state vector and set to zeroes
      HDSA::Ptr<Tpetra::MultiVector<> > du_ptr = assembler->createStateVector();     du_ptr->randomize();
      HDSA::Ptr<ROL::Vector<RealT> > dup = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(du_ptr,pde,assembler);
      
      // Create control vectors
      HDSA::Ptr<Tpetra::MultiVector<> > dz_ptr = assembler->createControlVector();     dz_ptr->randomize();
      HDSA::Ptr<ROL::Vector<RealT> > dzp = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(dz_ptr,pde,assembler);
      
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

    // Initialize reduced objective function
    HDSA::Opt_Problem_Objects<RealT>::fs_obj = HDSA::makePtr<ROL_FS_Objective_Model_Error<RealT> >(obj);
    HDSA::Opt_Problem_Objects<RealT>::con = HDSA::makePtr<ROL_Constraint<RealT> >(con);
    HDSA::Opt_Problem_Objects<RealT>::u = HDSA::makePtr<ROL_PDEOPT_Tpetra_Vector<RealT> >(pde,assembler);
    HDSA::Opt_Problem_Objects<RealT>::rs_obj = HDSA::makePtr<ROL_RS_Objective_Model_Error<RealT> >(robj);
    HDSA::Opt_Problem_Objects<RealT>::z = HDSA::makePtr<ROL_PDEOPT_Tpetra_Vector<RealT> >(pde,assembler);
    HDSA::Opt_Problem_Objects<RealT>::theta = theta_in;
  }

  virtual ~Opt_Problem_Objects_CDR()
  { }

  HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > Construct_Opt_Problem_Objects(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<const HDSA::Comm<int> > & comm) 
  { 
    HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects =  HDSA::makePtr<Opt_Problem_Objects_CDR<RealT> >(parlist,theta,comm);
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
