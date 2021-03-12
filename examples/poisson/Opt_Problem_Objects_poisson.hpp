#ifndef OPT_PROBLEM_OBJECTS_POISSON_HPP
#define OPT_PROBLEM_OBJECTS_POISSON_HPP

// ROL instantiation of Opt_Problem_Objects

template <class RealT>
class Opt_Problem_Objects_poisson : public Opt_Problem_Objects_ROL<RealT> {

private:

  const HDSA::Ptr<HDSA::ParameterList > parlist;
  const HDSA::Ptr<const HDSA::Comm<int> > comm_;
  HDSA::Ptr<std::ostream> outStream;
  HDSA::Ptr<ROL::Objective_SimOpt<RealT> > obj;
  HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > con;
  HDSA::Ptr<Linear_PDE_Constraint<RealT> > pdecon;
  HDSA::Ptr<Assembler<RealT> > assembler;
  HDSA::Ptr<Tpetra::MultiVector<> > u_ptr, z_ptr, p_ptr, r_ptr;
  HDSA::Ptr<ROL::Vector<RealT> > up, zp, pp, rp;

  HDSA::Ptr<Tpetra::MultiVector<> > Wd_ptr;
  HDSA::Ptr<ROL::Vector<RealT> > Wdp;

public:

  Opt_Problem_Objects_poisson(const HDSA::Ptr<HDSA::ParameterList > & parlist_in, const HDSA::Ptr<const HDSA::Comm<int> > & comm): 
    Opt_Problem_Objects_ROL<RealT>(parlist_in,comm), parlist(parlist_in), comm_(comm)
  { }

  Opt_Problem_Objects_poisson(const HDSA::Ptr<HDSA::ParameterList > & parlist_in, const HDSA::Ptr<HDSA::Vector<RealT> > & theta_in, const HDSA::Ptr<const HDSA::Comm<int> > & comm):
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
    HDSA::Ptr<MeshManager<RealT> > meshMgr = HDSA::makePtr<MeshManager_Rectangle<RealT> >(*parlist);
    HDSA::Ptr<PDE<RealT> > pde = HDSA::makePtr<PDE_poisson<RealT> >(*parlist);
    HDSA::Ptr<PDE_poisson<RealT> > pde_poisson = HDSA::dynamicPtrCast<PDE_poisson<RealT> >(pde);
    con = HDSA::makePtr<Linear_PDE_Constraint<RealT> >(pde,meshMgr,comm->Get_Teuchos_Communicator(),*parlist,*outStream,true);
    // Cast the constraint and get the assembler.
    pdecon = HDSA::dynamicPtrCast<Linear_PDE_Constraint<RealT> >(con);
    assembler = pdecon->getAssembler();
    // Initialize quadratic objective function
    std::vector<HDSA::Ptr<QoI<RealT> > > qoi_vec(2,HDSA::nullPtr);
    qoi_vec[0] = HDSA::makePtr<QoI_L2Tracking_poisson<RealT> >(pde_poisson->getFE());
    qoi_vec[1] = HDSA::makePtr<QoI_L2Penalty_poisson<RealT> >(pde_poisson->getFE());
    RealT alpha = parlist->sublist("Problem").get("Control penalty parameter",1e-2);
    std::vector<RealT> wt(2); wt[0] = static_cast<RealT>(1); wt[1] = alpha;
    obj = HDSA::makePtr<PDE_Objective<RealT> >(qoi_vec,wt,assembler);
    //pdecon->setParameter(*(dynamic_cast<const Std_Vector<RealT>&>(*theta_in).get_std_vec()));
    // Create state vector and set to zeroes
    u_ptr  = assembler->createStateVector();   u_ptr->putScalar(0.0);
    z_ptr  = assembler->createControlVector(); z_ptr->putScalar(0.0);
    p_ptr  = assembler->createStateVector();   p_ptr->putScalar(0.0);
    r_ptr  = assembler->createStateVector();   r_ptr->putScalar(0.0);
    up  = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(u_ptr,pde,assembler);
    zp  = HDSA::makePtr<PDE_PrimalOptVector<RealT> >(z_ptr,pde,assembler);
    pp  = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(p_ptr,pde,assembler);
    rp  = HDSA::makePtr<PDE_DualSimVector<RealT> >(r_ptr,pde,assembler);

    Wd_ptr  = assembler->createStateVector();   Wd_ptr->putScalar(0.0);
    Wdp  = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(Wd_ptr,pde,assembler);
    RealT tol = 1.e-8;
    obj->gradient_1(*Wdp, *up, *zp, tol);

    // Initialize reduced objective function
    HDSA::Opt_Problem_Objects<RealT>::fs_obj = HDSA::makePtr<ROL_FS_Objective_Model_Error<RealT> >(obj);
    HDSA::Opt_Problem_Objects<RealT>::con = HDSA::makePtr<ROL_Constraint<RealT> >(con);
    HDSA::Opt_Problem_Objects<RealT>::u = HDSA::makePtr<ROL_PDEOPT_Tpetra_Vector<RealT> >(pde,assembler);
    HDSA::Ptr<ROL::Objective<RealT> > robj = HDSA::makePtr<ROL::Reduced_Objective_SimOpt<RealT> >(obj, con, up, zp, pp);
    HDSA::Opt_Problem_Objects<RealT>::rs_obj = HDSA::makePtr<ROL_RS_Objective_Model_Error<RealT> >(robj);
    HDSA::Opt_Problem_Objects<RealT>::z = HDSA::makePtr<ROL_PDEOPT_Tpetra_Vector<RealT> >(pde,assembler);
    HDSA::Opt_Problem_Objects<RealT>::theta = theta_in;
  }

  virtual ~Opt_Problem_Objects_poisson()
  { }

  HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > Construct_Opt_Problem_Objects(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<const HDSA::Comm<int> > & comm) 
  { 
    HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects =  HDSA::makePtr<Opt_Problem_Objects_poisson<RealT> >(parlist,theta,comm);
    return OP_Objects;
  }

  void Load_Optimal_Solution() 
  { 
    // read in solution and write to Opt_Problem_Objects<RealT>::z
    std::ifstream inputFile("source_read.txt");          
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
    pdecon->outputTpetraVector(z_ptr,"source.txt");
    pdecon->outputTpetraVector(Wd_ptr,"Wd.txt");
    pdecon->outputTpetraData();
    assembler->printMeshData(*outStream);
  }

  const HDSA::Ptr<const Teuchos::Comm<int> > Get_comm(void) const
  {
    return comm_->Get_Teuchos_Communicator();
  }

};


#endif
