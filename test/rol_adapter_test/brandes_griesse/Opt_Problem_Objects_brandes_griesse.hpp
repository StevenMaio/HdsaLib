#ifndef OPT_PROBLEM_OBJECTS_BRANDES_GRIESSE_HPP
#define OPT_PROBLEM_OBJECTS_BRANDES_GRIESSE_HPP

#include "pde_brandes_griesse_parameter.hpp"
#include "ROL_B_Transpose_brandes_griesse.hpp"

// ROL instantiation of Opt_Problem_Objects

template <class RealT>
class Opt_Problem_Objects_brandes_griesse : public Opt_Problem_Objects_ROL<RealT> {

private:

  const HDSA::Ptr<HDSA::ParameterList > parlist;
  HDSA::Ptr<std::ostream> outStream;
  HDSA::Ptr<ROL::Objective_SimOpt<RealT> > obj;
  HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > con;
  HDSA::Ptr<Linear_PDE_Constraint<RealT> > pdecon;
  HDSA::Ptr<Assembler<RealT> > assembler;
  HDSA::Ptr<Tpetra::MultiVector<> > u_ptr, z_ptr, p_ptr, r_ptr;
  HDSA::Ptr<ROL::Vector<RealT> > up, zp, pp, rp;

public:

  Opt_Problem_Objects_brandes_griesse(const HDSA::Ptr<HDSA::ParameterList > & parlist_in, const HDSA::Ptr<const HDSA::Comm<int> > & comm): 
    Opt_Problem_Objects_ROL<RealT>(parlist_in,comm), parlist(parlist_in)
  { }

  Opt_Problem_Objects_brandes_griesse(const HDSA::Ptr<HDSA::ParameterList > & parlist_in, const HDSA::Ptr<HDSA::Vector<RealT> > & theta_in, const HDSA::Ptr<const HDSA::Comm<int> > & comm):
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
    HDSA::Ptr<PDE<RealT> > pde = HDSA::makePtr<PDE_Brandes_Griesse<RealT> >(*parlist);
    HDSA::Ptr<PDE_Brandes_Griesse<RealT> > pde_brandes_griesse = HDSA::dynamicPtrCast<PDE_Brandes_Griesse<RealT> >(pde);
    con = HDSA::makePtr<Linear_PDE_Constraint<RealT> >(pde,meshMgr,comm->Get_Teuchos_Communicator(),*parlist,*outStream,true);
    // Cast the constraint and get the assembler.
    pdecon = HDSA::dynamicPtrCast<Linear_PDE_Constraint<RealT> >(con);
    assembler = pdecon->getAssembler();
    // Initialize quadratic objective function
    std::vector<HDSA::Ptr<QoI<RealT> > > qoi_vec(2,HDSA::nullPtr);
    qoi_vec[0] = HDSA::makePtr<QoI_L2Tracking_Brandes_Griesse<RealT> >(pde_brandes_griesse->getFE());
    qoi_vec[1] = HDSA::makePtr<QoI_L2Penalty_Brandes_Griesse<RealT> >(pde_brandes_griesse->getFE());
    RealT alpha = parlist->sublist("Problem").get("Control penalty parameter",1e-2);
    std::vector<RealT> wt(2); wt[0] = static_cast<RealT>(1); wt[1] = alpha;
    obj = HDSA::makePtr<PDE_Objective<RealT> >(qoi_vec,wt,assembler);
    pdecon->setParameter(*(dynamic_cast<const Std_Vector<RealT>&>(*theta_in).get_std_vec()));
    // Create state vector and set to zeroes
    u_ptr  = assembler->createStateVector();   u_ptr->putScalar(0.0);
    z_ptr  = assembler->createControlVector(); z_ptr->putScalar(0.0);
    p_ptr  = assembler->createStateVector();   p_ptr->putScalar(0.0);
    r_ptr  = assembler->createStateVector();   r_ptr->putScalar(0.0);
    up  = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(u_ptr,pde,assembler);
    zp  = HDSA::makePtr<PDE_PrimalOptVector<RealT> >(z_ptr,pde,assembler);
    pp  = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(p_ptr,pde,assembler);
    rp  = HDSA::makePtr<PDE_DualSimVector<RealT> >(r_ptr,pde,assembler);

    bool  use_Full_Space = parlist->sublist("Problem").get("Full Space",false);

    if(use_Full_Space)
      {
	HDSA::Opt_Problem_Objects<RealT>::fs_obj = HDSA::makePtr<ROL_FS_Objective<RealT> >(obj);
	HDSA::Opt_Problem_Objects<RealT>::con = HDSA::makePtr<ROL_Constraint<RealT> >(con);
	HDSA::Opt_Problem_Objects<RealT>::u = HDSA::makePtr<ROL_PDEOPT_Tpetra_Vector<RealT> >(pde,assembler);
	HDSA::Opt_Problem_Objects<RealT>::z = HDSA::makePtr<ROL_PDEOPT_Tpetra_Vector<RealT> >(pde,assembler);
	HDSA::Opt_Problem_Objects<RealT>::lambda = HDSA::makePtr<ROL_PDEOPT_Tpetra_Vector<RealT> >(pde,assembler);
	HDSA::Opt_Problem_Objects<RealT>::theta = theta_in;
      }
    else
      {
	// Initialize reduced objective function
	HDSA::Ptr<ROL::Objective<RealT> > robj = HDSA::makePtr<ROL::Reduced_Objective_SimOpt<RealT> >(obj, con, up, zp, pp);
  
	HDSA::Ptr<ROL_B_Transpose<RealT> > rol_Bt = HDSA::makePtr<ROL_B_Transpose_brandes_griesse<RealT> >(meshMgr, assembler,parlist, outStream, comm->Get_Teuchos_Communicator(), obj, con, up, zp, theta_in);
	
	HDSA::Opt_Problem_Objects<RealT>::rs_obj = HDSA::makePtr<ROL_RS_Objective<RealT> >(robj, rol_Bt);
	HDSA::Opt_Problem_Objects<RealT>::z = HDSA::makePtr<ROL_PDEOPT_Tpetra_Vector<RealT> >(pde,assembler);
	HDSA::Opt_Problem_Objects<RealT>::theta = theta_in;
      }
  }

  virtual ~Opt_Problem_Objects_brandes_griesse()
  { }

  HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > Construct_Opt_Problem_Objects(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<const HDSA::Comm<int> > & comm) 
  { 
    HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects =  HDSA::makePtr<Opt_Problem_Objects_brandes_griesse<RealT> >(parlist,theta,comm);
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

    bool  use_Full_Space = parlist->sublist("Problem").get("Full Space",false);
    if(use_Full_Space)
      {
        std::ifstream inputFile_SR("state_read.txt"); 
	value = 0.0;
	count = 0;
	// read the elements in the file into a vector  
	// test file open   
	if (inputFile_SR) {   
	  while ( inputFile_SR >> value ) {
	    HDSA::Opt_Problem_Objects<RealT>::u->Replace_Element(count,value);
	    count += 1;
	  }
	}
	else
	  {
	    std::cout << "Error loading the optimal state solution" << std::endl;
	  }  

	std::ifstream inputFile_AR("adjoint_read.txt"); 
	value = 0.0;
	count = 0;
	// read the elements in the file into a vector  
	// test file open   
	if (inputFile_AR) {   
	  while ( inputFile_AR >> value ) {
	    HDSA::Opt_Problem_Objects<RealT>::lambda->Replace_Element(count,value);
	    count += 1;
	  }
	}
	else
	  {
	    std::cout << "Error loading the optimal adjoint solution" << std::endl;
	  }  

      }
  
  }
  
  void Write_Optimal_Solution()
  {
    bool  use_Full_Space = parlist->sublist("Problem").get("Full Space",false);

    if(use_Full_Space)
      {
	up->set(*dynamic_cast<const ROL_Vector<RealT>&>(*HDSA::Opt_Problem_Objects<RealT>::u).get_rol_vec());
	zp->set(*dynamic_cast<const ROL_Vector<RealT>&>(*HDSA::Opt_Problem_Objects<RealT>::z).get_rol_vec());
	// Output.
	pdecon->outputTpetraVector(u_ptr,"state.txt");
	pdecon->outputTpetraVector(z_ptr,"control.txt");
	pdecon->outputTpetraData();
	assembler->printMeshData(*outStream);

	// Compute Adjoint
	RealT tol = 1.e-8;
	HDSA::Ptr<ROL::Vector<RealT> > dualstate = up->dual().clone();
	dualstate->set(*up);
        obj->gradient_1(*dualstate,*up,*zp,tol);
	con->applyInverseAdjointJacobian_1(*pp,*dualstate,*up,*zp,tol);
	pp->scale(-1.0);
	pdecon->outputTpetraVector(p_ptr,"adjoint.txt");
      }
    else
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
  }

};


#endif
