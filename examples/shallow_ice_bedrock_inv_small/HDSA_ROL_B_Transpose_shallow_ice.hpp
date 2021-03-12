#ifndef HDSA_ROL_B_TRANSPOSE_SHALLOW_ICE_HPP
#define HDSA_ROL_B_TRANSPOSE_SHALLOW_ICE_HPP

#include "dynpde_shallow_ice_param.hpp"
#include "pde_shallow_ice_param.hpp"

template <class RealT>
class ROL_B_Transpose_shallow_ice: public ROL_B_Transpose<RealT> {

private:
  HDSA::Ptr<ROL::Objective<RealT> > robj_;

  HDSA::Ptr<std::ostream> outStream;
  HDSA::Ptr<DynamicPDE_shallow_ice_param<RealT> > pde_shallow_ice_;
  HDSA::Ptr<DynConstraint<RealT> > dyn_con;
  HDSA::Ptr<Tpetra::MultiVector<> > u0_ptr, uo_ptr, un_ptr,z_ptr, ck_ptr;
  HDSA::Ptr<ROL::Vector<RealT> > u0, uo, un, ck, z;
  std::vector<ROL::TimeStamp<RealT> > timeStamp;
  std::vector<int> data_weight_id;
  HDSA::Ptr<Assembler<RealT> > assembler_;
  HDSA::Ptr<Intrepid::FieldContainer<RealT> > z_coeff_;
  HDSA::Ptr<std::vector<RealT> > theta_rol_ptr;
  HDSA::Ptr<ROL::Vector<RealT> > theta_rolp;
  HDSA::Ptr<HDSA::Vector<RealT> > theta_grad_nom_;

public:

  ROL_B_Transpose_shallow_ice(const HDSA::Ptr<HDSA::ParameterList > & parlist, const HDSA::Ptr<const HDSA::Comm<int> > & comm, 
			      const HDSA::Ptr<Assembler<RealT> > & assembler, const HDSA::Ptr<HDSA::Vector<RealT> > & theta)
  { 
    theta_grad_nom_ = theta->Clone();
    assembler_ = assembler;
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
    HDSA::Ptr<DynamicPDE<RealT> > pde = HDSA::makePtr<DynamicPDE_shallow_ice_param<RealT> >(*parlist);
    pde_shallow_ice_ = HDSA::dynamicPtrCast<DynamicPDE_shallow_ice_param<RealT> >(pde);  

    /*************************************************************************/
    /***************** BUILD CONSTRAINT **************************************/
    /*************************************************************************/
    dyn_con = HDSA::makePtr<DynConstraint<RealT> >(pde,meshMgr,comm->Get_Teuchos_Communicator(),*parlist,*outStream);

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

    // Create theta parameter vectors
    theta_rol_ptr  = HDSA::makePtr<std::vector<RealT> >(theta->dimension(),0.0);
    theta_rolp = HDSA::makePtr<PDE_OptVector<RealT> >(HDSA::makePtr<ROL::StdVector<RealT> >(theta_rol_ptr));

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
    HDSA::Ptr<Objective_SimOpt_TS<RealT> > misfit_obj = HDSA::makePtr<State_Cost_shallow_ice<RealT> >(timeStamp,data,data_weight,data_weight_id,u0);
    HDSA::Ptr<LTI_Objective_TS<RealT> > dyn_obj = HDSA::makePtr<LTI_Objective_TS<RealT> >(*parlist,misfit_obj);

    /*************************************************************************/
    /***************** BUILD REDUCED COST FUNCTIONAL *************************/
    /*************************************************************************/
    ROL::ParameterList &rpl = parlist->sublist("Reduced Dynamic Objective");
    robj_ = HDSA::makePtr<ROL::ReducedDynamicObjective_Stationary_Control<RealT> >(dyn_obj, dyn_con, u0, theta_rolp, ck, timeStamp, rpl);
  
    /*************************************************************************/
    /***************** RUN VECTOR AND DERIVATIVE CHECKS **********************/
    /*************************************************************************/
    // Run derivative checks
    //Load_Optimal_Solution();
    bool checkDeriv = false; // If set to true, uncomment Load_Optimal_Solution above, keep it commented when not running derivative check
    if ( checkDeriv ) 
      {
	std::cout << "Finite Difference Test for B^T" << std::endl;
	std::cout << "------------------------------------------------------------" << std::endl;
	
	// Create state vector
	HDSA::Ptr<Tpetra::MultiVector<> > du_ptr = dyn_con->getAssembler()->createStateVector();
	HDSA::Ptr<ROL::Vector<RealT> > du = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(du_ptr,pde,*dyn_con->getAssembler());

	// Create control vector
	HDSA::Ptr<std::vector<RealT> > dtheta_rol_ptr  = HDSA::makePtr<std::vector<RealT> >(theta->dimension(),0.0);
	HDSA::Ptr<ROL::Vector<RealT> > dtheta_rolp = HDSA::makePtr<PDE_OptVector<RealT> >(HDSA::makePtr<ROL::StdVector<RealT> >(dtheta_rol_ptr));
	
	un->randomize();
	uo->randomize();
	du->randomize();
	theta_rolp->randomize();
	dtheta_rolp->randomize();
	misfit_obj->Update_current_TS(timeStamp[1].t.at(0));

	*outStream << std::endl << "Check State Gradient of Full Space Objective Function" << std::endl;
	misfit_obj->checkGradient_1(*uo,*theta_rolp,*du);
	*outStream << std::endl << "Check Control Gradient of Full Space Objective Function" << std::endl;
	misfit_obj->checkGradient_2(*uo,*theta_rolp,*dtheta_rolp);
  
	ROL::ValidateFunction<RealT> validate(1,13,20,11,true,*outStream);
	ROL::DynamicConstraintCheck<RealT>::check(*dyn_con,validate,*uo,*un,*theta_rolp);       

	*outStream << std::endl << "Check Gradient of Reduced Objective Function" << std::endl;
	robj_->checkGradient(*theta_rolp,*dtheta_rolp,true,*outStream);
	*outStream << std::endl << "Check Hessian of Reduced Objective Function" << std::endl;
	robj_->checkHessVec(*theta_rolp,*dtheta_rolp,true,*outStream);

	std::cout << "------------------------------------------------------------" << std::endl;

      }
    
  }

  virtual ~ROL_B_Transpose_shallow_ice( )
  { }

  // evaluate the theta,z hessian vector product, i.e. -B^T in HDSA
  void hessVec_theta_z(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, 
		       const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr)
  { 
    RealT tol = 1.e-8;
    const Std_Vector<RealT> &etheta = dynamic_cast<const Std_Vector<RealT>&>(theta);
    HDSA::Ptr<std::vector<RealT> > theta_std = etheta.get_std_vec();
    HDSA::Ptr<ROL::Vector<RealT> > theta_rol = HDSA::makePtr<PDE_OptVector<RealT> >(HDSA::makePtr<ROL::StdVector<RealT> >(theta_std));

    if((z_coeff_ == HDSA::nullPtr))
      {
	z_coeff_ = get_z_Field_Container(*dynamic_cast<const ROL_Vector<RealT>&>(z).get_rol_vec());
	pde_shallow_ice_->Update_Z_input(z_coeff_);
	const Std_Vector<RealT> &egradnom = dynamic_cast<const Std_Vector<RealT>&>(*theta_grad_nom_);
	HDSA::Ptr<std::vector<RealT> > gradnom_std = egradnom.get_std_vec();
	HDSA::Ptr<ROL::Vector<RealT> > grad_nom_rol = HDSA::makePtr<PDE_OptVector<RealT> >(HDSA::makePtr<ROL::StdVector<RealT> >(gradnom_std));
	robj_->update(*theta_rol);
	robj_->gradient(*grad_nom_rol,*theta_rol,tol);
      }
    
    HDSA::Ptr<HDSA::Vector<RealT> > theta_grad_pert = theta_grad_nom_->Clone();
    const Std_Vector<RealT> &egradpert = dynamic_cast<const Std_Vector<RealT>&>(*theta_grad_pert);
    HDSA::Ptr<std::vector<RealT> > gradpert_std = egradpert.get_std_vec();
    HDSA::Ptr<ROL::Vector<RealT> > grad_pert_rol = HDSA::makePtr<PDE_OptVector<RealT> >(HDSA::makePtr<ROL::StdVector<RealT> >(gradpert_std));
    
    HDSA::Ptr<HDSA::Vector<RealT> > z_pert = z.Clone();
    z_pert->set(z);
    z_pert->axpy(1.e-4,v);
    HDSA::Ptr<Intrepid::FieldContainer<RealT> > z_coeff_pert = get_z_Field_Container(*dynamic_cast<const ROL_Vector<RealT>&>(*z_pert).get_rol_vec());
    pde_shallow_ice_->Update_Z_input(z_coeff_pert);
    robj_->update(*theta_rol);

    robj_->gradient(*grad_pert_rol,*theta_rol,tol);
    hv.set(*theta_grad_pert);
    hv.axpy(-1.0,*theta_grad_nom_);
    hv.scale(1.e4);
  }

private:
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

  HDSA::Ptr<Intrepid::FieldContainer<RealT> > get_z_Field_Container(const ROL::Vector<RealT> & z)
  {
    HDSA::Ptr<Intrepid::FieldContainer<RealT> > z_coeff = assembler_->get_z_field_container(dynamic_cast<const ROL::TpetraMultiVector<RealT>&>(z).getVector());
    return z_coeff;
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
	z_ptr->replaceGlobalValue(count,0,value);
	count += 1;
      }
    }
    else
      {
	std::cout << "Error loading the optimal solution" << std::endl;
      } 
    z_coeff_ =  get_z_Field_Container(*z);
    pde_shallow_ice_->Update_Z_input(z_coeff_);
  }

};


#endif
