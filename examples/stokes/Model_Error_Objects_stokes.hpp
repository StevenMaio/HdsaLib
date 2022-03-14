#ifndef BAYES_MODEL_ERROR_OBJECTS_STOKES_HPP
#define BAYES_MODEL_ERROR_OBJECTS_STOKES_HPP

// Instantiation of Bayes_Model_Error_Objects

template <class RealT>
class Bayes_Model_Error_Objects_stokes : public HDSA::Bayes_Model_Error_Objects<RealT> {

private:
  HDSA::Ptr<HDSA::ParameterList> parlist_;
  HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > con_AG_;
  HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > con_AL_;
  HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > con_M_;
  HDSA::Ptr<HDSA::Vector<RealT> > M_row_sum_;

public:

  Bayes_Model_Error_Objects_stokes(const HDSA::Ptr<HDSA::ParameterList> & parlist_sensitivity, const HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > & OP_Objects_Factory, 
				   const HDSA::Ptr<HDSA::Weight_Matrices<RealT> > & weight_matrices_factory, const HDSA::Ptr<HDSA::ParameterList> & parlist):
    HDSA::Bayes_Model_Error_Objects<RealT>(parlist_sensitivity,OP_Objects_Factory,weight_matrices_factory)
  {
    parlist_ = parlist;
  }

  virtual ~Bayes_Model_Error_Objects_stokes()
  { }

  void Construct_Objects(const HDSA::Ptr<const HDSA::Comm<int> > & comm)
  {
    HDSA::Ptr<std::ostream> outStream;
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
    HDSA::Ptr<MeshManager<RealT> > meshMgr = HDSA::makePtr<MeshManager_BackwardFacingStepChannel<RealT> >(*parlist_);

    HDSA::Ptr<PDE<RealT> > pde_AG = HDSA::makePtr<PDE_Elliptic<RealT> >(*parlist_,1.e-5,1.0);
    con_AG_ = HDSA::makePtr<PDE_Constraint<RealT> >(pde_AG,meshMgr,comm->Get_Teuchos_Communicator(),*parlist_,*outStream);
    con_AG_->setSolveParameters(*parlist_);

    HDSA::Ptr<PDE<RealT> > pde_AL = HDSA::makePtr<PDE_Elliptic<RealT> >(*parlist_,1.0,1.0);
    con_AL_ = HDSA::makePtr<PDE_Constraint<RealT> >(pde_AL,meshMgr,comm->Get_Teuchos_Communicator(),*parlist_,*outStream);
    con_AL_->setSolveParameters(*parlist_);

    HDSA::Ptr<PDE<RealT> > pde_M = HDSA::makePtr<PDE_Elliptic<RealT> >(*parlist_,0.0,1.0);
    con_M_ = HDSA::makePtr<PDE_Constraint<RealT> >(pde_M,meshMgr,comm->Get_Teuchos_Communicator(),*parlist_,*outStream);
    con_M_->setSolveParameters(*parlist_);

    M_row_sum_ = HDSA::Model_Error_Objects<RealT>::OP_Objects_->u->Clone();
    HDSA::Ptr<HDSA::Vector<RealT> > tmp = M_row_sum_->Clone();
    RealT tol = 1.e-8;
    HDSA::Ptr<ROL::Vector<RealT> > v1 = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(M_row_sum_)->get_rol_vec();
    HDSA::Ptr<ROL::Vector<RealT> > v2 = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(tmp)->get_rol_vec();
    v2->setScalar(1.0);
    con_M_->applyJacobian_1(*v1,*v2,*v2,*v2,tol);

    // Define Gamma = A_G^{-1}*M*A_G^{-1}
    // Define L = A_L*M^{-1}*A_L
  }
  
  void Apply_A_G_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const
  {
    HDSA::Ptr<ROL::Vector<RealT> > z_in_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(z_in)->get_rol_vec();
    HDSA::Ptr<ROL::Vector<RealT> > z_out_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(z_out)->get_rol_vec();
    RealT tol = 1.e-8;
    con_AG_->applyJacobian_1(*z_out_rol,*z_in_rol,*z_in_rol,*z_in_rol,tol);
  }

  void Apply_A_G_Mat_Inv(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const
  {
    HDSA::Ptr<ROL::Vector<RealT> > z_in_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(z_in)->get_rol_vec();
    HDSA::Ptr<ROL::Vector<RealT> > z_out_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(z_out)->get_rol_vec();
    RealT tol = 1.e-8;
    con_AG_->applyInverseJacobian_1(*z_out_rol,*z_in_rol,*z_in_rol,*z_in_rol,tol);
  }

  void Apply_A_L_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in) const
  {
    HDSA::Ptr<ROL::Vector<RealT> > u_in_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(u_in)->get_rol_vec();
    HDSA::Ptr<ROL::Vector<RealT> > u_out_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(u_out)->get_rol_vec();
    RealT tol = 1.e-8;
    con_AL_->applyJacobian_1(*u_out_rol,*u_in_rol,*u_in_rol,*u_in_rol,tol);
  }

  void Apply_A_L_Mat_Inv(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in) const
  {
    HDSA::Ptr<ROL::Vector<RealT> > u_in_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(u_in)->get_rol_vec();
    HDSA::Ptr<ROL::Vector<RealT> > u_out_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(u_out)->get_rol_vec();
    RealT tol = 1.e-8;
    con_AL_->applyInverseJacobian_1(*u_out_rol,*u_in_rol,*u_in_rol,*u_in_rol,tol);
  }

  void Apply_Gamma_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const 
  {
    HDSA::Ptr<HDSA::Vector<RealT> > tmp = z_in->Clone();
    Apply_A_G_Mat_Inv(tmp,z_in);
    for(int k = 0; k < tmp->dimension(); k++)
      {
	RealT val = (*tmp)(k)*(*M_row_sum_)(k);
	tmp->Replace_Element(k,val);
      }
    Apply_A_G_Mat_Inv(z_out,tmp); 
  }

  void Apply_Gamma_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const 
  {
    HDSA::Ptr<HDSA::Vector<RealT> > tmp = z_in->Clone();
    Apply_A_G_Mat(tmp,z_in);
    for(int k = 0; k < tmp->dimension(); k++)
      {
	RealT val = (*tmp)(k)/(*M_row_sum_)(k);
	tmp->Replace_Element(k,val);
      }
    Apply_A_G_Mat(z_out,tmp); 
  }

  void Apply_L_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in) const 
  {
    HDSA::Ptr<HDSA::Vector<RealT> > tmp = u_in->Clone();
    Apply_A_L_Mat(tmp,u_in);
    for(int k = 0; k < tmp->dimension(); k++)
      {
	RealT val = (*tmp)(k)/(*M_row_sum_)(k);
	tmp->Replace_Element(k,val);
      }
    Apply_A_L_Mat(u_out,tmp); 
  }

  void Apply_L_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in) const 
  {
    HDSA::Ptr<HDSA::Vector<RealT> > tmp = u_in->Clone();
    Apply_A_L_Mat_Inv(tmp,u_in);
    for(int k = 0; k < tmp->dimension(); k++)
      {
	RealT val = (*tmp)(k)*(*M_row_sum_)(k);
	tmp->Replace_Element(k,val);
      }
    Apply_A_L_Mat_Inv(u_out,tmp); 
  }

  void Apply_Sqrt_Gamma_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const 
  {
    for(int k = 0; k < z_in->dimension(); k++)
      {
	RealT val = (*z_in)(k)*std::sqrt((*M_row_sum_)(k));
	z_in->Replace_Element(k,val);
      }
    Apply_A_G_Mat_Inv(z_out,z_in);
  }

  void Apply_Sqrt_Gamma_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const 
  {
    for(int k = 0; k < z_in->dimension(); k++)
      {
	RealT val = (*z_in)(k)/std::sqrt((*M_row_sum_)(k));
	z_in->Replace_Element(k,val);
      }
    Apply_A_G_Mat(z_out,z_in);
  }
  
  void Apply_Sqrt_L_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in) const 
  {
    for(int k = 0; k < u_in->dimension(); k++)
      {
	RealT val = (*u_in)(k)*std::sqrt((*M_row_sum_)(k));
	u_in->Replace_Element(k,val);
      }
    Apply_A_L_Mat_Inv(u_out,u_in);
  }
  
  // NEED TO IMPLEMENT
  void Apply_L_Plus_Shift_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in, const RealT beta) const 
  {
    u_out->set(*u_in);
  }

  // NEED TO IMPLEMENT
  void Apply_Sqrt_L_Plus_Shift_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in, const RealT beta) const 
  {
    u_out->set(*u_in);
  }

  HDSA::Ptr<HDSA::MultiVector<RealT> > Load_Z_Data(void) const 
  {
    int num_vecs = 1;
    int dim = HDSA::Model_Error_Objects<RealT>::OP_Objects_->z->dimension();
    HDSA::Ptr<HDSA::MultiVector<RealT> > Z = HDSA::makePtr<HDSA::MultiVector<RealT> >(num_vecs,HDSA::Model_Error_Objects<RealT>::OP_Objects_->z);
    // read in data
    std::ifstream in_Z("Z.txt");           
    // read the elements in the file into a vector  
    // test file open   
    RealT val;
    if (in_Z) {   
      for(int i = 0; i < dim; i++)
	{
	  for(int j = 0; j < num_vecs; j++)
	    {
	      in_Z >> val;
	      (*Z)[j]->Replace_Element(i,val);
	    }
	}   
    }
    else
      {
	std::cout << "Error loading the data from Z.txt" << std::endl;
      } 
    return Z;
  }

  HDSA::Ptr<HDSA::MultiVector<RealT> > Load_Y_Data(void) const 
  {
    int num_vecs = 1;
    int dim = HDSA::Model_Error_Objects<RealT>::OP_Objects_->u->dimension();
    HDSA::Ptr<HDSA::MultiVector<RealT> > Y = HDSA::makePtr<HDSA::MultiVector<RealT> >(num_vecs,HDSA::Model_Error_Objects<RealT>::OP_Objects_->u);
    // read in data
    std::ifstream in_Y("Y.txt");           
    // read the elements in the file into a vector  
    // test file open   
    RealT val;
    if (in_Y) {   
      for(int i = 0; i < dim; i++)
	{
	  for(int j = 0; j < num_vecs; j++)
	    {
	      in_Y >> val;
	      (*Y)[j]->Replace_Element(i,val);
	    }
	}   
    }
    else
      {
	std::cout << "Error loading the data from Y.txt" << std::endl;
      } 
    return Y;
  }

};


#endif
