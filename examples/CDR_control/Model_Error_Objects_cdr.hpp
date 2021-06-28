#ifndef MODEL_ERROR_OBJECTS_CDR_HPP
#define MODEL_ERROR_OBJECTS_CDR_HPP

#include "elliptic_op.hpp"

// Instantiation of Model_Error_Objects

template <class RealT>
class Model_Error_Objects_CDR : public HDSA::Model_Error_Objects<RealT> {

private:
  HDSA::Ptr<HDSA::ParameterList> parlist_;
  HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > con_;
  HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > mass_con_;
  bool constructed_elliptic_op_;
  RealT epsilon_;
  int nx_,ny_;
  std::vector<RealT> con_weights_;
  RealT z_cov_scale_;
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Con_Mat_;
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Con_Mat_Q_;
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Con_Mat_R_;

public:

  Model_Error_Objects_CDR(const HDSA::Ptr<HDSA::ParameterList> & parlist_sensitivity, const HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > & OP_Objects_Factory, 
			      const HDSA::Ptr<HDSA::Weight_Matrices<RealT> > & weight_matrices, const HDSA::Ptr<HDSA::ParameterList> & parlist):
    HDSA::Model_Error_Objects<RealT>(parlist_sensitivity,OP_Objects_Factory,weight_matrices), parlist_(parlist)
  {
    constructed_elliptic_op_ = false;
    epsilon_ = parlist_sensitivity->sublist("Model Error").get("Smoothing Factor", 1.e-6);
    nx_ = parlist->sublist("Geometry").get("NX", 10);
    ny_ = parlist->sublist("Geometry").get("NY", 10);
    RealT cw = parlist->sublist("Problem").get("Constraint Weight", 1.0);
    con_weights_ = std::vector<RealT>(nx_+1,cw);
    z_cov_scale_ = parlist->sublist("Problem").get("Control Variance", 1.0);
  }

  virtual ~Model_Error_Objects_CDR()
  { }

  void Construct_Objects(const HDSA::Ptr<const HDSA::Comm<int> > & comm) 
  {
    // Initialize PDE
    HDSA::Ptr<MeshManager<RealT> > meshMgr = HDSA::makePtr<MeshManager_Rectangle<RealT> >(*parlist_);
    HDSA::Ptr<PDE<RealT> > pde = HDSA::makePtr<Elliptic_Op<RealT> >(*parlist_,epsilon_);
    con_ = HDSA::makePtr<Linear_PDE_Constraint<RealT> >(pde,meshMgr,comm->Get_Teuchos_Communicator(),*parlist_);
    HDSA::Ptr<PDE<RealT> > pde_mass = HDSA::makePtr<PDE_Mass_Mat<RealT> >(*parlist_);
    mass_con_ = HDSA::makePtr<Linear_PDE_Constraint<RealT> >(pde_mass,meshMgr,comm->Get_Teuchos_Communicator(),*parlist_);

    // Construct constraint matrix
    Con_Mat_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nx_+1,nx_+1);
    Con_Mat_Q_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nx_+1,nx_+1);
    Con_Mat_R_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nx_+1,nx_+1);
   
    HDSA::Ptr<Linear_PDE_Constraint<RealT> > pdecon = HDSA::dynamicPtrCast<Linear_PDE_Constraint<RealT> >(con_);
    HDSA::Ptr<Assembler<RealT> > assembler = pdecon->getAssembler();
    HDSA::Ptr<Tpetra::MultiVector<> > u1_ptr  = assembler->createStateVector();   u1_ptr->putScalar(0.0);
    HDSA::Ptr<ROL::Vector<RealT> > u1p  = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(u1_ptr,pde,assembler);
    HDSA::Ptr<Tpetra::MultiVector<> > u2_ptr  = assembler->createStateVector();   u2_ptr->putScalar(0.0);
    HDSA::Ptr<ROL::Vector<RealT> > u2p  = HDSA::makePtr<PDE_PrimalSimVector<RealT> >(u2_ptr,pde,assembler);
    RealT tol = 1.e-8;
    for(int j = 0; j < nx_+1; j++)
      {
	u1p->set(*u2p->basis(j));
	con_->applyInverseJacobian_1(*u2p,*u1p,*u1p,*u1p,tol);
	mass_con_->applyJacobian_1(*u1p,*u2p,*u1p,*u1p,tol);
	con_->applyInverseJacobian_1(*u2p,*u1p,*u1p,*u1p,tol);
	for(int i = 0; i < nx_+1; i++)
	  {
	    Con_Mat_->Replace_Element(i,j,u2p->dot(*u1p->basis(i)));
	  }
	Con_Mat_->Replace_Element(j,j,(*Con_Mat_)(j,j)+1.0/con_weights_[j]);
      }
    HDSA::Linear_Algebra::QR_Factorization<RealT>(Con_Mat_, Con_Mat_Q_, Con_Mat_R_);
  }

  std::vector<RealT> Set_z_cov(void) const
  {
    int dim = (nx_+1)*(ny_+1);
    std::vector<RealT> z_cov = std::vector<RealT>(dim,z_cov_scale_);
    return z_cov;
  }
  
  void Apply_K_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in) const 
  {
    HDSA::Ptr<ROL::Vector<RealT> > u_in_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(u_in)->get_rol_vec();
    HDSA::Ptr<ROL::Vector<RealT> > u_out_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(u_out)->get_rol_vec();
    HDSA::Ptr<ROL::Vector<RealT> > u_tmp = u_out_rol->clone();
    RealT tol = 1.e-8;
    con_->applyJacobian_1(*u_out_rol,*u_in_rol,*u_in_rol,*u_in_rol,tol);
    //mass_con_->applyInverseJacobian_1(*u_tmp,*u_out_rol,*u_in_rol,*u_in_rol,tol);
    //con_->applyJacobian_1(*u_out_rol,*u_tmp,*u_in_rol,*u_in_rol,tol);
  }

  void Apply_K_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in) const 
  {
    HDSA::Ptr<ROL::Vector<RealT> > u_in_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(u_in)->get_rol_vec();
    HDSA::Ptr<ROL::Vector<RealT> > u_out_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(u_out)->get_rol_vec();
    HDSA::Ptr<ROL::Vector<RealT> > u_tmp = u_out_rol->clone();
    RealT tol = 1.e-8;
    con_->applyInverseJacobian_1(*u_out_rol,*u_in_rol,*u_in_rol,*u_in_rol,tol);
    //mass_con_->applyJacobian_1(*u_tmp,*u_out_rol,*u_in_rol,*u_in_rol,tol);
    //con_->applyInverseJacobian_1(*u_out_rol,*u_tmp,*u_in_rol,*u_in_rol,tol);
  }

  void Apply_L_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in) const 
  {
    Apply_K_Mat(u_out,u_in);
    for(int k = 0; k < nx_+1; k++)
      {
    	u_out->Replace_Element(k,(*u_out)(k)+con_weights_[k]*(*u_in)(k));
      }
  }

  void Apply_L_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in) const 
  {
   
    RealT tol = parlist_->sublist("Problem").get("L inverse Tolerance", 1.e-8);
    std::string solver = parlist_->sublist("Problem").get("L inverse Solver", "CG");
    bool verbose = parlist_->sublist("Problem").get("L inverse verbose", false);
    HDSA::Ptr<HDSA::Linear_Operator<RealT> > A = HDSA::makePtr<L_Mat<RealT> >(*this);
    HDSA::Linear_Algebra::Iterative_Linear_Solve<RealT>(*u_out, *u_in, A, tol, solver, verbose);

    // HDSA::Ptr<HDSA::Vector<RealT> > u_tmp = u_out->Clone();
    // Apply_K_Mat_Inverse(u_tmp,u_in);

    // HDSA::Ptr<HDSA::Dense_Matrix<RealT> > v = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nx_+1,1);
    // HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Qv = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nx_+1,1);
    // for(int k = 0; k < nx_+1; k++)
    //   {
    // 	v->Replace_Element(k,0,(*u_tmp)(k));
    //   }
    // Con_Mat_Q_->Multiply(Qv,v,true,false);
    // HDSA::Ptr<HDSA::Vector<RealT> > qv = HDSA::makePtr<Std_Vector<RealT> >(nx_+1);
    // Qv->Write_Column_to_Vector(0,qv);
    // HDSA::Ptr<HDSA::Vector<RealT> > x = HDSA::makePtr<Std_Vector<RealT> >(nx_+1);
    // HDSA::Linear_Algebra::Upper_Tri_Solve<RealT>( x, qv, Con_Mat_R_);
    
    // u_tmp->set(*u_in);
    // for(int k = 0; k < nx_+1; k++)
    //   {
    // 	u_tmp->Replace_Element(k,(*u_in)(k)-(*x)(k));
    //   }

    // Apply_K_Mat_Inverse(u_out,u_tmp);
  }

  // Overload HDSA::Linear_Operator to take matrix vector products
  template <class ScalarType>
  class L_Mat : public HDSA::Linear_Operator<ScalarType>
  {
    Model_Error_Objects_CDR<ScalarType> model_error_objects_;
      
    public:
      
    L_Mat(const Model_Error_Objects_CDR<ScalarType> & model_error_objects): model_error_objects_(model_error_objects)
      { }
      
      //! Dtor
      ~L_Mat()
      {}
      
      void matvec(HDSA::Ptr<HDSA::Vector<ScalarType> > & y, const HDSA::Ptr<HDSA::Vector<ScalarType> > & x) const 
      {
        model_error_objects_.Apply_L_Mat(y,x);  
      }
      
  };


};


#endif
