#ifndef MODEL_ERROR_OBJECTS_CDR_HPP
#define MODEL_ERROR_OBJECTS_CDR_HPP

#include "elliptic_op.hpp"

// Instantiation of Model_Error_Objects

template <class RealT>
class Model_Error_Objects_CDR : public HDSA::Model_Error_Objects<RealT> {

private:
  HDSA::Ptr<HDSA::ParameterList> parlist_;
  HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > con_;
  bool constructed_elliptic_op_;
  RealT epsilon_;

public:

  Model_Error_Objects_CDR(const HDSA::Ptr<HDSA::ParameterList> & parlist_sensitivity, const HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > & OP_Objects_Factory, 
			      const HDSA::Ptr<HDSA::Weight_Matrices<RealT> > & weight_matrices, const HDSA::Ptr<HDSA::ParameterList> & parlist):
    HDSA::Model_Error_Objects<RealT>(parlist_sensitivity,OP_Objects_Factory,weight_matrices), parlist_(parlist)
  {
    constructed_elliptic_op_ = false;
    epsilon_ = parlist_sensitivity->sublist("Model Error").get("Smoothing Factor", 1.e-6);
  }

  virtual ~Model_Error_Objects_CDR()
  { }

  void Construct_Objects(const HDSA::Ptr<const HDSA::Comm<int> > & comm) 
  {
    // Initialize PDE
    HDSA::Ptr<MeshManager<RealT> > meshMgr = HDSA::makePtr<MeshManager_Rectangle<RealT> >(*parlist_);
    HDSA::Ptr<PDE<RealT> > pde = HDSA::makePtr<Elliptic_Op<RealT> >(*parlist_,epsilon_);
    con_ = HDSA::makePtr<Linear_PDE_Constraint<RealT> >(pde,meshMgr,comm->Get_Teuchos_Communicator(),*parlist_);
  }

  std::vector<RealT> Set_z_cov(void) const
  {
    int nx = parlist_->sublist("Geometry").get("NX", 10);
    int ny = parlist_->sublist("Geometry").get("NY", 10);
    int dim = (nx+1)*(ny+1);
    std::vector<RealT> z_cov = std::vector<RealT>(dim,1.0);
    return z_cov;
  }
  

  void Apply_L_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in) const 
  {
    HDSA::Ptr<ROL::Vector<RealT> > u_in_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(u_in)->get_rol_vec();
    HDSA::Ptr<ROL::Vector<RealT> > u_out_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(u_out)->get_rol_vec();
    RealT tol = 1.e-8;
    con_->applyJacobian_1(*u_out_rol,*u_in_rol,*u_in_rol,*u_in_rol,tol);
  }

  void Apply_L_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in) const 
  {
    HDSA::Ptr<ROL::Vector<RealT> > u_in_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(u_in)->get_rol_vec();
    HDSA::Ptr<ROL::Vector<RealT> > u_out_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(u_out)->get_rol_vec();
    RealT tol = 1.e-8;
    con_->applyInverseJacobian_1(*u_out_rol,*u_in_rol,*u_in_rol,*u_in_rol,tol);
  }

};


#endif
