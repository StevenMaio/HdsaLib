#ifndef WEIGHT_MATRICES_CDR_HPP
#define WEIGHT_MATRICES_CDR_HPP

#include "elliptic_op.hpp"

template <class RealT>
class Weight_Matrices_CDR : public HDSA::Weight_Matrices<RealT> {

private:
  HDSA::Ptr<HDSA::ParameterList> parlist_;
  HDSA::Ptr<HDSA::ParameterList> parlist_sensitivity_;
  HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > con_;

public:

  Weight_Matrices_CDR(const HDSA::Ptr<HDSA::ParameterList> & parlist, const HDSA::Ptr<HDSA::ParameterList> & parlist_sensitivity):
    HDSA::Weight_Matrices<RealT>(parlist_sensitivity), parlist_(parlist), parlist_sensitivity_(parlist_sensitivity)
  { }

  Weight_Matrices_CDR(const HDSA::Ptr<HDSA::ParameterList> & parlist, const HDSA::Ptr<const HDSA::Comm<int> > & comm, const HDSA::Ptr<HDSA::ParameterList> & parlist_sensitivity):
    HDSA::Weight_Matrices<RealT>(parlist_sensitivity), parlist_(parlist), parlist_sensitivity_(parlist_sensitivity)
  { 
    HDSA::Ptr<MeshManager<RealT> > meshMgr = HDSA::makePtr<MeshManager_Rectangle<RealT> >(*parlist);
    HDSA::Ptr<PDE<RealT> > elliptic_pde = HDSA::makePtr<Elliptic_Op<RealT> >(*parlist,0.0);
    con_ = HDSA::makePtr<Linear_PDE_Constraint<RealT> >(elliptic_pde,meshMgr,comm->Get_Teuchos_Communicator(),*parlist);
  }

  HDSA::Ptr<HDSA::Weight_Matrices<RealT> > Construct_Weight_Matrices(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<const HDSA::Comm<int> > & comm) 
  {
    HDSA::Ptr<HDSA::Weight_Matrices<RealT> > weight_matrices = HDSA::makePtr<Weight_Matrices_CDR<RealT> >(parlist_,comm,parlist_sensitivity_);
    return weight_matrices;
  }

  void Apply_theta_Weight_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & theta_out, const HDSA::Ptr<HDSA::Vector<RealT> > & theta_in) const 
  {
    theta_out->set(*theta_in);
  }
    
  void Apply_z_Weight_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const
  {
    HDSA::Ptr<ROL::Vector<RealT> > z_in_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(z_in)->get_rol_vec();
    HDSA::Ptr<ROL::Vector<RealT> > z_out_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(z_out)->get_rol_vec();
    RealT tol = 1.e-8;
    con_->applyJacobian_1(*z_out_rol,*z_in_rol,*z_in_rol,*z_in_rol,tol);
  }

  void Apply_z_Weight_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const
  {
    HDSA::Ptr<ROL::Vector<RealT> > z_in_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(z_in)->get_rol_vec();
    HDSA::Ptr<ROL::Vector<RealT> > z_out_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(z_out)->get_rol_vec();
    RealT tol = 1.e-8;
    con_->applyInverseJacobian_1(*z_out_rol,*z_in_rol,*z_in_rol,*z_in_rol,tol);
  }
    
};


#endif
