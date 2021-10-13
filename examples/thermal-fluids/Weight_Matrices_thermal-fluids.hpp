#ifndef WEIGHT_MATRICES_THERMAL_FLUIDS_HPP
#define WEIGHT_MATRICES_THERMAL_FLUIDS_HPP

#include "../../../PDE-OPT/TOOLS/linearpdeconstraint.hpp"
#include "control_mass_mat.hpp"

template <class RealT>
class Weight_Matrices_thermal_fluids : public HDSA::Weight_Matrices<RealT> {

private:
  HDSA::Ptr<HDSA::ParameterList> parlist_;
  HDSA::Ptr<HDSA::ParameterList> parlist_sensitivity_;
  HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > con_mass_;

public:

  Weight_Matrices_thermal_fluids(const HDSA::Ptr<HDSA::ParameterList> & parlist, const HDSA::Ptr<HDSA::ParameterList> & parlist_sensitivity):
    HDSA::Weight_Matrices<RealT>(parlist_sensitivity), parlist_(parlist), parlist_sensitivity_(parlist_sensitivity)
  {  }

  Weight_Matrices_thermal_fluids(const HDSA::Ptr<HDSA::ParameterList> & parlist, const HDSA::Ptr<HDSA::ParameterList> & parlist_sensitivity, const HDSA::Ptr<const HDSA::Comm<int> > & comm):
    HDSA::Weight_Matrices<RealT>(parlist_sensitivity), parlist_(parlist), parlist_sensitivity_(parlist_sensitivity)
  {  
    /*** Initialize main data structure. ***/
    HDSA::Ptr<MeshManager<RealT> > meshMgr = HDSA::makePtr<MeshManager_ThermalFluids<RealT> >(*parlist_);
    HDSA::Ptr<PDE<RealT> > pde = HDSA::makePtr<Control_Mass_Mat<RealT> >(*parlist_,1.0,0.0);
    con_mass_ = HDSA::makePtr<Linear_PDE_Constraint<RealT> >(pde,meshMgr,comm->Get_Teuchos_Communicator(),*parlist_);
  }

  HDSA::Ptr<HDSA::Weight_Matrices<RealT> > Construct_Weight_Matrices(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<const HDSA::Comm<int> > & comm) 
  {
    HDSA::Ptr<HDSA::Weight_Matrices<RealT> > weight_matrices = HDSA::makePtr<Weight_Matrices_thermal_fluids<RealT> >(parlist_,parlist_sensitivity_,comm);
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
    con_mass_->applyJacobian_2(*z_out_rol,*z_in_rol,*z_in_rol,*z_in_rol,tol);
    z_out->Set_Zeros();
  }

};


#endif
