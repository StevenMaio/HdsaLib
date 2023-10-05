#ifndef MODEL_ERROR_OBJECTS_THERMAL_FLUIDS_HPP
#define MODEL_ERROR_OBJECTS_THERMAL_FLUIDS_HPP

// Instantiation of Model_Error_Objects

#include "../../../PDE-OPT/TOOLS/linearpdeconstraint.hpp"
#include "elliptic_operator.hpp"
#include "control_mass_mat.hpp"
#include "div_operator.hpp"
#include "pde_thermal-fluids_weighted.hpp"

template <class RealT>
class Model_Error_Objects_thermal_fluids : public HDSA::Model_Error_Objects<RealT> {

private:
  HDSA::Ptr<HDSA::ParameterList> parlist_;
  HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > con_state_mass_;
  HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > con_gamma_;
  HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > con_div_;
  HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > con_state_res_weighted_;
  RealT z_var_;
  RealT z_var_smooth_;
  RealT div_penalty_;
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Gamma_Inv_;
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Gamma_Inv_Q_;
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Gamma_Inv_R_;

public:

  Model_Error_Objects_thermal_fluids(const HDSA::Ptr<HDSA::ParameterList> & parlist_sensitivity, const HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > & OP_Objects_Factory, 
			      const HDSA::Ptr<HDSA::Weight_Matrices<RealT> > & weight_matrices, const HDSA::Ptr<HDSA::ParameterList> & parlist):
    HDSA::Model_Error_Objects<RealT>(parlist_sensitivity,OP_Objects_Factory,weight_matrices), parlist_(parlist)
  { 
    z_var_ = parlist->sublist("Problem").get("Control Variance", 1.0);
    z_var_smooth_ = parlist->sublist("Problem").get("Control Smoothing", 1.e-6);
    div_penalty_ = parlist->sublist("Problem").get("Divergence Penalty", 1.0);
  }

  virtual ~Model_Error_Objects_thermal_fluids()
  { }

  void Construct_Objects(const HDSA::Ptr<const HDSA::Comm<int> > & comm) 
  {
    /*** Initialize main data structure. ***/
    HDSA::Ptr<MeshManager<RealT> > meshMgr = HDSA::makePtr<MeshManager_ThermalFluids<RealT> >(*parlist_);
   
    HDSA::Ptr<PDE<RealT> > pde_state_mass = HDSA::makePtr<Elliptic_Operator<RealT> >(*parlist_,0.0);
    con_state_mass_ = HDSA::makePtr<Linear_PDE_Constraint<RealT> >(pde_state_mass,meshMgr,comm->Get_Teuchos_Communicator(),*parlist_);

    HDSA::Ptr<PDE<RealT> > pde_div = HDSA::makePtr<PDE_Div_Op<RealT> >(*parlist_);
    con_div_ = HDSA::makePtr<Linear_PDE_Constraint<RealT> >(pde_div,meshMgr,comm->Get_Teuchos_Communicator(),*parlist_);

    HDSA::Ptr<PDE<RealT> > pde_gamma = HDSA::makePtr<Control_Mass_Mat<RealT> >(*parlist_,1.0/(z_var_*z_var_),z_var_smooth_);
    con_gamma_ = HDSA::makePtr<Linear_PDE_Constraint<RealT> >(pde_gamma,meshMgr,comm->Get_Teuchos_Communicator(),*parlist_);

    std::vector<RealT> weights = std::vector<RealT>(4,1.0);
    weights[2] = div_penalty_;
    HDSA::Ptr<PDE<RealT> > pde_state_res = HDSA::makePtr<PDE_ThermalFluids_Weighted<RealT> >(*parlist_,weights);
    con_state_res_weighted_ = HDSA::makePtr<PDE_Constraint<RealT> >(pde_state_res,meshMgr,comm->Get_Teuchos_Communicator(),*parlist_);


    HDSA::Ptr<HDSA::Vector<RealT> > z_in = HDSA::Model_Error_Objects<RealT>::OP_Objects_->z->Clone();
    HDSA::Ptr<HDSA::Vector<RealT> > z_out = HDSA::Model_Error_Objects<RealT>::OP_Objects_->z->Clone();
    HDSA::Ptr<ROL::Vector<RealT> > z_in_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(z_in)->get_rol_vec();
    HDSA::Ptr<ROL::Vector<RealT> > z_out_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(z_out)->get_rol_vec();
    int n = z_in->Get_nonzero_dim();
    Gamma_Inv_ =  HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(n,n);
    Gamma_Inv_Q_ =  HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(n,n);
    Gamma_Inv_R_ =  HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(n,n);
    RealT tol = 1.e-8;
    for(int k = 0; k < n; k++)
      {
	z_in->basis(z_in->Get_map_reduced_to_full(k));
	z_out->zero();
	con_gamma_->applyJacobian_2(*z_out_rol,*z_in_rol,*z_in_rol,*z_in_rol,tol);
	z_out->Set_Zeros();
	Gamma_Inv_->Write_Vector_to_Column(k,z_out);
      }
    HDSA::Linear_Algebra::QR_Factorization<RealT>(Gamma_Inv_, Gamma_Inv_Q_, Gamma_Inv_R_);
  }

  void Apply_Gamma_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const 
  {
    int n = z_in->Get_nonzero_dim();
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > v = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(n,1);
    v->Write_Vector_to_Column(0,z_in);
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Qv = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(n,1);
    Gamma_Inv_Q_->Multiply(Qv,v,true,false);
    HDSA::Ptr<HDSA::Vector<RealT> > qv = HDSA::makePtr<Std_Vector<RealT> >(n);
    Qv->Write_Column_to_Vector(0,qv);
    HDSA::Ptr<HDSA::Vector<RealT> > x = HDSA::makePtr<Std_Vector<RealT> >(n);
    HDSA::Linear_Algebra::Upper_Tri_Solve<RealT>( x, qv, Gamma_Inv_R_);
    for(int k = 0; k < n; k++)
      {
        z_out->Replace_Element(z_out->Get_map_reduced_to_full(k),(*x)(k));
      }
  }

  void Apply_Gamma_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const 
  {
    int n = z_in->Get_nonzero_dim();
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > v = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(n,1);
    v->Write_Vector_to_Column(0,z_in);
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Av = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(n,1);
    Gamma_Inv_->Multiply(Av,v,true,false);
    for(int k = 0; k < n; k++)
      {
        z_out->Replace_Element(z_out->Get_map_reduced_to_full(k),(*Av)(k,0));
      }
  }
  
  void Apply_L_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in) const 
  {
    HDSA::Ptr<ROL::Vector<RealT> > u_in_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(u_in)->get_rol_vec();
    HDSA::Ptr<ROL::Vector<RealT> > u_out_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(u_out)->get_rol_vec();
    HDSA::Ptr<HDSA::Vector<RealT> > u_tmp1 = u_out->Clone();
    HDSA::Ptr<ROL::Vector<RealT> > u_tmp1_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(u_tmp1)->get_rol_vec();
    HDSA::Ptr<HDSA::Vector<RealT> > u_tmp2 = u_out->Clone();
    HDSA::Ptr<ROL::Vector<RealT> > u_tmp2_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(u_tmp2)->get_rol_vec();
    RealT tol = 1.e-8;
    u_out_rol->zero();

    HDSA::Ptr<HDSA::Vector<RealT> > u_nom = u_out->Clone();
    u_nom->set(*HDSA::Model_Error_Objects<RealT>::OP_Objects_->u);
    HDSA::Ptr<ROL::Vector<RealT> > u_nom_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(u_nom)->get_rol_vec();
    HDSA::Ptr<HDSA::Vector<RealT> > z_nom = HDSA::Model_Error_Objects<RealT>::OP_Objects_->z->Clone();
    z_nom->set(*HDSA::Model_Error_Objects<RealT>::OP_Objects_->z);
    HDSA::Ptr<ROL::Vector<RealT> > z_nom_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(z_nom)->get_rol_vec();

    con_state_res_weighted_->applyJacobian_1(*u_tmp1_rol,*u_in_rol,*u_nom_rol,*z_nom_rol,tol);
    con_state_mass_->applyInverseJacobian_1(*u_tmp2_rol,*u_tmp1_rol,*u_nom_rol,*z_nom_rol,tol);
    con_state_res_weighted_->applyAdjointJacobian_1(*u_out_rol,*u_tmp2_rol,*u_nom_rol,*z_nom_rol,tol);
  }

  void Apply_L_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in) const 
  {
    HDSA::Ptr<ROL::Vector<RealT> > u_in_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(u_in)->get_rol_vec();
    HDSA::Ptr<ROL::Vector<RealT> > u_out_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(u_out)->get_rol_vec();
    HDSA::Ptr<HDSA::Vector<RealT> > u_tmp1 = u_out->Clone();
    HDSA::Ptr<ROL::Vector<RealT> > u_tmp1_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(u_tmp1)->get_rol_vec();
    HDSA::Ptr<HDSA::Vector<RealT> > u_tmp2 = u_out->Clone();
    HDSA::Ptr<ROL::Vector<RealT> > u_tmp2_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(u_tmp2)->get_rol_vec();
    RealT tol = 1.e-8;
    u_out_rol->zero();

    HDSA::Ptr<HDSA::Vector<RealT> > u_nom = u_out->Clone();
    u_nom->set(*HDSA::Model_Error_Objects<RealT>::OP_Objects_->u);
    HDSA::Ptr<ROL::Vector<RealT> > u_nom_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(u_nom)->get_rol_vec();
    HDSA::Ptr<HDSA::Vector<RealT> > z_nom = HDSA::Model_Error_Objects<RealT>::OP_Objects_->z->Clone();
    z_nom->set(*HDSA::Model_Error_Objects<RealT>::OP_Objects_->z);
    HDSA::Ptr<ROL::Vector<RealT> > z_nom_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(z_nom)->get_rol_vec();

    con_state_res_weighted_->applyInverseAdjointJacobian_1(*u_tmp1_rol,*u_in_rol,*u_nom_rol,*z_nom_rol,tol);
    con_state_mass_->applyJacobian_1(*u_tmp2_rol,*u_tmp1_rol,*u_nom_rol,*z_nom_rol,tol);
    con_state_res_weighted_->applyInverseJacobian_1(*u_out_rol,*u_tmp2_rol,*u_nom_rol,*z_nom_rol,tol);
  }

  void Optional_Postprocessing(HDSA::Ptr<HDSA::Vector<RealT> > S, HDSA::Ptr<HDSA::Dense_Matrix<RealT> > U, 
			       HDSA::Ptr<HDSA::Model_Error_Kronecker_Matrix<RealT> > V, HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Mz_V) const 
  { 
    HDSA::Ptr<HDSA::Vector<RealT> > u_vec1 = HDSA::Model_Error_Objects<RealT>::OP_Objects_->u->Clone();
    HDSA::Ptr<HDSA::Vector<RealT> > u_vec2 = HDSA::Model_Error_Objects<RealT>::OP_Objects_->u->Clone();
    HDSA::Ptr<ROL::Vector<RealT> > u_vec1_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(u_vec1)->get_rol_vec();
    HDSA::Ptr<ROL::Vector<RealT> > u_vec2_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(u_vec2)->get_rol_vec();
    RealT tol = 1.e-8;
    std::string name = "name_here.txt";
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > A = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(u_vec1->dimension(),S->dimension());

    for(int k = 0; k < S->dimension(); k++)
      {
	u_vec1->zero(); u_vec2->zero();
	V->uk->Write_Column_to_Vector(k, u_vec1);
	con_div_->applyJacobian_1(*u_vec2_rol,*u_vec1_rol,*u_vec1_rol,*u_vec1_rol,tol);
	A->Write_Vector_to_Column(k,u_vec2);
      }
    name = "divergence_matvec_with_uk.txt";
    A->Write_to_File(name);
  } 


};


#endif
