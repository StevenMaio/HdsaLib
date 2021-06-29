#ifndef MODEL_ERROR_OBJECTS_THERMAL_FLUIDS_HPP
#define MODEL_ERROR_OBJECTS_THERMAL_FLUIDS_HPP

// Instantiation of Model_Error_Objects

#include "../../../PDE-OPT/TOOLS/linearpdeconstraint.hpp"
#include "elliptic_operator.hpp"

template <class RealT>
class Model_Error_Objects_thermal_fluids : public HDSA::Model_Error_Objects<RealT> {

private:
  HDSA::Ptr<HDSA::ParameterList> parlist_;
  RealT epsilon_;
  HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > con_elliptic_;
  RealT z_cov_scale_;
  bool res_penalty_;
  RealT norm_penalty_weight_;

public:

  Model_Error_Objects_thermal_fluids(const HDSA::Ptr<HDSA::ParameterList> & parlist_sensitivity, const HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > & OP_Objects_Factory, 
			      const HDSA::Ptr<HDSA::Weight_Matrices<RealT> > & weight_matrices, const HDSA::Ptr<HDSA::ParameterList> & parlist):
    HDSA::Model_Error_Objects<RealT>(parlist_sensitivity,OP_Objects_Factory,weight_matrices), parlist_(parlist)
  { 
    z_cov_scale_ = parlist->sublist("Problem").get("Control Variance", 1.0);
    res_penalty_ = parlist->sublist("Problem").get("Penalize Residual", false);
    norm_penalty_weight_ = parlist->sublist("Problem").get("Norm Penalty Weight", 1.0);
    epsilon_ = parlist_sensitivity->sublist("Problem").get("Epsilon", 1.e-6);
  }

  virtual ~Model_Error_Objects_thermal_fluids()
  { }

  void Construct_Objects(const HDSA::Ptr<const HDSA::Comm<int> > & comm) 
  {
    /*** Initialize main data structure. ***/
    HDSA::Ptr<MeshManager<RealT> > meshMgr = HDSA::makePtr<MeshManager_ThermalFluids<RealT> >(*parlist_);
    HDSA::Ptr<PDE<RealT> > pde = HDSA::makePtr<Elliptic_Operator<RealT> >(*parlist_,epsilon_);
    con_elliptic_ = HDSA::makePtr<Linear_PDE_Constraint<RealT> >(pde,meshMgr,comm->Get_Teuchos_Communicator(),*parlist_);
  }

  std::vector<RealT> Set_z_cov(void) const
  {
    int nx = parlist_->sublist("Geometry").get("NX", 10);
    int ny = parlist_->sublist("Geometry").get("NY", 10);
    int dim = 3*(2*nx+1)*(2*ny+1) + (nx+1)*(ny+1);
    std::vector<RealT> z_cov = std::vector<RealT>(dim,z_cov_scale_);
    return z_cov;
  }
  
  void Apply_L_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in) const 
  {
    HDSA::Ptr<ROL::Vector<RealT> > u_in_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(u_in)->get_rol_vec();
    HDSA::Ptr<ROL::Vector<RealT> > u_out_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(u_out)->get_rol_vec();
    HDSA::Ptr<HDSA::Vector<RealT> > u_tmp1 = u_out->Clone();
    HDSA::Ptr<ROL::Vector<RealT> > u_tmp1_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(u_tmp1)->get_rol_vec();
    RealT tol = 1.e-8;
    u_out_rol->zero();

    HDSA::Model_Error_Objects<RealT>::OP_Objects_->con->jacobian_u(*u_tmp1,*u_in,*HDSA::Model_Error_Objects<RealT>::OP_Objects_->u,
    								   *HDSA::Model_Error_Objects<RealT>::OP_Objects_->z,*HDSA::Model_Error_Objects<RealT>::OP_Objects_->theta,false);
    HDSA::Model_Error_Objects<RealT>::OP_Objects_->con->jacobian_u_adjoint(*u_out,*u_tmp1,*HDSA::Model_Error_Objects<RealT>::OP_Objects_->u,
    									   *HDSA::Model_Error_Objects<RealT>::OP_Objects_->z,*HDSA::Model_Error_Objects<RealT>::OP_Objects_->theta,false);
    u_out->scale(norm_penalty_weight_);

    con_elliptic_->applyJacobian_1(*u_tmp1_rol,*u_in_rol,*u_in_rol,*u_in_rol,tol);
    u_out_rol->plus(*u_tmp1_rol);
  }

  void Apply_L_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in) const 
  {
    RealT tol = parlist_->sublist("Problem").get("L inverse Tolerance", 1.e-8);
    std::string solver = parlist_->sublist("Problem").get("L inverse Solver", "CG");
    bool verbose = parlist_->sublist("Problem").get("L inverse verbose", false);
    HDSA::Ptr<HDSA::Linear_Operator<RealT> > A = HDSA::makePtr<L_Mat<RealT> >(*this);
    HDSA::Linear_Algebra::Iterative_Linear_Solve<RealT>(*u_out, *u_in, A, tol, solver, verbose);
  }

  // Overload HDSA::Linear_Operator to take matrix vector products
  template <class ScalarType>
  class L_Mat : public HDSA::Linear_Operator<ScalarType>
  {
    Model_Error_Objects_thermal_fluids<ScalarType> model_error_objects_;
      
    public:
      
    L_Mat(const Model_Error_Objects_thermal_fluids<ScalarType> & model_error_objects): model_error_objects_(model_error_objects)
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
