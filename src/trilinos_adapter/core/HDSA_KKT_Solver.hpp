#ifndef HDSA_KKT_SOLVER_HPP
#define HDSA_KKT_SOLVER_HPP

#include "ROL_OptimizationProblem.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "ROL_Objective_SimOpt.hpp"
#include "ROL_Constraint_SimOpt.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"

namespace HDSA
{

namespace KKT_Solver
{

  template <class ScalarType>
  class HDSA_QP_Objective : public ROL::Objective_SimOpt<ScalarType> {
      
  private:
    HDSA::Ptr<HDSA::Opt_Problem_Objects<ScalarType> > OP_Objects_;
    HDSA::Ptr<HDSA::Vector<ScalarType> > rhs_u_;
    HDSA::Ptr<HDSA::Vector<ScalarType> > rhs_z_;
    HDSA::Ptr<HDSA::Vector<ScalarType> > rhs_lambda_;
    
  public:
    HDSA_QP_Objective(const HDSA::Ptr<HDSA::Opt_Problem_Objects<ScalarType> > & OP_Objects, const HDSA::Ptr<HDSA::Vector<ScalarType> > & rhs_u, 
		      const HDSA::Ptr<HDSA::Vector<ScalarType> > & rhs_z, const HDSA::Ptr<HDSA::Vector<ScalarType> > & rhs_lambda)
      : OP_Objects_(OP_Objects), rhs_u_(rhs_u), rhs_z_(rhs_z), rhs_lambda_(rhs_lambda)
    { }
    
    ScalarType value( const ROL::Vector<ScalarType> &u, const ROL::Vector<ScalarType> &z,  ScalarType &tol )
    {
      HDSA_Rol_Vector<ScalarType>* My_u;
      My_u = dynamic_cast<HDSA_Rol_Vector<ScalarType>*>(&const_cast<ROL::Vector<ScalarType> &>(u));
      HDSA_Rol_Vector<ScalarType>* My_z;
      My_z = dynamic_cast<HDSA_Rol_Vector<ScalarType>*>(&const_cast<ROL::Vector<ScalarType> &>(z));
      
      HDSA::Ptr<HDSA::Vector<ScalarType> > hv11 = OP_Objects_->u->Clone();
      OP_Objects_->fs_obj->hessVec_u_u(*hv11, *My_u->vec, *OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta,false);
      
      HDSA::Ptr<HDSA::Vector<ScalarType> > con_hv11 = OP_Objects_->u->Clone();
      OP_Objects_->con->hessian_u_u_adjoint(*con_hv11, *OP_Objects_->lambda, *My_u->vec, *OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta,false);	
      
      HDSA::Ptr<HDSA::Vector<ScalarType> > hv12 = OP_Objects_->u->Clone();
      OP_Objects_->fs_obj->hessVec_u_z(*hv12, *My_z->vec, *OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta,false);
      
      HDSA::Ptr<HDSA::Vector<ScalarType> > con_hv21 = OP_Objects_->u->Clone();
      OP_Objects_->con->hessian_u_z_adjoint(*con_hv21, *OP_Objects_->lambda, *My_z->vec, *OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta,false);	
      
      HDSA::Ptr<HDSA::Vector<ScalarType> > hv22 = OP_Objects_->z->Clone();
      OP_Objects_->fs_obj->hessVec_z_z(*hv22, *My_z->vec, *OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta,false);
      
      HDSA::Ptr<HDSA::Vector<ScalarType> > con_hv22 = OP_Objects_->z->Clone();
      OP_Objects_->con->hessian_z_z_adjoint(*con_hv22, *OP_Objects_->lambda, *My_z->vec, *OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta,false);	
      
      ScalarType val = 0.5*( hv11->dot(*My_u->vec) + con_hv11->dot(*My_u->vec) + hv22->dot(*My_z->vec) + con_hv22->dot(*My_z->vec) )
	+ hv12->dot(*My_u->vec) + con_hv21->dot(*My_u->vec) - rhs_u_->dot(*My_u->vec) - rhs_z_->dot(*My_z->vec);
      
      return val;
    }
    
    void gradient_1( ROL::Vector<ScalarType> &g, const ROL::Vector<ScalarType> &u, const ROL::Vector<ScalarType> &z, ScalarType &tol ) 
    {
      HDSA_Rol_Vector<ScalarType>* My_g;
      My_g = dynamic_cast<HDSA_Rol_Vector<ScalarType>*>(&const_cast<ROL::Vector<ScalarType> &>(g));
      HDSA_Rol_Vector<ScalarType>* My_u;
      My_u = dynamic_cast<HDSA_Rol_Vector<ScalarType>*>(&const_cast<ROL::Vector<ScalarType> &>(u));
      HDSA_Rol_Vector<ScalarType>* My_z;
      My_z = dynamic_cast<HDSA_Rol_Vector<ScalarType>*>(&const_cast<ROL::Vector<ScalarType> &>(z));
      
      HDSA::Ptr<HDSA::Vector<ScalarType> > hv11 = OP_Objects_->u->Clone();
      OP_Objects_->fs_obj->hessVec_u_u(*hv11, *My_u->vec, *OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta,false);
      
      HDSA::Ptr<HDSA::Vector<ScalarType> > con_hv11 = OP_Objects_->u->Clone();
      OP_Objects_->con->hessian_u_u_adjoint(*con_hv11, *OP_Objects_->lambda, *My_u->vec, *OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta,false);	
      
      HDSA::Ptr<HDSA::Vector<ScalarType> > hv12 = OP_Objects_->u->Clone();
      OP_Objects_->fs_obj->hessVec_u_z(*hv12, *My_z->vec, *OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta,false);
      
      HDSA::Ptr<HDSA::Vector<ScalarType> > con_hv21 = OP_Objects_->u->Clone();
      OP_Objects_->con->hessian_u_z_adjoint(*con_hv21, *OP_Objects_->lambda, *My_z->vec, *OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta,false);	
      
      My_g->vec->zero();
      My_g->vec->plus(*hv11);
      My_g->vec->plus(*con_hv11);
      My_g->vec->plus(*hv12);
      My_g->vec->plus(*con_hv21);
      My_g->vec->axpy(-1.0,*rhs_u_);
    }
    
    void gradient_2( ROL::Vector<ScalarType> &g, const ROL::Vector<ScalarType> &u, const ROL::Vector<ScalarType> &z, ScalarType &tol ) 
    {
      HDSA_Rol_Vector<ScalarType>* My_g;
      My_g = dynamic_cast<HDSA_Rol_Vector<ScalarType>*>(&const_cast<ROL::Vector<ScalarType> &>(g));
      HDSA_Rol_Vector<ScalarType>* My_u;
      My_u = dynamic_cast<HDSA_Rol_Vector<ScalarType>*>(&const_cast<ROL::Vector<ScalarType> &>(u));
      HDSA_Rol_Vector<ScalarType>* My_z;
      My_z = dynamic_cast<HDSA_Rol_Vector<ScalarType>*>(&const_cast<ROL::Vector<ScalarType> &>(z));
      
      HDSA::Ptr<HDSA::Vector<ScalarType> > hv22 = OP_Objects_->z->Clone();
      OP_Objects_->fs_obj->hessVec_z_z(*hv22, *My_z->vec, *OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta,false);
      
      HDSA::Ptr<HDSA::Vector<ScalarType> > con_hv22 = OP_Objects_->z->Clone();
      OP_Objects_->con->hessian_z_z_adjoint(*con_hv22, *OP_Objects_->lambda, *My_z->vec, *OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta,false);	
      
      HDSA::Ptr<HDSA::Vector<ScalarType> > hv21 = OP_Objects_->z->Clone();
      OP_Objects_->fs_obj->hessVec_z_u(*hv21, *My_u->vec, *OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta);
      
      HDSA::Ptr<HDSA::Vector<ScalarType> > con_hv12 = OP_Objects_->z->Clone();
      OP_Objects_->con->hessian_z_u_adjoint(*con_hv12, *OP_Objects_->lambda, *My_u->vec, *OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta,false);	
      
      My_g->vec->zero();
      My_g->vec->plus(*hv22);
      My_g->vec->plus(*con_hv22);
      My_g->vec->plus(*hv21);
      My_g->vec->plus(*con_hv12);
      My_g->vec->axpy(-1.0,*rhs_z_);
    }
    
    void hessVec_11( ROL::Vector<ScalarType> &hv, const ROL::Vector<ScalarType> &v, 
		     const ROL::Vector<ScalarType> &u,  const ROL::Vector<ScalarType> &z, ScalarType &tol ) 
    {
      HDSA_Rol_Vector<ScalarType>* My_hv;
      My_hv = dynamic_cast<HDSA_Rol_Vector<ScalarType>*>(&const_cast<ROL::Vector<ScalarType> &>(hv));
      HDSA_Rol_Vector<ScalarType>* My_v;
      My_v = dynamic_cast<HDSA_Rol_Vector<ScalarType>*>(&const_cast<ROL::Vector<ScalarType> &>(v));
      
      My_hv->vec->zero();
      OP_Objects_->fs_obj->hessVec_u_u(*My_hv->vec, *My_v->vec, *OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta,false);
      
      HDSA::Ptr<HDSA::Vector<ScalarType> > con_hv11 = OP_Objects_->u->Clone();
      OP_Objects_->con->hessian_u_u_adjoint(*con_hv11, *OP_Objects_->lambda, *My_v->vec, *OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta,false);	
      
      My_hv->vec->plus(*con_hv11);
    }
      
    
    void hessVec_12( ROL::Vector<ScalarType> &hv, const ROL::Vector<ScalarType> &v, 
		     const ROL::Vector<ScalarType> &u,  const ROL::Vector<ScalarType> &z, ScalarType &tol ) 
    {
      HDSA_Rol_Vector<ScalarType>* My_hv;
      My_hv = dynamic_cast<HDSA_Rol_Vector<ScalarType>*>(&const_cast<ROL::Vector<ScalarType> &>(hv));
      HDSA_Rol_Vector<ScalarType>* My_v;
      My_v = dynamic_cast<HDSA_Rol_Vector<ScalarType>*>(&const_cast<ROL::Vector<ScalarType> &>(v));
      
      My_hv->vec->zero();
      OP_Objects_->fs_obj->hessVec_u_z(*My_hv->vec, *My_v->vec, *OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta,false);

      HDSA::Ptr<HDSA::Vector<ScalarType> > con_hv21 = OP_Objects_->u->Clone();
      OP_Objects_->con->hessian_u_z_adjoint(*con_hv21, *OP_Objects_->lambda, *My_v->vec, *OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta,false);	
      
      My_hv->vec->plus(*con_hv21);	
    }

      
    void hessVec_21( ROL::Vector<ScalarType> &hv, const ROL::Vector<ScalarType> &v, 
		     const ROL::Vector<ScalarType> &u,  const ROL::Vector<ScalarType> &z, ScalarType &tol ) 
    {
      HDSA_Rol_Vector<ScalarType>* My_hv;
      My_hv = dynamic_cast<HDSA_Rol_Vector<ScalarType>*>(&const_cast<ROL::Vector<ScalarType> &>(hv));
      HDSA_Rol_Vector<ScalarType>* My_v;
      My_v = dynamic_cast<HDSA_Rol_Vector<ScalarType>*>(&const_cast<ROL::Vector<ScalarType> &>(v));
      
      My_hv->vec->zero();
      OP_Objects_->fs_obj->hessVec_z_u(*My_hv->vec, *My_v->vec, *OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta,false);
      
      HDSA::Ptr<HDSA::Vector<ScalarType> > con_hv12 = OP_Objects_->z->Clone();
      OP_Objects_->con->hessian_z_u_adjoint(*con_hv12, *OP_Objects_->lambda, *My_v->vec, *OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta,false);	
	
      My_hv->vec->plus(*con_hv12);
    }
      
    
    void hessVec_22( ROL::Vector<ScalarType> &hv, const ROL::Vector<ScalarType> &v, 
		     const ROL::Vector<ScalarType> &u,  const ROL::Vector<ScalarType> &z, ScalarType &tol ) 
    {
      HDSA_Rol_Vector<ScalarType>* My_hv;
      My_hv = dynamic_cast<HDSA_Rol_Vector<ScalarType>*>(&const_cast<ROL::Vector<ScalarType> &>(hv));
      HDSA_Rol_Vector<ScalarType>* My_v;
      My_v = dynamic_cast<HDSA_Rol_Vector<ScalarType>*>(&const_cast<ROL::Vector<ScalarType> &>(v));
      
      My_hv->vec->zero();
      OP_Objects_->fs_obj->hessVec_z_z(*My_hv->vec, *My_v->vec, *OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta,false);

      HDSA::Ptr<HDSA::Vector<ScalarType> > con_hv22 = OP_Objects_->z->Clone();
      OP_Objects_->con->hessian_z_z_adjoint(*con_hv22, *OP_Objects_->lambda, *My_v->vec, *OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta,false);	
      
      My_hv->vec->plus(*con_hv22);
    }
    
  };
  
  template <class ScalarType>
  class HDSA_QP_Constraint : public ROL::Constraint_SimOpt<ScalarType> {
    
  private:
    HDSA::Ptr<HDSA::Opt_Problem_Objects<ScalarType> > OP_Objects_;
    HDSA::Ptr<HDSA::Vector<ScalarType> > rhs_u_;
    HDSA::Ptr<HDSA::Vector<ScalarType> > rhs_z_;
    HDSA::Ptr<HDSA::Vector<ScalarType> > rhs_lambda_;
    
  public:
    HDSA_QP_Constraint(const HDSA::Ptr<HDSA::Opt_Problem_Objects<ScalarType> > & OP_Objects, const HDSA::Ptr<HDSA::Vector<ScalarType> > & rhs_u, 
		       const HDSA::Ptr<HDSA::Vector<ScalarType> > & rhs_z, const HDSA::Ptr<HDSA::Vector<ScalarType> > & rhs_lambda)
      : OP_Objects_(OP_Objects), rhs_u_(rhs_u), rhs_z_(rhs_z), rhs_lambda_(rhs_lambda)
    { }
    
    void value( ROL::Vector<ScalarType> &c, const ROL::Vector<ScalarType> &u, const ROL::Vector<ScalarType> &z,  ScalarType &tol ) 
    {
      HDSA_Rol_Vector<ScalarType>* My_c;
      My_c = dynamic_cast<HDSA_Rol_Vector<ScalarType>*>(&const_cast<ROL::Vector<ScalarType> &>(c));
      HDSA_Rol_Vector<ScalarType>* My_u;
      My_u = dynamic_cast<HDSA_Rol_Vector<ScalarType>*>(&const_cast<ROL::Vector<ScalarType> &>(u));
      HDSA_Rol_Vector<ScalarType>* My_z;
      My_z = dynamic_cast<HDSA_Rol_Vector<ScalarType>*>(&const_cast<ROL::Vector<ScalarType> &>(z));
      
      HDSA::Ptr<HDSA::Vector<ScalarType> > eu = OP_Objects_->lambda->Clone();
      OP_Objects_->con->jacobian_u(*eu,*My_u->vec,*OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta,false);
      
      HDSA::Ptr<HDSA::Vector<ScalarType> > ez = OP_Objects_->lambda->Clone();
      OP_Objects_->con->jacobian_z(*ez,*My_z->vec,*OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta,false);
      
      My_c->vec->set(*rhs_lambda_);
      My_c->vec->scale(-1.0);
      My_c->vec->plus(*eu);
      My_c->vec->plus(*ez);
    }
    
    void applyJacobian_1(ROL::Vector<ScalarType> &jv,
			 const ROL::Vector<ScalarType> &v,
			 const ROL::Vector<ScalarType> &u,
			 const ROL::Vector<ScalarType> &z,
			 ScalarType &tol) 
    {
      HDSA_Rol_Vector<ScalarType>* My_jv;
      My_jv = dynamic_cast<HDSA_Rol_Vector<ScalarType>*>(&const_cast<ROL::Vector<ScalarType> &>(jv));
      HDSA_Rol_Vector<ScalarType>* My_v;
      My_v = dynamic_cast<HDSA_Rol_Vector<ScalarType>*>(&const_cast<ROL::Vector<ScalarType> &>(v));
      
      My_jv->vec->zero();
      OP_Objects_->con->jacobian_u(*My_jv->vec,*My_v->vec,*OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta,false);
    }
  
    
    void applyJacobian_2(ROL::Vector<ScalarType> &jv,
			 const ROL::Vector<ScalarType> &v,
			 const ROL::Vector<ScalarType> &u,
			 const ROL::Vector<ScalarType> &z,
			 ScalarType &tol) 
    {
      HDSA_Rol_Vector<ScalarType>* My_jv;
      My_jv = dynamic_cast<HDSA_Rol_Vector<ScalarType>*>(&const_cast<ROL::Vector<ScalarType> &>(jv));
      HDSA_Rol_Vector<ScalarType>* My_v;
      My_v = dynamic_cast<HDSA_Rol_Vector<ScalarType>*>(&const_cast<ROL::Vector<ScalarType> &>(v));
      
      My_jv->vec->zero();
      OP_Objects_->con->jacobian_z(*My_jv->vec,*My_v->vec,*OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta,false); 
    }
    
    void applyAdjointJacobian_1(ROL::Vector<ScalarType> &ajv,
				const ROL::Vector<ScalarType> &v,
				const ROL::Vector<ScalarType> &u,
				const ROL::Vector<ScalarType> &z,
				ScalarType &tol) 
    {
      HDSA_Rol_Vector<ScalarType>* My_ajv;
      My_ajv = dynamic_cast<HDSA_Rol_Vector<ScalarType>*>(&const_cast<ROL::Vector<ScalarType> &>(ajv));
      HDSA_Rol_Vector<ScalarType>* My_v;
      My_v = dynamic_cast<HDSA_Rol_Vector<ScalarType>*>(&const_cast<ROL::Vector<ScalarType> &>(v));
      
      My_ajv->vec->zero();
      OP_Objects_->con->jacobian_u_adjoint(*My_ajv->vec,*My_v->vec,*OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta,false); 
    }
    
    void applyAdjointJacobian_2(ROL::Vector<ScalarType> &ajv,
				const ROL::Vector<ScalarType> &v,
				const ROL::Vector<ScalarType> &u,
				const ROL::Vector<ScalarType> &z,
				ScalarType &tol) 
    { 
      HDSA_Rol_Vector<ScalarType>* My_ajv;
      My_ajv = dynamic_cast<HDSA_Rol_Vector<ScalarType>*>(&const_cast<ROL::Vector<ScalarType> &>(ajv));
      HDSA_Rol_Vector<ScalarType>* My_v;
      My_v = dynamic_cast<HDSA_Rol_Vector<ScalarType>*>(&const_cast<ROL::Vector<ScalarType> &>(v));
      
      My_ajv->vec->zero();
      OP_Objects_->con->jacobian_z_adjoint(*My_ajv->vec,*My_v->vec,*OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta,false); 
    }
    
    void applyInverseJacobian_1(ROL::Vector<ScalarType> &ijv,
				const ROL::Vector<ScalarType> &v,
				const ROL::Vector<ScalarType> &u,
				const ROL::Vector<ScalarType> &z,
				ScalarType &tol) 
    {
      HDSA_Rol_Vector<ScalarType>* My_ijv;
      My_ijv = dynamic_cast<HDSA_Rol_Vector<ScalarType>*>(&const_cast<ROL::Vector<ScalarType> &>(ijv));
      HDSA_Rol_Vector<ScalarType>* My_v;
      My_v = dynamic_cast<HDSA_Rol_Vector<ScalarType>*>(&const_cast<ROL::Vector<ScalarType> &>(v));
      
      My_ijv->vec->zero();
      OP_Objects_->con->jacobian_u_inverse(*My_ijv->vec,*My_v->vec,*OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta,false); 
    }
    
    void applyInverseAdjointJacobian_1(ROL::Vector<ScalarType> &iajv,
				       const ROL::Vector<ScalarType> &v,
				       const ROL::Vector<ScalarType> &u,
				       const ROL::Vector<ScalarType> &z,
				       ScalarType &tol) 
    {
      HDSA_Rol_Vector<ScalarType>* My_iajv;
      My_iajv = dynamic_cast<HDSA_Rol_Vector<ScalarType>*>(&const_cast<ROL::Vector<ScalarType> &>(iajv));
      HDSA_Rol_Vector<ScalarType>* My_v;
      My_v = dynamic_cast<HDSA_Rol_Vector<ScalarType>*>(&const_cast<ROL::Vector<ScalarType> &>(v));
      
      My_iajv->vec->zero();
      OP_Objects_->con->jacobian_u_adjoint_inverse(*My_iajv->vec,*My_v->vec,*OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta,false); 
    }
    
    void applyAdjointHessian_11(ROL::Vector<ScalarType> &ahwv,
				const ROL::Vector<ScalarType> &w,
				const ROL::Vector<ScalarType> &v,
				const ROL::Vector<ScalarType> &u,
				const ROL::Vector<ScalarType> &z,
				ScalarType &tol)
    {
      ahwv.zero();
    }
    
    void applyAdjointHessian_12(ROL::Vector<ScalarType> &ahwv,
				const ROL::Vector<ScalarType> &w,
				const ROL::Vector<ScalarType> &v,
				const ROL::Vector<ScalarType> &u,
				const ROL::Vector<ScalarType> &z,
				ScalarType &tol)
    {
      ahwv.zero();
    }
    
    void applyAdjointHessian_21(ROL::Vector<ScalarType> &ahwv,
				const ROL::Vector<ScalarType> &w,
				const ROL::Vector<ScalarType> &v,
				const ROL::Vector<ScalarType> &u,
				const ROL::Vector<ScalarType> &z,
				ScalarType &tol)
    {
      ahwv.zero();
    }
    
    void applyAdjointHessian_22(ROL::Vector<ScalarType> &ahwv,
				const ROL::Vector<ScalarType> &w,
				const ROL::Vector<ScalarType> &v,
				const ROL::Vector<ScalarType> &u,
				const ROL::Vector<ScalarType> &z,
				ScalarType &tol)
    {
      ahwv.zero();
    }
    
  };

  template <class RealT>
  // Solve the KKT system
  void QP_Solve(const HDSA::Ptr<HDSA::Vector<RealT> > & sol_u, const HDSA::Ptr<HDSA::Vector<RealT> > & sol_z, const HDSA::Ptr<HDSA::Vector<RealT> > & sol_lambda,
		 const HDSA::Ptr<HDSA::Vector<RealT> > & rhs_u, const HDSA::Ptr<HDSA::Vector<RealT> > & rhs_z, const HDSA::Ptr<HDSA::Vector<RealT> > & rhs_lambda,
		 const HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > & OP_Objects, int maxits, RealT tol, HDSA::Ptr<std::ostream> & outStream)
  {
    // Instantiate the objective and constraint objects
    HDSA::Ptr<ROL::Objective_SimOpt<RealT> > obj = HDSA::makePtr<HDSA_QP_Objective<RealT> >(OP_Objects,rhs_u,rhs_z,rhs_lambda);
    HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > con = HDSA::makePtr<HDSA_QP_Constraint<RealT> >(OP_Objects,rhs_u,rhs_z,rhs_lambda);
      
    HDSA::Ptr<ROL::Vector<RealT> > up = HDSA::makePtr<HDSA_Rol_Vector<RealT> >(sol_u);
    HDSA::Ptr<ROL::Vector<RealT> > zp = HDSA::makePtr<HDSA_Rol_Vector<RealT> >(sol_z);
    HDSA::Ptr<ROL::Vector<RealT> > pp = HDSA::makePtr<HDSA_Rol_Vector<RealT> >(sol_lambda);
    HDSA::Ptr<ROL::Vector<RealT> > rp = HDSA::makePtr<HDSA_Rol_Vector<RealT> >(sol_lambda);
    
    // Instantiate the reduced objective
    HDSA::Ptr<ROL::Reduced_Objective_SimOpt<RealT> > robj = HDSA::makePtr<ROL::Reduced_Objective_SimOpt<RealT> >(obj, con, up, zp, pp);     
    
    Teuchos::ParameterList solver_parlist;
    Teuchos::ParameterList& general_list = solver_parlist.sublist("General");
    Teuchos::ParameterList& general_list_krylov = general_list.sublist("Krylov");
    general_list_krylov.set("Type","Conjugate Gradients");
    general_list_krylov.set("Relative Tolerance",tol);
    general_list_krylov.set("Iteration Limit",maxits);
    
    Teuchos::ParameterList& step_list = solver_parlist.sublist("Step");
    Teuchos::ParameterList& step_list_tr = step_list.sublist("Trust Region");
    step_list_tr.set("Subproblem Solver","Truncated CG");
    step_list_tr.set("Output Level",0);
    step_list_tr.set("Initial Radius",1.e14);
    step_list_tr.set("Maximum Radius",1.e16);
    
    Teuchos::ParameterList& status_test_list = solver_parlist.sublist("Status Test");
    status_test_list.set("Gradient Tolerance",tol);
    status_test_list.set("Constraint Tolerance",1.e-11);
    status_test_list.set("Step Tolerance",1.e-10);
    status_test_list.set("Iteration Limit",1);
    
    // Build optimization problem and check derivatives
    ROL::OptimizationProblem<RealT> optProb(robj,zp);
    // Build optimization solver and solve
    ROL::OptimizationSolver<RealT> optSolver(optProb,solver_parlist);
    optSolver.solve(*outStream);
    
    // Compute State
    RealT tols(1.e-8);
    con->solve(*rp,*up,*zp,tols);
    // Compute Adjoint
    HDSA::Ptr<ROL::Vector<RealT> > dualstate = up->dual().clone();
    dualstate->set(up->dual());
    obj->gradient_1(*dualstate,*up,*zp,tols);
    con->applyInverseAdjointJacobian_1(*pp,*dualstate,*up,*zp,tols);
    pp->scale(-1.0);
    
    sol_u->set(*HDSA::dynamicPtrCast<HDSA_Rol_Vector<RealT> >(up)->vec);
    sol_z->set(*HDSA::dynamicPtrCast<HDSA_Rol_Vector<RealT> >(zp)->vec);
    sol_lambda->set(*HDSA::dynamicPtrCast<HDSA_Rol_Vector<RealT> >(pp)->vec);
  }

}

}
    
#endif
