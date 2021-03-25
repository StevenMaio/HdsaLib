#ifndef HDSA_RS_OBJECTIVE_MODEL_ERROR_HPP
#define HDSA_RS_OBJECTIVE_MODEL_ERROR_HPP

namespace HDSA
{

template <class RealT>
class RS_Objective_Model_Error : public RS_Objective<RealT>{

  const HDSA::Ptr<HDSA::Model_Error_Objects<RealT> > model_error_objects_;

public:

  RS_Objective_Model_Error(const HDSA::Ptr<HDSA::Model_Error_Objects<RealT> > & model_error_objects):
    model_error_objects_(model_error_objects)
  {}

  virtual ~RS_Objective_Model_Error()
  { }

  // evaluate objective function
  RealT value(const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    RealT val = model_error_objects_->OP_Objects_->rs_obj->value(z,theta,update);
    return val;
  }

  // evaluate the gradient with respect to z
  void gradient_z(HDSA::Vector<RealT> & grad, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true) 
  {
    model_error_objects_->OP_Objects_->rs_obj->gradient_z(grad,z,theta,update);
  } 
 
  // evaluate the z,z hessian vector product
  void hessVec_z_z(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta,
		   const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr)
  {
    model_error_objects_->OP_Objects_->rs_obj->hessVec_z_z(hv,v,z,theta,update,grad_at_input);
  }

  // evaluate the z,theta hessian vector product
  void hessVec_z_theta(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta,
		       const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr)
  {
    HDSA::Ptr<HDSA::Vector<RealT> > u_vec_1 = model_error_objects_->OP_Objects_->u->Clone();
    model_error_objects_->Apply_A(*u_vec_1,v);
    HDSA::Ptr<HDSA::Vector<RealT> > u_vec_2 = u_vec_1->Clone();
    model_error_objects_->OP_Objects_->fs_obj->hessVec_u_u(*u_vec_2, *u_vec_1, *model_error_objects_->OP_Objects_->u, *model_error_objects_->OP_Objects_->z, *model_error_objects_->OP_Objects_->theta, false);
    model_error_objects_->Apply_Solution_Operator_z_Jacobian_Transpose(hv,*u_vec_2);
    HDSA::Ptr<HDSA::Vector<RealT> > z_vec_1 = z.Clone();
    model_error_objects_->Apply_X(*z_vec_1,v);
    hv.plus(*z_vec_1);
  }

  // evaluate the theta,z hessian vector product
  void hessVec_theta_z(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, 
		       const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr)
  {
    HDSA::Ptr<HDSA::Vector<RealT> > u_vec_1 = model_error_objects_->OP_Objects_->u->Clone();
    model_error_objects_->Apply_Solution_Operator_z_Jacobian(*u_vec_1,v);
    HDSA::Ptr<HDSA::Vector<RealT> > u_vec_2 = u_vec_1->Clone();
    model_error_objects_->OP_Objects_->fs_obj->hessVec_u_u(*u_vec_2, *u_vec_1, *model_error_objects_->OP_Objects_->u, *model_error_objects_->OP_Objects_->z, *model_error_objects_->OP_Objects_->theta, false);
    model_error_objects_->Apply_A_Transpose(hv,*u_vec_2);
    HDSA::Ptr<HDSA::Vector<RealT> > theta_vec_1 = hv.Clone();
    model_error_objects_->Apply_X_Transpose(*theta_vec_1,v);
    hv.plus(*theta_vec_1);
  }

  // evaluate the misfit z,z hessian vector product (for computed likelihood informed subspaces only)
  void Misfit_hessVec_z_z(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, 
			  const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr)
  { 
    model_error_objects_->OP_Objects_->rs_obj->Misfit_hessVec_z_z(hv,v,z,theta,update,grad_at_input);
  }

  // evaluate the regularization z,z hessian vector product (for computed likelihood informed subspaces only)
  void Regularization_hessVec_z_z(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, 
				  const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr)
  { 
    model_error_objects_->OP_Objects_->rs_obj->Regularization_hessVec_z_z(hv,v,z,theta,update,grad_at_input);
  }
 
 
};

}

#endif




// #ifndef HDSA_RS_OBJECTIVE_MODEL_ERROR_HPP
// #define HDSA_RS_OBJECTIVE_MODEL_ERROR_HPP

// namespace HDSA
// {

// template <class RealT>
// class RS_Objective_Model_Error : public RS_Objective<RealT>{

//   const HDSA::Ptr<HDSA::Model_Error_Objects<RealT> > model_error_objects_;

// public:

//   RS_Objective_Model_Error(const HDSA::Ptr<HDSA::Model_Error_Objects<RealT> > & model_error_objects):
//     model_error_objects_(model_error_objects)
//   {}

//   virtual ~RS_Objective_Model_Error()
//   { }

//   // evaluate objective function
//   RealT value(const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true)
//   {
//     RealT val = model_error_objects_->OP_Objects_->rs_obj->value(z,theta,update);
//     return val;
//   }

//   // evaluate the gradient with respect to z
//   void gradient_z(HDSA::Vector<RealT> & grad, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true) 
//   {
//     model_error_objects_->OP_Objects_->rs_obj->gradient_z(grad,z,theta,update);
//   } 
 
//   // evaluate the z,z hessian vector product
//   void hessVec_z_z(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta,
// 		   const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr)
//   {
//     model_error_objects_->OP_Objects_->rs_obj->hessVec_z_z(hv,v,z,theta,update,grad_at_input);
//   }

//   // evaluate the z,theta hessian vector product
//   void hessVec_z_theta(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta,
// 		       const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr)
//   {
//     HDSA::Ptr<HDSA::Vector<RealT> > u_vec_1 = model_error_objects_->OP_Objects_->u->Clone();
//     model_error_objects_->Apply_AY(*u_vec_1,v);
//     HDSA::Ptr<HDSA::Vector<RealT> > u_vec_2 = u_vec_1->Clone();
//     model_error_objects_->Apply_P(*u_vec_2,*u_vec_1);
//     model_error_objects_->OP_Objects_->fs_obj->hessVec_u_u(*u_vec_1, *u_vec_2, *model_error_objects_->OP_Objects_->u, *model_error_objects_->OP_Objects_->z, *model_error_objects_->OP_Objects_->theta, false);
//     model_error_objects_->Apply_Solution_Operator_z_Jacobian_Tranpose(hv,*u_vec_1);
//     HDSA::Ptr<HDSA::Vector<RealT> > z_vec_1 = z.Clone();
//     model_error_objects_->Apply_X(*z_vec_1,v);
//     hv.plus(*z_vec_1);
//   }

//   // evaluate the theta,z hessian vector product
//   void hessVec_theta_z(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, 
// 		       const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr)
//   {
//     HDSA::Ptr<HDSA::Vector<RealT> > u_vec_1 = model_error_objects_->OP_Objects_->u->Clone();
//     model_error_objects_->Apply_Solution_Operator_z_Jacobian(*u_vec_1,v);
//     HDSA::Ptr<HDSA::Vector<RealT> > u_vec_2 = u_vec_1->Clone();
//     model_error_objects_->OP_Objects_->fs_obj->hessVec_u_u(*u_vec_2, *u_vec_1, *model_error_objects_->OP_Objects_->u, *model_error_objects_->OP_Objects_->z, *model_error_objects_->OP_Objects_->theta, false);
//     model_error_objects_->Apply_P(*u_vec_1,*u_vec_2);
//     model_error_objects_->Apply_AY_Transpose(hv,*u_vec_1);
//     HDSA::Ptr<HDSA::Vector<RealT> > theta_vec_1 = hv.Clone();
//     model_error_objects_->Apply_X_Transpose(*theta_vec_1,v);
//     hv.plus(*theta_vec_1);
//   }

//   // evaluate the misfit z,z hessian vector product (for computed likelihood informed subspaces only)
//   void Misfit_hessVec_z_z(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, 
// 			  const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr)
//   { 
//     model_error_objects_->OP_Objects_->rs_obj->Misfit_hessVec_z_z(hv,v,z,theta,update,grad_at_input);
//   }

//   // evaluate the regularization z,z hessian vector product (for computed likelihood informed subspaces only)
//   void Regularization_hessVec_z_z(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, 
// 				  const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr)
//   { 
//     model_error_objects_->OP_Objects_->rs_obj->Regularization_hessVec_z_z(hv,v,z,theta,update,grad_at_input);
//   }
 
 
// };

// }

// #endif
