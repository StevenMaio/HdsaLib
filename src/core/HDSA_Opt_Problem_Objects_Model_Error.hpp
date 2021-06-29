#ifndef HDSA_OPT_PROBLEM_OBJECTS_MODEL_ERROR_HPP
#define HDSA_OPT_PROBLEM_OBJECTS_MODEL_ERROR_HPP

namespace HDSA
{

template <class RealT>
class Opt_Problem_Objects_Model_Error : public Opt_Problem_Objects<RealT>{

  HDSA::Ptr<HDSA::Model_Error_Objects<RealT> > model_error_objects_;

public:

  Opt_Problem_Objects_Model_Error(HDSA::Ptr<HDSA::Model_Error_Objects<RealT> > & model_error_objects): model_error_objects_(model_error_objects) 
  {}

  virtual ~Opt_Problem_Objects_Model_Error()
  { }

  HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > Construct_Opt_Problem_Objects(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<const HDSA::Comm<int> > & comm) 
  {
    HDSA::Ptr<HDSA::Opt_Problem_Objects_Model_Error<RealT> > OP_Objects_model_error = HDSA::makePtr<HDSA::Opt_Problem_Objects_Model_Error<RealT> >(model_error_objects_);
    OP_Objects_model_error->model_error_objects_->Instantiate_Objects(theta,comm);
    OP_Objects_model_error->u = OP_Objects_model_error->model_error_objects_->OP_Objects_->u;
    OP_Objects_model_error->z = OP_Objects_model_error->model_error_objects_->OP_Objects_->z;
    OP_Objects_model_error->fs_obj = OP_Objects_model_error->model_error_objects_->OP_Objects_->fs_obj;
    OP_Objects_model_error->rs_obj = OP_Objects_model_error->model_error_objects_->OP_Objects_->rs_obj;
    OP_Objects_model_error->theta = HDSA::makePtr<Std_Vector<RealT> >(1);
    return OP_Objects_model_error;
  }

  void Solve_Optimization_Problem() 
  {
    model_error_objects_->OP_Objects_->Solve_Optimization_Problem();
  }

  void Load_Optimal_Solution() 
  {
    model_error_objects_->OP_Objects_->Load_Optimal_Solution();
    bool enforce_z_zeros_ = model_error_objects_->parlist_sensitivity_->sublist("Formulation").get("Enforce z Zeros",false);
    if(enforce_z_zeros_)
      {
	model_error_objects_->OP_Objects_->z->Enforce_Zeros();
      }
    model_error_objects_->Precompute_Model_Error_Objects_Data();
  }

  void Write_Optimal_Solution() 
  {
    model_error_objects_->OP_Objects_->Write_Optimal_Solution();
  }

  void Construct_Model_Error_Objects_Test()
  {
    model_error_objects_->Construct_Model_Error_Objects_Test();
  }

  HDSA::Ptr<HDSA::Model_Error_Objects<RealT> > Get_Model_Error_Objects() const
  {
    return model_error_objects_;
  }

};

}

#endif
