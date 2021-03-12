template <class RealT>
class Matlab_RS_Objective: public HDSA::RS_Objective<RealT> {

private:
  std::unique_ptr<matlab::engine::MATLABEngine> matlabPtr_;
  std::vector<matlab::data::Array> rs_obj_;

public:

  Matlab_RS_Objective(void)
  { 
    std::cout << "Starting Matlab session" << std::endl;
    matlabPtr_ = matlab::engine::startMATLAB(); 
    std::vector<matlab::data::Array> args;
    rs_obj_ = matlabPtr_->feval(u"Matlab_RS_Objective_Interface",1,args);
    std::cout << "Finished starting Matlab session" << std::endl;
  }

  virtual ~Matlab_RS_Objective()
  { }

  // evaluate objective function
  RealT value(const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true)
  { 
    // Need to map HDSA::Vector to matlab data type to pass into value
    long unsigned int z_dim = z.dimension();
    long unsigned int theta_dim = theta.dimension();
    matlab::data::ArrayFactory factory;
    matlab::data::TypedArray<RealT> Z = factory.createArray<RealT>({z_dim,1});
    for(int k = 0; k < z_dim; k++)
    {
     Z[k] = z(k);   
    }
    matlab::data::TypedArray<RealT> Theta = factory.createArray<RealT>({theta_dim,1});
    for(int k = 0; k < theta_dim; k++)
    {
     Theta[k] = theta(k);   
    }
    std::vector<matlab::data::Array> args;
    args.resize(4);
    args[0] = rs_obj_[0];
    args[1] = Z;
    args[2] = Theta;
    matlab::data::TypedArray<int> update_int = factory.createArray<int>({1,1});
    if(update)
    {
       update_int[0] = 1;
    }
    else
    {
       update_int[1] = 0;
    }
    args[3] = update_int;
    
    // Call value function in matlab
    std::vector<matlab::data::Array> result = matlabPtr_->feval(u"value",1,args);
    return result[0][0];
  }

  // evaluate the gradient with respect to z
  void gradient_z(HDSA::Vector<RealT> & grad, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const bool & update = true)
  {
    // Need to map HDSA::Vector to matlab data type to pass into value
    long unsigned int z_dim = z.dimension();
    long unsigned int theta_dim = theta.dimension();
    matlab::data::ArrayFactory factory;
    matlab::data::TypedArray<RealT> Z = factory.createArray<RealT>({z_dim,1});
    for(int k = 0; k < z_dim; k++)
    {
     Z[k] = z(k);   
    }
    matlab::data::TypedArray<RealT> Theta = factory.createArray<RealT>({theta_dim,1});
    for(int k = 0; k < theta_dim; k++)
    {
     Theta[k] = theta(k);   
    }
    std::vector<matlab::data::Array> args;
    args.resize(4);
    args[0] = rs_obj_[0];
    args[1] = Z;
    args[2] = Theta;
    matlab::data::TypedArray<int> update_int = factory.createArray<int>({1,1});
    if(update)
    {
       update_int[0] = 1;
    }
    else
    {
       update_int[0] = 0;
    }
    args[3] = update_int;
    // Call gradient_z function in matlab
    std::vector<matlab::data::Array> result = matlabPtr_->feval(u"gradient_z",1,args);
    // Map data from result to grad
    for(int k = 0; k < z_dim; k++)
    {
        grad.Replace_Element(k,result[0][k]);
    }
  }

  #if Is_hessVec_z_z_Implemented
  // evaluate the z,z hessian vector product
  void hessVec_z_z(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta,
		   const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr)
  {
    // Need to map HDSA::Vector to matlab data type to pass into value
    long unsigned int z_dim = z.dimension();
    long unsigned int theta_dim = theta.dimension();
    matlab::data::ArrayFactory factory;
    matlab::data::TypedArray<RealT> V = factory.createArray<RealT>({z_dim,1});
    for(int k = 0; k < z_dim; k++)
      {
	V[k] = v(k);   
      }
    matlab::data::TypedArray<RealT> Z = factory.createArray<RealT>({z_dim,1});
    for(int k = 0; k < z_dim; k++)
      {
	Z[k] = z(k);   
      }
    matlab::data::TypedArray<RealT> Theta = factory.createArray<RealT>({theta_dim,1});
    for(int k = 0; k < theta_dim; k++)
      {
	Theta[k] = theta(k);   
      }
    std::vector<matlab::data::Array> args;
    args.resize(5);
    args[0] = rs_obj_[0];
    args[1] = V;
    args[2] = Z;
    args[3] = Theta;
    matlab::data::TypedArray<int> update_int = factory.createArray<int>({1,1});
    if(update)
      {
	update_int[0] = 1;
      }
    else
      {
	update_int[0] = 0;
      }
    args[4] = update_int;
    // Call hessVec_z_z function in matlab
    std::vector<matlab::data::Array> result = matlabPtr_->feval(u"hessVec_z_z",1,args);
    // Map data from result to grad
    for(int k = 0; k < z_dim; k++)
      {
	hv.Replace_Element(k,result[0][k]);
      }
  }
  #endif

  #if Is_hessVec_z_theta_Implemented
  // evaluate the z,theta hessian vector product
  void hessVec_z_theta(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta,
		       const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr)
  {
    // Need to map HDSA::Vector to matlab data type to pass into value
    long unsigned int z_dim = z.dimension();
    long unsigned int theta_dim = theta.dimension();
    matlab::data::ArrayFactory factory;
    matlab::data::TypedArray<RealT> V = factory.createArray<RealT>({theta_dim,1});
    for(int k = 0; k < theta_dim; k++)
      {
	V[k] = v(k);   
      }
    matlab::data::TypedArray<RealT> Z = factory.createArray<RealT>({z_dim,1});
    for(int k = 0; k < z_dim; k++)
      {
	Z[k] = z(k);   
      }
    matlab::data::TypedArray<RealT> Theta = factory.createArray<RealT>({theta_dim,1});
    for(int k = 0; k < theta_dim; k++)
      {
	Theta[k] = theta(k);   
      }
    std::vector<matlab::data::Array> args;
    args.resize(5);
    args[0] = rs_obj_[0];
    args[1] = V;
    args[2] = Z;
    args[3] = Theta;
    matlab::data::TypedArray<int> update_int = factory.createArray<int>({1,1});
    if(update)
      {
	update_int[0] = 1;
      }
    else
      {
	update_int[0] = 0;
      }
    args[4] = update_int;
    // Call hessVec_z_theta function in matlab
    std::vector<matlab::data::Array> result = matlabPtr_->feval(u"hessVec_z_theta",1,args);
    // Map data from result to grad
    for(int k = 0; k < z_dim; k++)
      {
	hv.Replace_Element(k,result[0][k]);
      }
  }
  #endif

  #if Is_hessVec_theta_z_Implemented
  // evaluate the theta,z hessian vector product
  void hessVec_theta_z(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, 
                       const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr)
  {
    // Need to map HDSA::Vector to matlab data type to pass into value
    long unsigned int z_dim = z.dimension();
    long unsigned int theta_dim = theta.dimension();
    matlab::data::ArrayFactory factory;
    matlab::data::TypedArray<RealT> V = factory.createArray<RealT>({z_dim,1});
    for(int k = 0; k < z_dim; k++)
      {
	V[k] = v(k);   
      }
    matlab::data::TypedArray<RealT> Z = factory.createArray<RealT>({z_dim,1});
    for(int k = 0; k < z_dim; k++)
      {
	Z[k] = z(k);   
      }
    matlab::data::TypedArray<RealT> Theta = factory.createArray<RealT>({theta_dim,1});
    for(int k = 0; k < theta_dim; k++)
      {
	Theta[k] = theta(k);   
      }
    std::vector<matlab::data::Array> args;
    args.resize(5);
    args[0] = rs_obj_[0];
    args[1] = V;
    args[2] = Z;
    args[3] = Theta;
    matlab::data::TypedArray<int> update_int = factory.createArray<int>({1,1});
    if(update)
      {
	update_int[0] = 1;
      }
    else
      {
	update_int[0] = 0;
      }
    args[4] = update_int;
    // Call hessVec_theta_z function in matlab
    std::vector<matlab::data::Array> result = matlabPtr_->feval(u"hessVec_theta_z",1,args);
    // Map data from result to grad
    for(int k = 0; k < theta_dim; k++)
      {
	hv.Replace_Element(k,result[0][k]);
      }
  }
  #endif

  #if Is_Misfit_hessVec_z_z_Implemented
  // evaluate the misfit z,z hessian vector product (for computed likelihood informed subspaces only)
  void Misfit_hessVec_z_z(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, 
			  const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr)
  { 
    // Need to map HDSA::Vector to matlab data type to pass into value
    long unsigned int z_dim = z.dimension();
    long unsigned int theta_dim = theta.dimension();
    matlab::data::ArrayFactory factory;
    matlab::data::TypedArray<RealT> V = factory.createArray<RealT>({z_dim,1});
    for(int k = 0; k < z_dim; k++)
      {
	V[k] = v(k);   
      }
    matlab::data::TypedArray<RealT> Z = factory.createArray<RealT>({z_dim,1});
    for(int k = 0; k < z_dim; k++)
      {
	Z[k] = z(k);   
      }
    matlab::data::TypedArray<RealT> Theta = factory.createArray<RealT>({theta_dim,1});
    for(int k = 0; k < theta_dim; k++)
      {
	Theta[k] = theta(k);   
      }
    std::vector<matlab::data::Array> args;
    args.resize(5);
    args[0] = rs_obj_[0];
    args[1] = V;
    args[2] = Z;
    args[3] = Theta;
    matlab::data::TypedArray<int> update_int = factory.createArray<int>({1,1});
    if(update)
      {
	update_int[0] = 1;
      }
    else
      {
	update_int[0] = 0;
      }
    args[4] = update_int;
    // Call Misfit_hessVec_z_z function in matlab
    std::vector<matlab::data::Array> result = matlabPtr_->feval(u"Misfit_hessVec_z_z",1,args);
    // Map data from result to grad
    for(int k = 0; k < z_dim; k++)
      {
	hv.Replace_Element(k,result[0][k]);
      }
  }
  #endif

  #if Is_Regularization_hessVec_z_z_Implemented
  // evaluate the regularization z,z hessian vector product (for computed likelihood informed subspaces only)
  void Regularization_hessVec_z_z(HDSA::Vector<RealT> & hv, const HDSA::Vector<RealT> & v, const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, 
				  const bool & update = true, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_at_input = HDSA::nullPtr)
  { 
    // Need to map HDSA::Vector to matlab data type to pass into value
    long unsigned int z_dim = z.dimension();
    long unsigned int theta_dim = theta.dimension();
    matlab::data::ArrayFactory factory;
    matlab::data::TypedArray<RealT> V = factory.createArray<RealT>({z_dim,1});
    for(int k = 0; k < z_dim; k++)
      {
	V[k] = v(k);   
      }
    matlab::data::TypedArray<RealT> Z = factory.createArray<RealT>({z_dim,1});
    for(int k = 0; k < z_dim; k++)
      {
	Z[k] = z(k);   
      }
    matlab::data::TypedArray<RealT> Theta = factory.createArray<RealT>({theta_dim,1});
    for(int k = 0; k < theta_dim; k++)
      {
	Theta[k] = theta(k);   
      }
    std::vector<matlab::data::Array> args;
    args.resize(5);
    args[0] = rs_obj_[0];
    args[1] = V;
    args[2] = Z;
    args[3] = Theta;
    matlab::data::TypedArray<int> update_int = factory.createArray<int>({1,1});
    if(update)
      {
	update_int[0] = 1;
      }
    else
      {
	update_int[0] = 0;
      }
    args[4] = update_int;
    // Call Regularization_hessVec_z_z function in matlab
    std::vector<matlab::data::Array> result = matlabPtr_->feval(u"Regularization_hessVec_z_z",1,args);
    // Map data from result to grad
    for(int k = 0; k < z_dim; k++)
      {
	hv.Replace_Element(k,result[0][k]);
      }
  }
  #endif

};
