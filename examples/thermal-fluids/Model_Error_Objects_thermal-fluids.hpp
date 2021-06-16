#ifndef MODEL_ERROR_OBJECTS_THERMAL_FLUIDS_HPP
#define MODEL_ERROR_OBJECTS_THERMAL_FLUIDS_HPP

// Instantiation of Model_Error_Objects

template <class RealT>
class Model_Error_Objects_thermal_fluids : public HDSA::Model_Error_Objects<RealT> {

private:
  HDSA::Ptr<HDSA::ParameterList> parlist_;
  
public:

  Model_Error_Objects_thermal_fluids(const HDSA::Ptr<HDSA::ParameterList> & parlist_sensitivity, const HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > & OP_Objects_Factory, 
			      const HDSA::Ptr<HDSA::Weight_Matrices<RealT> > & weight_matrices, const HDSA::Ptr<HDSA::ParameterList> & parlist):
    HDSA::Model_Error_Objects<RealT>(parlist_sensitivity,OP_Objects_Factory,weight_matrices), parlist_(parlist)
  { }

  virtual ~Model_Error_Objects_thermal_fluids()
  { }

  void Construct_Objects(const HDSA::Ptr<const HDSA::Comm<int> > & comm) 
  {
    // Need to implement
  }

  std::vector<RealT> Set_z_cov(void) const
  {
    int nx = parlist_->sublist("Geometry").get("NX", 10);
    int ny = parlist_->sublist("Geometry").get("NY", 10);
    int dim = 3*(2*nx+1)*(2*ny+1) + (nx+1)*(ny+1);
    std::vector<RealT> z_cov = std::vector<RealT>(dim,1.0);
    return z_cov;
  }
  
  void Apply_K_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in) const 
  {
    // Need to implement
    u_out->set(*u_in);
  }

  void Apply_K_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in) const 
  {
    // Need to implement
    u_out->set(*u_in);
  }

  void Apply_L_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in) const 
  {
    // Need to implement
    u_out->set(*u_in);
  }

  void Apply_L_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in) const 
  {
    // Need to implement
    u_out->set(*u_in);
  }

};


#endif
