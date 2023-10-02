#ifndef HDSA_ROL_MD_INTERFACE_HDSA_TEST_PROBLEM_HPP
#define HDSA_ROL_MD_INTERFACE_HDSA_TEST_PROBLEM_HPP

template <class RealT>
class ROL_Model_Discrepancy_Interface_hdsa_test_problem : public HDSA::ROL_Model_Discrepancy_Interface<RealT> {

private:
  int m_;

public:

  ROL_Model_Discrepancy_Interface_hdsa_test_problem(ROL::Ptr<ROL::Objective_SimOpt<RealT> > & obj_simopt, ROL::Ptr<ROL::Constraint_SimOpt<RealT> > & con_simopt, 
						    ROL::Ptr<ROL::Vector<RealT> > & u, ROL::Ptr<ROL::Vector<RealT> > & z, int m): 
    HDSA::ROL_Model_Discrepancy_Interface<RealT>(obj_simopt,con_simopt,u,z) 
  {  
    m_ = m;
  }

  virtual ~ROL_Model_Discrepancy_Interface_hdsa_test_problem()
  { }


  void Apply_u_Elliptic_Operator_Inverse(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const 
  {
    u_out.set(u_in);
  }

  void Apply_u_Elliptic_Operator_Inverse_Transpose(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const 
  {
    u_out.set(u_in);
  }

  void Apply_u_Mass_Mat(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const 
  {
    u_out.set(u_in);
  }

  void Apply_u_Mass_Mat_Inverse(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const 
  {
    u_out.set(u_in);
  }

  void Apply_Gamma_Mat_Inverse(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const 
  {
    z_out.set(z_in);
  }

  HDSA::Ptr<HDSA::Vector<RealT> > Load_Optimal_u(void) const 
  {
    ROL::Ptr<std::vector<RealT> > u_ptr    = ROL::makePtr<std::vector<RealT> >(m_, 0.0);
    ROL::StdVector<RealT> u(u_ptr);
    ROL::Ptr<ROL::Vector<RealT> > up  = ROL::makePtrFromRef(u);
    for (int i=0;i<m_;i++){
      (*u_ptr)[i] = std::pow(1.0 + (double)i/double(m_-1),3.0);
    }
    HDSA::Ptr<HDSA::Vector<RealT> > u_opt = HDSA::makePtr<HDSA::ROL_Vector<RealT> >(*up);
    HDSA::ROL_Vector<RealT>& u_opt_rol = dynamic_cast<HDSA::ROL_Vector<RealT>&>(*u_opt);
    u_opt_rol.rol_vec->set(*up);
    return u_opt;
  }

  HDSA::Ptr<HDSA::Vector<RealT> > Load_Optimal_z(void) const 
  {
    ROL::Ptr<std::vector<RealT> > z_ptr    = ROL::makePtr<std::vector<RealT> >(m_, 0.0);
    ROL::StdVector<RealT> z(z_ptr);
    ROL::Ptr<ROL::Vector<RealT> > zp  = ROL::makePtrFromRef(z);
    for (int i=0;i<m_;i++){
      (*z_ptr)[i] = 1.0 + (double)i/double(m_-1);
    }
    HDSA::Ptr<HDSA::Vector<RealT> > z_opt = HDSA::makePtr<HDSA::ROL_Vector<RealT> >(*zp);
    HDSA::ROL_Vector<RealT>& z_opt_rol = dynamic_cast<HDSA::ROL_Vector<RealT>&>(*z_opt);
    z_opt_rol.rol_vec->set(*zp);
    return z_opt;
  }

  HDSA::Ptr<HDSA::MultiVector<RealT> > Load_Z_Data(void) const 
  {
    HDSA::Ptr<HDSA::Vector<RealT> > z_opt = Load_Optimal_z();
    HDSA::Ptr<HDSA::MultiVector<RealT> > Z = HDSA::makePtr<HDSA::MultiVector<RealT> >(1,*z_opt);
    (*Z)[0]->set(*z_opt);
    return Z;
  }

  HDSA::Ptr<HDSA::MultiVector<RealT> > Load_Y_Data(void) const 
  {
    ROL::Ptr<std::vector<RealT> > u_ptr    = ROL::makePtr<std::vector<RealT> >(m_, 0.0);
    ROL::StdVector<RealT> u(u_ptr);
    ROL::Ptr<ROL::Vector<RealT> > up  = ROL::makePtrFromRef(u);
    for (int i=0;i<m_;i++){
      (*u_ptr)[i] = std::pow(1.0+0.2*(double)i/double(m_-1),2.0);
    }
    HDSA::Ptr<HDSA::Vector<RealT> > u_hdsa = HDSA::makePtr<HDSA::ROL_Vector<RealT> >(*up);
    HDSA::ROL_Vector<RealT>& u_hdsa_rol = dynamic_cast<HDSA::ROL_Vector<RealT>&>(*u_hdsa);
    u_hdsa_rol.rol_vec->set(*up);
    HDSA::Ptr<HDSA::MultiVector<RealT> > Y = HDSA::makePtr<HDSA::MultiVector<RealT> >(1,*u_hdsa);
    (*Y)[0]->set(*u_hdsa);
    return Y;
  }

};

#endif


