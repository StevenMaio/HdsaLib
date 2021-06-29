#ifndef HDSA_MODEL_ERROR_OBJECTS_HPP
#define HDSA_MODEL_ERROR_OBJECTS_HPP

namespace HDSA
{

template <class RealT>
class Model_Error_Objects {

public:

  HDSA::Ptr<HDSA::ParameterList> parlist_sensitivity_;
  HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects_Factory_;
  std::vector<RealT> z_cov_;
  HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects_;
  int m_;
  int n_;
  HDSA::Ptr<HDSA::Weight_Matrices<RealT> > weight_matrices_factory_;
  HDSA::Ptr<HDSA::Weight_Matrices<RealT> > weight_matrices_;
  HDSA::Ptr<HDSA::Vector<RealT> > theta_pde_;
  HDSA::Ptr<HDSA::Vector<RealT> > Mz_star_;
  HDSA::Ptr<HDSA::Vector<RealT> > g_;
  HDSA::Ptr<HDSA::Vector<RealT> > gamma_inv_z_star_;
  HDSA::Ptr<HDSA::Vector<RealT> > Einv_Mz_star_;
  HDSA::Ptr<HDSA::Vector<RealT> > coeff_Linv_g_;
  HDSA::Ptr<HDSA::Vector<RealT> > N_min_Einv_Mz_star_;
  RealT z_star_gamma_inv_z_star_;
  RealT coeff_g_Linv_g_;
  RealT zstar_Mz_zstar_;
  

  Model_Error_Objects(const HDSA::Ptr<HDSA::ParameterList> & parlist_sensitivity, const HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > & OP_Objects_Factory, 
		      const HDSA::Ptr<HDSA::Weight_Matrices<RealT> > & weight_matrices_factory):
    parlist_sensitivity_(parlist_sensitivity), OP_Objects_Factory_(OP_Objects_Factory), weight_matrices_factory_(weight_matrices_factory)
  {  }

  virtual ~Model_Error_Objects()
  { }

  virtual void Construct_Objects(const HDSA::Ptr<const HDSA::Comm<int> > & comm)
  { }

  virtual std::vector<RealT> Set_z_cov(void) const =0;

  virtual void Apply_L_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in) const = 0;

  virtual void Apply_L_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in) const = 0;
 
  void Construct_Model_Error_Objects_Test(void)
  {
    HDSA::Ptr<HDSA::Vector<RealT> > v_in = OP_Objects_->u->Clone();
    HDSA::Ptr<HDSA::Vector<RealT> > v_out = OP_Objects_->u->Clone();
    HDSA::Ptr<HDSA::Vector<RealT> > z_in = OP_Objects_->z->Clone();
    HDSA::Ptr<HDSA::Vector<RealT> > z_out = OP_Objects_->z->Clone();

    // Test L^{-1}
    v_in->randomize();
    std::string name1 = "v_in.txt";
    v_in->Write_to_File(name1);
    HDSA::Ptr<HDSA::Vector<RealT> > v = OP_Objects_->u->Clone();
    Apply_L_Mat_Inverse(v,v_in);
    std::string name2 = "v.txt";
    v->Write_to_File(name2);
    Apply_L_Mat(v_out,v);
    std::string name3 = "v_out.txt";
    v_out->Write_to_File(name3);

    v_out->axpy(-1.0,*v_in);
    std::cout << "Norm of L*L^{-1}*v - v = " << v_out->norm() << std::endl;

    std::string name;
    std::ofstream fout;

    // State weighting matrix
    std::vector<std::vector<RealT> > L;
    L.resize(m_);
    for(int i = 0; i < m_; i++)
      {
	L[i].resize(m_);
      }
    for(int j = 0; j < m_; j++)
      {
	std::cout << "Computing column " << j+1 << " out of " << m_ << " for L matrix." << std::endl;
	v_in->basis(j);
	v_out->zero();
	Apply_L_Mat(v_out,v_in);
	for(int i = 0; i < m_; i++)
	  {
	    L[i][j] = (*v_out)(i);
	  }
      }

      // Write Solutions to text files
      name = "L.txt";
      fout.open(name);
      for(int i = 0; i < m_; i++)
	{
	  for(int j = 0; j < m_; j++)
	    {
	      fout << std::setprecision(16) << L[i][j] << "  ";
	    }
	  fout << "  " << std::endl;
	}
      fout.close(); 

    // Solution operator jacobian
    std::vector<std::vector<RealT> > J;
    J.resize(m_);
    for(int i = 0; i < m_; i++)
      {
	J[i].resize(n_);
      }
    for(int j = 0; j < n_; j++)
      {
	std::cout << "Computing column " << j+1 << " out of " << n_ << " for J matrix." << std::endl;
	z_in->basis(j);
	v_out->zero();
	Apply_Solution_Operator_z_Jacobian(*v_out,*z_in);
	for(int i = 0; i < m_; i++)
	  {
	    J[i][j] = (*v_out)(i);
	  }
      }

      // Write Solutions to text files
      name = "J.txt";
      fout.open(name);
      for(int i = 0; i < m_; i++)
	{
	  for(int j = 0; j < n_; j++)
	    {
	      fout << std::setprecision(16) << J[i][j] << "  ";
	    }
	  fout << "  " << std::endl;
	}
      fout.close(); 

    // Solution operator jacobian
    std::vector<std::vector<RealT> > Mz;
    Mz.resize(n_);
    for(int i = 0; i < n_; i++)
      {
	Mz[i].resize(n_);
      }
    for(int j = 0; j < n_; j++)
      {
	std::cout << "Computing column " << j+1 << " out of " << n_ << " for Mz matrix." << std::endl;
	z_in->basis(j);
	z_out->zero();
        weight_matrices_->Apply_z_Weight_Mat(z_out,z_in);
	for(int i = 0; i < n_; i++)
	  {
	    Mz[i][j] = (*z_out)(i);
	  }
      }

      // Write Solutions to text files
      name = "Mz.txt";
      fout.open(name);
      for(int i = 0; i < n_; i++)
	{
	  for(int j = 0; j < n_; j++)
	    {
	      fout << std::setprecision(16) << Mz[i][j] << "  ";
	    }
	  fout << "  " << std::endl;
	}
      fout.close(); 

    // Full space hessian
    std::vector<std::vector<RealT> > H;
    H.resize(m_);
    for(int i = 0; i < m_; i++)
      {
	H[i].resize(m_);
      }
    for(int j = 0; j < m_; j++)
      {
	std::cout << "Computing column " << j+1 << " out of " << m_ << " for H matrix." << std::endl;
	v_in->basis(j);
	v_out->zero();
	OP_Objects_->fs_obj->hessVec_u_u(*v_out,*v_in,*OP_Objects_->u,*OP_Objects_->z,*theta_pde_,false);
	for(int i = 0; i < m_; i++)
	  {
	    H[i][j] = (*v_out)(i);
	  }
      }

      // Write Solutions to text files
      name = "H.txt";
      fout.open(name);
      for(int i = 0; i < m_; i++)
	{
	  for(int j = 0; j < m_; j++)
	    {
	      fout << std::setprecision(16) << H[i][j] << "  ";
	    }
	  fout << "  " << std::endl;
	}
      fout.close(); 

      // z_cov
      // Write Solutions to text files
      name = "z_cov.txt";
      fout.open(name);
      for(int k = 0; k < n_; k++)
	{
	  fout << std::setprecision(16) << z_cov_[k] << "  ";
	}
      fout.close();

      // Write Solutions to text files
      name = "g.txt";
      fout.open(name);
      for(int k = 0; k < m_; k++)
	{
	  fout << std::setprecision(16) << (*g_)(k) << "  ";
	}
      fout.close();
  }

  void Instantiate_Objects(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<const HDSA::Comm<int> > & comm)
  {
    theta_pde_ = theta;
    OP_Objects_ = OP_Objects_Factory_->Construct_Opt_Problem_Objects(theta,comm);
    z_cov_ = Set_z_cov();
    m_ = OP_Objects_->u->dimension();
    n_ = z_cov_.size();
    weight_matrices_ = weight_matrices_factory_->Construct_Weight_Matrices(theta,comm);
    Construct_Objects(comm);
  }

  void Precompute_Model_Error_Objects_Data(void)
  {
    // Compute Mass matrix product with nominal solution
    Mz_star_ = OP_Objects_->z->Clone();
    weight_matrices_->Apply_z_Weight_Mat(Mz_star_,OP_Objects_->z);

    // Compute zstar^T*Mz*zstar
    zstar_Mz_zstar_ = Mz_star_->dot(*OP_Objects_->z);

    // Compute misfit gradient at nominal solution
    g_ = OP_Objects_->u->Clone();
    OP_Objects_->fs_obj->gradient_u(*g_,*OP_Objects_->u,*OP_Objects_->z,*theta_pde_,false);

    // Compute Gamma^{-1} z^star and z^star^T*Gamma^{-1}*z^star
    gamma_inv_z_star_ = OP_Objects_->z->Clone();
    for(int k = 0; k < n_; k++)
      {
        gamma_inv_z_star_->Replace_Element(k,(*OP_Objects_->z)(k)/z_cov_[k]);
      }
    z_star_gamma_inv_z_star_ = gamma_inv_z_star_->dot(*OP_Objects_->z);

    // Compute E^{-1}*M*z_star
    Einv_Mz_star_ = OP_Objects_->z->Clone();
    HDSA::Ptr<HDSA::Vector<RealT> > z_tmp1 = OP_Objects_->z->Clone();
    HDSA::Ptr<HDSA::Vector<RealT> > z_tmp2 = OP_Objects_->z->Clone();
    weight_matrices_->Apply_z_Weight_Mat_Inverse(z_tmp1,Mz_star_);
    for(int k = 0; k < n_; k++)
      {
	RealT val = z_tmp1->Get_Element(k)/z_cov_[k];
	z_tmp2->Replace_Element(k,val);
      }
    RealT val_coeff = -z_tmp2->dot(*OP_Objects_->z)/(1.0+z_star_gamma_inv_z_star_);
    z_tmp2->axpy(val_coeff,*gamma_inv_z_star_);
    weight_matrices_->Apply_z_Weight_Mat_Inverse(Einv_Mz_star_,z_tmp2);
    
    // Compute (1+beta)*g^T*L^{-1}*g
    coeff_Linv_g_ = OP_Objects_->u->Clone();
    Apply_L_Mat_Inverse(coeff_Linv_g_,g_);
    coeff_Linv_g_->scale(1.0+z_star_gamma_inv_z_star_);
    coeff_g_Linv_g_ = g_->dot(*coeff_Linv_g_);

    N_min_Einv_Mz_star_ = OP_Objects_->z->Clone();
    Apply_N(N_min_Einv_Mz_star_,Mz_star_);
    N_min_Einv_Mz_star_->axpy(-1.0,*Einv_Mz_star_);

  }

  void Apply_Solution_Operator_z_Jacobian(HDSA::Vector<RealT> & S_z, const HDSA::Vector<RealT> & z)
  {
    HDSA::Ptr<HDSA::Vector<RealT> > Jz_z = OP_Objects_->u->Clone();
    OP_Objects_->con->jacobian_z(*Jz_z, z, *OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta, false);
    Jz_z->scale(-1.0);
    OP_Objects_->con->jacobian_u_inverse(S_z,*Jz_z, *OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta, false);
  }

  void Apply_Solution_Operator_z_Jacobian_Transpose(HDSA::Vector<RealT> & S_u, const HDSA::Vector<RealT> & u)
  {
    HDSA::Ptr<HDSA::Vector<RealT> > Jzt_u = OP_Objects_->u->Clone();
    OP_Objects_->con->jacobian_u_adjoint_inverse(*Jzt_u, u, *OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta, false);
    Jzt_u->scale(-1.0);
    OP_Objects_->con->jacobian_z_adjoint(S_u, *Jzt_u, *OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta, false);
  }

  void Apply_N(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const
  {
    HDSA::Ptr<HDSA::Vector<RealT> > z_tmp = z_in->Clone();
    weight_matrices_->Apply_z_Weight_Mat_Inverse(z_tmp,z_in);
    
    RealT val_coeff = z_tmp->dot(*gamma_inv_z_star_);
    z_tmp->axpy(val_coeff,*OP_Objects_->z);

    val_coeff = z_tmp->dot(*gamma_inv_z_star_)/(1.0+z_star_gamma_inv_z_star_);
    for(int k = 0; k < n_; k++)
      {
	RealT val = (*z_tmp)(k)/(z_cov_[k]) - val_coeff*(*gamma_inv_z_star_)(k);
    	z_tmp->Replace_Element(k,val);
      }

    z_tmp->scale(1.0/(1.0+z_star_gamma_inv_z_star_));

    weight_matrices_->Apply_z_Weight_Mat_Inverse(z_out,z_tmp);
  }

};

}

#endif


