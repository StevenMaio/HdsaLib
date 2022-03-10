#ifndef HDSA_BAYES_MODEL_ERROR_OBJECTS_HPP
#define HDSA_BAYES_MODEL_ERROR_OBJECTS_HPP

namespace HDSA
{

  template <class RealT>
  class Bayes_Model_Error_Objects: public HDSA::Model_Error_Objects<RealT>{

  private:

  public:
    
    Bayes_Model_Error_Objects(const HDSA::Ptr<HDSA::ParameterList> & parlist_sensitivity, const HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > & OP_Objects_Factory, 
			      const HDSA::Ptr<HDSA::Weight_Matrices<RealT> > & weight_matrices_factory): Model_Error_Objects<RealT>(parlist_sensitivity, OP_Objects_Factory, weight_matrices_factory)
    {
       
    }
    
    ~Bayes_Model_Error_Objects()
    { }
    
    virtual void Apply_Sqrt_Gamma_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const = 0;

    virtual void Apply_Sqrt_Gamma_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const = 0;
    
    virtual void Apply_Sqrt_L_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in) const = 0;

    virtual void Apply_L_Plus_Shift_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in, const RealT beta) const = 0;

    virtual void Apply_Sqrt_L_Plus_Shift_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in, const RealT beta) const = 0;

    virtual HDSA::Ptr<HDSA::MultiVector<RealT> > Load_Z_Data(void) const = 0;

    virtual HDSA::Ptr<HDSA::MultiVector<RealT> > Load_Y_Data(void) const = 0;

    void Construct_Model_Error_Objects_Test(void)
    {
      HDSA::Model_Error_Objects<RealT>::Construct_Model_Error_Objects_Test();

      HDSA::Ptr<HDSA::Vector<RealT> > v_in = HDSA::Model_Error_Objects<RealT>::OP_Objects_->u->Clone();
      HDSA::Ptr<HDSA::Vector<RealT> > v_out = HDSA::Model_Error_Objects<RealT>::OP_Objects_->u->Clone();
      HDSA::Ptr<HDSA::Vector<RealT> > z_in = HDSA::Model_Error_Objects<RealT>::OP_Objects_->z->Clone();
      HDSA::Ptr<HDSA::Vector<RealT> > z_out = HDSA::Model_Error_Objects<RealT>::OP_Objects_->z->Clone();

      std::string name;
      std::ofstream fout;
      // State weighting matrix
      std::vector<std::vector<RealT> > Gamma_Inv_Sqrt;
      Gamma_Inv_Sqrt.resize(HDSA::Model_Error_Objects<RealT>::n_);
      for(int i = 0; i < HDSA::Model_Error_Objects<RealT>::n_; i++)
	{
	  Gamma_Inv_Sqrt[i].resize(HDSA::Model_Error_Objects<RealT>::n_);
	}
      for(int j = 0; j < HDSA::Model_Error_Objects<RealT>::n_; j++)
	{
	  std::cout << "Computing column " << j+1 << " out of " << HDSA::Model_Error_Objects<RealT>::n_ << " for Gamma_Inv_Sqrt matrix." << std::endl;
	  z_in->basis(j);
	  z_in->Set_Zeros();
	  z_out->zero();
	  Apply_Sqrt_Gamma_Mat_Inverse(z_out,z_in);
	  for(int i = 0; i < HDSA::Model_Error_Objects<RealT>::n_; i++)
	    {
	      Gamma_Inv_Sqrt[i][j] = (*z_out)(i);
	    }
	}
      
      // Write Solutions to text files
      name = "Gamma_Inv_Sqrt.txt";
      fout.open(name);
      for(int i = 0; i < HDSA::Model_Error_Objects<RealT>::n_; i++)
	{
	  for(int j = 0; j < HDSA::Model_Error_Objects<RealT>::n_; j++)
	    {
	      fout << std::setprecision(16) << Gamma_Inv_Sqrt[i][j] << "  ";
	    }
	  fout << "  " << std::endl;
	}
      fout.close(); 

      std::vector<std::vector<RealT> > L_inv_Sqrt;
      L_inv_Sqrt.resize(HDSA::Model_Error_Objects<RealT>::m_);
      for(int i = 0; i < HDSA::Model_Error_Objects<RealT>::m_; i++)
	{
	  L_inv_Sqrt[i].resize(HDSA::Model_Error_Objects<RealT>::m_);
	}
      for(int j = 0; j < HDSA::Model_Error_Objects<RealT>::m_; j++)
	{
	  std::cout << "Computing column " << j+1 << " out of " << HDSA::Model_Error_Objects<RealT>::m_ << " for L_Inv_Sqrt matrix." << std::endl;
	  v_in->basis(j);
	  v_in->Set_Zeros();
	  v_out->zero();
	  Apply_Sqrt_L_Mat_Inverse(v_out,v_in);
	  for(int i = 0; i < HDSA::Model_Error_Objects<RealT>::m_; i++)
	    {
	      L_inv_Sqrt[i][j] = (*v_out)(i);
	    }
	}
      
      // Write Solutions to text files
      name = "L_Inv_Sqrt.txt";
      fout.open(name);
      for(int i = 0; i < HDSA::Model_Error_Objects<RealT>::m_; i++)
	{
	  for(int j = 0; j < HDSA::Model_Error_Objects<RealT>::m_; j++)
	    {
	      fout << std::setprecision(16) << L_inv_Sqrt[i][j] << "  ";
	    }
	  fout << "  " << std::endl;
	}
      fout.close(); 

      std::vector<std::vector<RealT> > L_shift_inv;
      L_shift_inv.resize(HDSA::Model_Error_Objects<RealT>::m_);
      for(int i = 0; i < HDSA::Model_Error_Objects<RealT>::m_; i++)
	{
	  L_shift_inv[i].resize(HDSA::Model_Error_Objects<RealT>::m_);
	}
      for(int j = 0; j < HDSA::Model_Error_Objects<RealT>::m_; j++)
	{
	  std::cout << "Computing column " << j+1 << " out of " << HDSA::Model_Error_Objects<RealT>::m_ << " for L_shift_inv matrix." << std::endl;
	  v_in->basis(j);
	  v_in->Set_Zeros();
	  v_out->zero();
	  Apply_L_Plus_Shift_Mat_Inverse(v_out,v_in,1.0);
	  for(int i = 0; i < HDSA::Model_Error_Objects<RealT>::m_; i++)
	    {
	      L_shift_inv[i][j] = (*v_out)(i);
	    }
	}
      
      // Write Solutions to text files
      name = "L_shift_inv.txt";
      fout.open(name);
      for(int i = 0; i < HDSA::Model_Error_Objects<RealT>::m_; i++)
	{
	  for(int j = 0; j < HDSA::Model_Error_Objects<RealT>::m_; j++)
	    {
	      fout << std::setprecision(16) << L_shift_inv[i][j] << "  ";
	    }
	  fout << "  " << std::endl;
	}
      fout.close(); 

      std::vector<std::vector<RealT> > L_shift_inv_Sqrt;
      L_shift_inv_Sqrt.resize(HDSA::Model_Error_Objects<RealT>::m_);
      for(int i = 0; i < HDSA::Model_Error_Objects<RealT>::m_; i++)
	{
	  L_shift_inv_Sqrt[i].resize(HDSA::Model_Error_Objects<RealT>::m_);
	}
      for(int j = 0; j < HDSA::Model_Error_Objects<RealT>::m_; j++)
	{
	  std::cout << "Computing column " << j+1 << " out of " << HDSA::Model_Error_Objects<RealT>::m_ << " for L_shift_inv_Sqrt matrix." << std::endl;
	  v_in->basis(j);
	  v_in->Set_Zeros();
	  v_out->zero();
	  Apply_Sqrt_L_Plus_Shift_Mat_Inverse(v_out,v_in,1.0);
	  for(int i = 0; i < HDSA::Model_Error_Objects<RealT>::m_; i++)
	    {
	      L_shift_inv_Sqrt[i][j] = (*v_out)(i);
	    }
	}
      
      // Write Solutions to text files
      name = "L_shift_inv_sqrt.txt";
      fout.open(name);
      for(int i = 0; i < HDSA::Model_Error_Objects<RealT>::m_; i++)
	{
	  for(int j = 0; j < HDSA::Model_Error_Objects<RealT>::m_; j++)
	    {
	      fout << std::setprecision(16) << L_shift_inv_Sqrt[i][j] << "  ";
	    }
	  fout << "  " << std::endl;
	}
      fout.close(); 

    }
    
  };

}

#endif
