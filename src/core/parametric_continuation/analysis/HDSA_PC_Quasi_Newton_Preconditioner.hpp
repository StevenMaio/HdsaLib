#ifndef HDSA_PC_QUASI_NEWTON_PRECONDITIONER_HPP
#define HDSA_PC_QUASI_NEWTON_PRECONDITIONER_HPP

namespace HDSA
{

  template <class RealT>
  class PC_Quasi_Newton_Preconditioner{

  private:
    RealT tau_;
    int param_current_data_step_;
    int block_current_data_step_;
    std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > s_;
    std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > y_;    
    std::vector<RealT> rho_;
    std::vector<HDSA::Ptr<HDSA::Dense_Matrix<RealT> > > Lr_;
    std::vector<HDSA::Ptr<HDSA::Dense_Matrix<RealT> > > Dr_;
    std::vector<HDSA::Ptr<HDSA::MultiVector<RealT> > > Pr_;
    std::vector<HDSA::Ptr<HDSA::MultiVector<RealT> > > Wr_;
    
  protected:

    // Overload this function if a better initialization is available
    virtual void Apply_Initial_Inverse_Hessian_Approximation(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const
    {
      z_out.set(z_in);
    }

  public:
    
    PC_Quasi_Newton_Preconditioner()
    {
      tau_ = 1.e-6;
    }

    virtual ~PC_Quasi_Newton_Preconditioner()
    { }

    void Set_N(const int & N)
    {
      s_.resize(N);
      y_.resize(N);
      rho_.resize(N);
      Lr_.resize(N);
      Dr_.resize(N);
      Pr_.resize(N);
      Wr_.resize(N);
      param_current_data_step_ = 0;
      block_current_data_step_ = 0;
    }

    void Add_Parametric_Quasi_Newton_Data(const HDSA::Vector<RealT> & s_k, const HDSA::Vector<RealT> & y_k)
    {
      s_[param_current_data_step_] = s_k.clone();
      s_[param_current_data_step_]->set(s_k);
      y_[param_current_data_step_] = y_k.clone();
      y_[param_current_data_step_]->set(y_k);
      rho_[param_current_data_step_] = 1.0/(s_k.dot(y_k));
      if(rho_[param_current_data_step_] < 0.0)
	{
	  std::cout << "Error: Negative Curvature Parameter" << std::endl;
	}
      param_current_data_step_ += 1;
    }
    
    void Add_Block_Quasi_Newton_Data(const std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > & P, const std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > & W)
    {
      int m = P.size();
      std::vector<RealT> D;
      std::vector<int> vecs_to_retain;
      for(int i = 0; i < m; i++)
	{
	  RealT D_val = P[i]->dot(*W[i]);
	  RealT val = P[i]->norm();
	  if( D_val > val*tau_)
	    {
	      vecs_to_retain.push_back(i);
	      D.push_back(D_val);
	    }
	}
      int num_to_retain = vecs_to_retain.size();
      Pr_[block_current_data_step_] = HDSA::makePtr<HDSA::MultiVector<RealT> >(num_to_retain,*P[0]);
      Wr_[block_current_data_step_] = HDSA::makePtr<HDSA::MultiVector<RealT> >(num_to_retain,*W[0]);
      Dr_[block_current_data_step_] = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(num_to_retain,1);
      for(int i = 0; i < num_to_retain; i++)
	{
	  (*Pr_[block_current_data_step_])[i]->set(*P[vecs_to_retain[i]]);
	  (*Wr_[block_current_data_step_])[i]->set(*W[vecs_to_retain[i]]);
	  Dr_[block_current_data_step_]->Replace_Element(i,0,D[i]);
	}
      
      block_current_data_step_ += 1;
    }

    void Apply_Inverse_Hessian_Approximation(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const
    {
      Apply_QN_Inverse_Hessian_Approximation(z_out,z_in,param_current_data_step_,block_current_data_step_);
    }

    void Apply_QN_Inverse_Hessian_Approximation(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in, const int & param_counter, const int & block_counter) const
    {
      if(param_counter == 0 & block_counter == 0)
	{
	  Apply_Initial_Inverse_Hessian_Approximation(z_out,z_in);
	}
      else if(param_counter == block_counter)
	{
	  RealT alpha = s_[param_counter-1]->dot(z_in);
	  HDSA::Ptr<HDSA::Vector<RealT> > tmp = z_out.clone();
	  tmp->set(z_in);
	  tmp->axpy(-rho_[param_counter-1]*alpha,*y_[param_counter-1]);
	  HDSA::Ptr<HDSA::Vector<RealT> > tmp_out = z_out.clone();
	  Apply_QN_Inverse_Hessian_Approximation(*tmp_out,*tmp,param_counter-1,block_counter);
	  z_out.set(*tmp_out);
	  RealT coeff = rho_[param_counter-1] * (alpha - y_[param_counter-1]->dot(*tmp_out));
	  z_out.axpy(coeff,*s_[param_counter-1]);
	}
      else
	{
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Pr_tmp = Pr_[block_counter-1]->MatVec(z_in);
	  int m = Pr_tmp->numRows();
	  HDSA::Ptr<HDSA::Vector<RealT> > tmp = z_out.clone();
	  tmp->set(z_in);
	  for(int i = 0; i < m; i++)
	    {
	      RealT val = (*Pr_tmp)(i,0)/(*Dr_[block_counter-1])(i,0);
	      tmp->axpy(-val,*(*Wr_[block_counter-1])[i]);
	    }
	  
	  Apply_QN_Inverse_Hessian_Approximation(z_out,*tmp,param_counter,block_counter-1);
	  
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Wr_tmp = Wr_[block_counter-1]->MatVec(z_out);
	  for(int i = 0; i < m; i++)
	    {
	      RealT val = (*Wr_tmp)(i,0)/(*Dr_[block_counter-1])(i,0);
	      z_out.axpy(-val,*(*Pr_[block_counter-1])[i]);
	      val = (*Pr_tmp)(i,0)/(*Dr_[block_counter-1])(i,0);
	      z_out.axpy(val,*(*Pr_[block_counter-1])[i]);
	    }
	}
    }
    
  
  };

}

#endif
