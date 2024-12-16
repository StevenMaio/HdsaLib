#ifndef HDSA_PC_QUASI_NEWTON_PRECONDITIONER_HPP
#define HDSA_PC_QUASI_NEWTON_PRECONDITIONER_HPP

namespace HDSA
{

  template <class RealT>
  class PC_Quasi_Newton_Preconditioner{

  private:
    bool use_block_update_;
    RealT tau_;
    int max_storage_;
    int param_current_data_step_;
    int block_current_data_step_;
    std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > s_;
    std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > y_;    
    std::vector<RealT> rho_;
    std::vector<HDSA::Ptr<HDSA::Dense_Matrix<RealT> > > Dr_;
    std::vector<HDSA::Ptr<HDSA::MultiVector<RealT> > > Pr_;
    std::vector<HDSA::Ptr<HDSA::MultiVector<RealT> > > Wr_;
    
  protected:

    int total_vecs_stored_;

    // Overload this function if a better initialization is available
    virtual void Apply_Initial_Inverse_Hessian_Approximation(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const
    {
      z_out.set(z_in);
    }

  public:
    
    int getStep(void) const
    {
      return block_current_data_step_;
    }

    PC_Quasi_Newton_Preconditioner(bool use_block_update = true, RealT tau = 1.e-6, int max_storage = 10): use_block_update_(use_block_update), tau_(tau), max_storage_(max_storage)
    { }

    virtual ~PC_Quasi_Newton_Preconditioner()
    { }

    int Get_Number_Vecs_Stored(void)
    {
      return total_vecs_stored_;
    }

    void Set_N(const int & N)
    {
      s_.resize(N);
      y_.resize(N);
      rho_.resize(N);
      if(use_block_update_)
      {
        Dr_.resize(N);
        Pr_.resize(N);
        Wr_.resize(N);
        Pr_[0] = HDSA::makePtr<HDSA::MultiVector<RealT> >();
        Wr_[0] = HDSA::makePtr<HDSA::MultiVector<RealT> >();
      }
      param_current_data_step_ = 0;
      block_current_data_step_ = 0;
      total_vecs_stored_ = 0;
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
      total_vecs_stored_ += 2;
    }
    
    void Add_Block_Quasi_Newton_Step(const HDSA::Ptr<HDSA::Vector<RealT> > & p, const HDSA::Ptr<HDSA::Vector<RealT> > & w)
    {
      if(use_block_update_)
      {
        int num_vecs = Pr_[block_current_data_step_]->Number_of_Vectors() + 1;
        HDSA::Ptr<HDSA::MultiVector<RealT> >  Ptmp_ = HDSA::makePtr<HDSA::MultiVector<RealT> >(num_vecs,*p);
        HDSA::Ptr<HDSA::MultiVector<RealT> >  Wtmp_ = HDSA::makePtr<HDSA::MultiVector<RealT> >(num_vecs,*w);
        for(int k = 0; k < num_vecs-1; k++)
        {
          (*Ptmp_)[k]->set(*(*Pr_[block_current_data_step_])[k]);
          (*Wtmp_)[k]->set(*(*Wr_[block_current_data_step_])[k]);
        }

        RealT normalization = p->dot(*p);
        (*Ptmp_)[num_vecs-1]->set(*p);
        (*Ptmp_)[num_vecs-1]->scale(1.0/std::sqrt(normalization));
        (*Wtmp_)[num_vecs-1]->set(*w);
        (*Wtmp_)[num_vecs-1]->scale(1.0/std::sqrt(normalization));

        HDSA::Ptr<HDSA::Dense_Matrix<RealT> >  D = Ptmp_->MatMat(*Wtmp_);
        HDSA::Ptr<HDSA::Dense_Matrix<RealT> >  V = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(num_vecs,num_vecs);
        HDSA::Ptr<HDSA::Dense_Matrix<RealT> >  S = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(num_vecs,1);
        HDSA::Linear_Algebra::Symmetric_Eig_Decomposition<RealT>(*D,*V,*S);

        int modes_to_retain = 0;
        for(int k = 0; k < num_vecs; k++)
        {
         if( ((*S)(k,0) > tau_) & (modes_to_retain < max_storage_))
         {
          modes_to_retain += 1;
         } 
        }

        Pr_[block_current_data_step_]->Clear();
        Wr_[block_current_data_step_]->Clear();
        for(int k = 0; k < modes_to_retain; k++)
        {
          HDSA::Ptr<HDSA::Vector<RealT> > p_tmp = p->clone();
          HDSA::Ptr<HDSA::Vector<RealT> > w_tmp = w->clone();
          for(int j = 0; j < num_vecs; j++)
          {
            p_tmp->axpy((*V)(j,k),*(*Ptmp_)[j]);
            w_tmp->axpy((*V)(j,k),*(*Wtmp_)[j]);
          }
          Pr_[block_current_data_step_]->Push_Back(p_tmp);
          Wr_[block_current_data_step_]->Push_Back(w_tmp);
        }

      }

    }

    void Add_Block_Quasi_Newton_Data(void)
    {
      if(use_block_update_)
      {
        Dr_[block_current_data_step_] = Pr_[block_current_data_step_]->MatMat(*Wr_[block_current_data_step_]);
        int m = Pr_[block_current_data_step_]->Number_of_Vectors();
        total_vecs_stored_ += 2*m;
        Pr_[block_current_data_step_+1] = HDSA::makePtr<HDSA::MultiVector<RealT> >();
        Wr_[block_current_data_step_+1] = HDSA::makePtr<HDSA::MultiVector<RealT> >();
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
          if(use_block_update_)
            {
              HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Pr_z = Pr_[block_counter-1]->MatVec(z_in);
              int m = Pr_z->numRows();
              HDSA::Ptr<HDSA::Vector<RealT> > tmp = z_out.clone();
              tmp->set(z_in);
              for(int i = 0; i < m; i++)
                {
                  tmp->axpy(-(*Pr_z)(i,0)/(*Dr_[block_counter-1])(i,i),*(*Wr_[block_counter-1])[i]);
                }
              
              Apply_QN_Inverse_Hessian_Approximation(z_out,*tmp,param_counter,block_counter-1);
              
              HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Wr_z = Wr_[block_counter-1]->MatVec(z_out);
              for(int i = 0; i < m; i++)
                {
                  z_out.axpy(-(*Wr_z)(i,0)/(*Dr_[block_counter-1])(i,i),*(*Pr_[block_counter-1])[i]);
                  z_out.axpy((*Pr_z)(i,0)/(*Dr_[block_counter-1])(i,i),*(*Pr_[block_counter-1])[i]);
                }
            }
          else
            {
              Apply_QN_Inverse_Hessian_Approximation(z_out,z_in,param_counter,block_counter-1);
            }
	      }
    }
    
  
  };

}

#endif
