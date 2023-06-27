#ifndef HDSA_BAYES_POSTERIOR_DATA_HPP
#define HDSA_BAYES_POSTERIOR_DATA_HPP

namespace HDSA
{

  template <class RealT>
  class Bayes_Posterior_Data{

  public:
    RealT alpha;
    int N;
    HDSA::Ptr<HDSA::MultiVector<RealT> > Z;
    HDSA::Ptr<HDSA::MultiVector<RealT> > Y;
    HDSA::Ptr<HDSA::MultiVector<RealT> > Gamma_inv_Z;
    HDSA::Ptr<HDSA::Vector<RealT> > Gamma_inv_z_opt;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > G;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > g_vecs;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Lambda;
    HDSA::Ptr<HDSA::Vector<RealT> > state_grad;
    HDSA::Ptr<HDSA::MultiVector<RealT> > u_ell;
    std::vector<HDSA::Ptr<HDSA::MultiVector<RealT> > > u_i_ell;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > a_ell;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > b_i_ell;

    Bayes_Posterior_Data(void) 
    { }

    ~Bayes_Posterior_Data(void)
    { }
  
    void Compute_Posterior_Data(const HDSA::Model_Discrepancy_Interface<RealT> & md_interface, RealT alpha_in, const HDSA::Vector<RealT> & u_opt, const HDSA::Vector<RealT> & z_opt)
    {
      alpha = alpha_in;
      Z = md_interface.Load_Z_Data();
      Y = md_interface.Load_Y_Data();
      N = Z->Number_of_Vectors();
      state_grad = u_opt.clone();
      md_interface.Misfit_Gradient(*state_grad,u_opt,z_opt);
      
      Gamma_inv_Z = HDSA::makePtr<HDSA::MultiVector<RealT> >(N,z_opt);
      for(int k = 0; k < N; k++)
	{
	  HDSA::Ptr<HDSA::Vector<RealT> > zk = (*Z)[k];
	  HDSA::Ptr<HDSA::Vector<RealT> > gzk = (*Gamma_inv_Z)[k];
	  md_interface.Apply_Gamma_Mat_Inverse(*gzk,*zk);
	}

      G = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(N,N);
      Gamma_inv_z_opt = z_opt.clone(); 
      md_interface.Apply_Gamma_Mat_Inverse(*Gamma_inv_z_opt,z_opt);
      RealT z_opt_Gamma_inv_z_opt = z_opt.dot(*Gamma_inv_z_opt);
      for(int i = 0; i < N; i++)
	{
	  HDSA::Ptr<HDSA::Vector<RealT> > zi = (*Z)[i];
	  HDSA::Ptr<HDSA::Vector<RealT> > gzi = (*Gamma_inv_Z)[i];
	  RealT vali = 1.0 + z_opt_Gamma_inv_z_opt - zi->dot(*Gamma_inv_z_opt);
	  for(int j = 0; j < i+1; j++)
	    {
	      HDSA::Ptr<HDSA::Vector<RealT> > zj = (*Z)[j];
	      RealT val = vali;
	      val -= zj->dot(*Gamma_inv_z_opt);
	      val += zj->dot(*gzi);
	      G->Replace_Element(i,j,val);
	    }
	}
      for(int i = 0; i < N; i++)
	{
	  for(int j = i+1; j < N; j++)
	    {
	      G->Replace_Element(i,j,(*G)(j,i));
	    }
	}

      g_vecs = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(N,N);
      Lambda = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(N,1);
      HDSA::Linear_Algebra::Symmetric_Eig_Decomposition<RealT>(*G, *g_vecs, *Lambda);

      u_ell = HDSA::makePtr<HDSA::MultiVector<RealT> >(N,*(*Y)[0]);      
      for(int ell = 0; ell < N; ell++)
	{
	  HDSA::Ptr<HDSA::Vector<RealT> > yl = (*Y)[ell];
	  HDSA::Ptr<HDSA::Vector<RealT> > ul = (*u_ell)[ell];
	  md_interface.Apply_L_Mat_Inverse(*ul,*yl);
	}

      u_i_ell.resize(N);
      for(int i = 0; i < N; i++)
	{
	  u_i_ell[i] = HDSA::makePtr<HDSA::MultiVector<RealT> >(N,*(*Y)[0]);
	  for(int ell = 0; ell < N; ell++)
	    {
	      HDSA::Ptr<HDSA::Vector<RealT> > uil = (*u_i_ell[i])[ell];
	      HDSA::Ptr<HDSA::Vector<RealT> > ul = (*u_ell)[ell];
	      md_interface.Apply_L_Plus_beta_Identity_Mat_Inverse(*uil,*ul,(*Lambda)(i,0)/alpha);
	      uil->scale(1.0/alpha);
	    }
	}

      a_ell = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(N,1);
      b_i_ell = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(N,N);
      for(int ell = 0; ell < N; ell++)
	{
	  HDSA::Ptr<HDSA::Vector<RealT> > zl = (*Z)[ell];
	  RealT val_a = 1.0 - zl->dot(*Gamma_inv_z_opt) + z_opt_Gamma_inv_z_opt;
	  a_ell->Replace_Element(ell,0,val_a);
	  for(int i = 0; i < N; i++)
	    {
	      RealT val_b = 0.0;
	      for(int k = 0; k < N; k++)
		{
		  HDSA::Ptr<HDSA::Vector<RealT> > gzk = (*Gamma_inv_Z)[k];
		  val_b += (*g_vecs)(k,i)*(zl->dot(*gzk) - gzk->dot(z_opt) + (*a_ell)(ell,0));
		}
	      b_i_ell->Replace_Element(i,ell,val_b);
	    }
	}
      
    }

  };

}

#endif
