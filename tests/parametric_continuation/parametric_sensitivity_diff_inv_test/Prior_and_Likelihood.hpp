#ifndef PRIOR_AND_LIKELIHOOD_HPP
#define PRIOR_AND_LIKELIHOOD_HPP

template <class RealT>
class Prior_and_Likelihood {

public:
  HDSA::Ptr<Adv_Diff_Constraint<RealT> > con_;
  int m_;
  std::vector<int> obs_locations_;
  RealT noise_var_;
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > L_;
  HDSA::Ptr<HDSA::Random_Number_Generator<RealT> > rand_num_gen_;
  std::vector<RealT> data_;

  Prior_and_Likelihood(HDSA::Ptr<Adv_Diff_Constraint<RealT> > & con): con_(con)
  {
    m_ = con->m_;
    obs_locations_ = {0, 14, 28, 42, 57, 71, 85, 99};
    noise_var_ = 9.e-6;
    L_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);
    for(int i = 0; i < m_; i++)
      {
	for(int j = 0; j < m_; j++)
	  {
	    RealT val = (1.0/0.75) * ( (1.e-3)*(*con_->S_)(i,j) + (*con_->M_)(i,j) );
	    L_->Replace_Element(i,j,val);
	  }
      }

    int num_random_numbers = 1.e5;
    std::string random_number_file = "random_numbers.txt";
    rand_num_gen_ = HDSA::makePtr<HDSA::Random_Number_Generator<RealT> >(num_random_numbers,random_number_file);

    data_ = std::vector<RealT>(8);
    RealT val = 0.0;
    // read in data
    std::ifstream in("obs_data.txt");           
    // read the elements in the file into a vector  
    // test file open   
    if (in) {   
      for(int i = 0; i < 8; i++)
	{
	  in >> val;
	  data_[i] = val;
	}   
    }
    else
      {
	std::cout << "Error loading the data from obs_data.txt" << std::endl;
      }  
    
  }

  virtual ~Prior_and_Likelihood()
  { }

  // Objective function J(u,z) = (1/2)*(O*u-d)^T*Sigma^{-1}*(O*u-d) + (1/2)*z^T*Gamma^{-1}*z
  // Where O is the observation operator, and Sigma is the noise covariance, and Gamma is the prior covariance
  
  void Apply_Observation_Operator(std::vector<RealT> & d_out, const HDSA::Vector<RealT> & u_in) const
  {
    const Std_Vector<RealT>& u_in_std = dynamic_cast<const Std_Vector<RealT>&>(u_in);
    for(int k = 0; k < 8; k++)
      {
	d_out[k] = u_in_std(obs_locations_[k]);
      }
  }

  void Apply_Observation_Operator_Transpose(HDSA::Vector<RealT> & u_out, const std::vector<RealT> & d_in) const
  {
    Std_Vector<RealT>& u_out_std = dynamic_cast<Std_Vector<RealT>&>(u_out);
    for(int k = 0; k < 8; k++)
      {
	u_out_std.Replace_Element(obs_locations_[k],d_in[k]);
      }
  }

  void Apply_Noise_Precision(std::vector<RealT> & d_out, const std::vector<RealT> & d_in) const
  {
    for(int k = 0; k < 8; k++)
      {
	d_out[k] = (1.0/noise_var_)*d_in[k];
      }
  }
  
  void Apply_Prior_Precision(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const
  {
    Std_Vector<RealT>& z_out_std = dynamic_cast<Std_Vector<RealT>&>(z_out);
    const Std_Vector<RealT>& z_in_std = dynamic_cast<const Std_Vector<RealT>&>(z_in);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > z_tmp1 = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    for(int k = 0; k < m_; k++)
      {
	z_tmp1->Replace_Element(k,0,z_in_std(k));
      }

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > z_tmp2 = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    L_->Multiply(*z_tmp2,*z_tmp1);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > z_tmp3 = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*con_->M_,*z_tmp3,*z_tmp2);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > z_tmp4 = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    L_->Multiply(*z_tmp4,*z_tmp3);
    
    for(int k = 0; k < m_; k++)
      {
	z_out_std.Replace_Element(k,(*z_tmp4)(k,0));
      }
  }

  void Apply_Prior_Covariance(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const
  {
    Std_Vector<RealT>& z_out_std = dynamic_cast<Std_Vector<RealT>&>(z_out);
    const Std_Vector<RealT>& z_in_std = dynamic_cast<const Std_Vector<RealT>&>(z_in);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > z_tmp1 = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    for(int k = 0; k < m_; k++)
      {
	z_tmp1->Replace_Element(k,0,z_in_std(k));
      }

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > z_tmp2 = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*L_,*z_tmp2,*z_tmp1);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > z_tmp3 = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    con_->M_->Multiply(*z_tmp3,*z_tmp2);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > z_tmp4 = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*L_,*z_tmp4,*z_tmp3);
    
    for(int k = 0; k < m_; k++)
      {
	z_out_std.Replace_Element(k,(*z_tmp4)(k,0));
      }
  }

  void Generate_Prior_Samples(HDSA::MultiVector<RealT> & samples) const
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > R = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);
    HDSA::Linear_Algebra::Cholesky_Factorization(*con_->M_,*R);

    int num_samples = samples.Number_of_Vectors();
    for(int i = 0; i < num_samples; i++)
      {

	HDSA::Ptr<HDSA::Dense_Matrix<RealT> > b = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
        for(int k = 0; k < m_; k++)
          {
            b->Replace_Element(k,0,rand_num_gen_->Generate_Standard_Normal_Sample());
          }
	HDSA::Ptr<HDSA::Dense_Matrix<RealT> > x = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
	
	// Should be E_z^{-1}*R^T*b
	HDSA::Ptr<HDSA::Dense_Matrix<RealT> > tmp = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
	R->Multiply(*tmp, *b, true);
	HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*L_,*x,*tmp);

	Std_Vector<RealT>& vec_out_std = dynamic_cast<Std_Vector<RealT>&>(*samples[i]);
        for(int k = 0; k < m_; k++)
          {
            vec_out_std.Replace_Element(k,(*x)(k,0));
          }
      }
  }
  
};

#endif
