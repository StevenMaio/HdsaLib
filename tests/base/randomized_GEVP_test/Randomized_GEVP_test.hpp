#ifndef RANDOMIZED_GEVP_TEST_HPP
#define RANDOMIZED_GEVP_TEST_HPP


template <class RealT>
class Randomized_GEVP_test : public HDSA::Randomized_GEVP<RealT> {

private:
  int m_;
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > S_;
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > M_;
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > A_;
 const HDSA::Ptr<HDSA::Random_Number_Generator<RealT> > random_number_generator_;
  
public:

  Randomized_GEVP_test(HDSA::Vector<RealT> & vec, const HDSA::Ptr<HDSA::Random_Number_Generator<RealT> > & random_number_generator)
    : HDSA::Randomized_GEVP<RealT>(vec), random_number_generator_(random_number_generator)
  {
    m_ = vec.dimension();
    RealT h = 1.0/static_cast<RealT>(m_-1);

    S_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);
    M_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);
    A_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);

    S_->Replace_Element(0,0,1.0/h);
    S_->Replace_Element(0,1,-1.0/h);
    for(int i = 1; i < m_-1; i++)
      {
        S_->Replace_Element(i,i,2.0/h);
        S_->Replace_Element(i,i-1,-1.0/h);
        S_->Replace_Element(i,i+1,-1.0/h);
      }
    S_->Replace_Element(m_-1,m_-2,-1.0/h);
    S_->Replace_Element(m_-1,m_-1,1.0/h);

    M_->Replace_Element(0,0,(1.0/3.0)*h);
    M_->Replace_Element(0,1,(1.0/6.0)*h);
    for(int i = 1; i < m_-1; i++)
      {
	M_->Replace_Element(i,i,(2.0/3.0)*h);
        M_->Replace_Element(i,i-1,(1.0/6.0)*h);
	M_->Replace_Element(i,i+1,(1.0/6.0)*h);
      }
    M_->Replace_Element(m_-1,m_-2,(1.0/6.0)*h);
    M_->Replace_Element(m_-1,m_-1,(1.0/3.0)*h);

    for(int i = 0; i < m_; i++)
      {
	for(int j = 0; j < m_; j++)
	  {
	    RealT val = (1.e-2)*(*S_)(i,j) + (*M_)(i,j);
	    A_->Replace_Element(i,j,val);
	  }
      }

  }
  
  void Apply_Operator(HDSA::Vector<RealT> & vec_out, const HDSA::Vector<RealT> & vec_in) const 
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > b = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    const Std_Vector<RealT>& vec_in_std = dynamic_cast<const Std_Vector<RealT>&>(vec_in);
    Std_Vector<RealT>& vec_out_std = dynamic_cast<Std_Vector<RealT>&>(vec_out);
    for(int k = 0; k < m_; k++)
      {
        b->Replace_Element(k,0,vec_in_std(k));
      }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > x = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*A_,*x,*b);
    for(int k = 0; k < m_; k++)
      {
	vec_out_std.Replace_Element(k,(*x)(k,0));
      }
  }
 
  void Apply_Weighting_Operator(HDSA::Vector<RealT> & vec_out, const HDSA::Vector<RealT> & vec_in) const 
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > b = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    const Std_Vector<RealT>& vec_in_std = dynamic_cast<const Std_Vector<RealT>&>(vec_in);
    Std_Vector<RealT>& vec_out_std = dynamic_cast<Std_Vector<RealT>&>(vec_out);
    for(int k = 0; k < m_; k++)
      {
        b->Replace_Element(k,0,vec_in_std(k));
      }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > x = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    M_->Multiply(*x,*b);
    for(int k = 0; k < m_; k++)
      {
	vec_out_std.Replace_Element(k,(*x)(k,0));
      }
  }
  
  void Apply_Weighting_Operator_Inverse(HDSA::Vector<RealT> & vec_out, const HDSA::Vector<RealT> & vec_in) const 
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > b = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    const Std_Vector<RealT>& vec_in_std = dynamic_cast<const Std_Vector<RealT>&>(vec_in);
    Std_Vector<RealT>& vec_out_std = dynamic_cast<Std_Vector<RealT>&>(vec_out);
    for(int k = 0; k < m_; k++)
      {
        b->Replace_Element(k,0,vec_in_std(k));
      }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > x = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*M_,*x,*b);
    for(int k = 0; k < m_; k++)
      {
        vec_out_std.Replace_Element(k,(*x)(k,0));
      }
  }

  void Generate_Random_Samples(HDSA::MultiVector<RealT> & samples) const 
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > R = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,m_);
    HDSA::Linear_Algebra::Cholesky_Factorization(*M_,*R);

    int num_samples = samples.Number_of_Vectors();
    for(int i = 0; i < num_samples; i++)
      {

	HDSA::Ptr<HDSA::Dense_Matrix<RealT> > b = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
	HDSA::Ptr<Std_Vector<RealT> > vec_in_std = HDSA::makePtr<Std_Vector<RealT> >(m_,random_number_generator_);
	vec_in_std->randomize_standard_normal();
	for(int k = 0; k < m_; k++)
	  {
	    b->Replace_Element(k,0,(*vec_in_std)(k));
	  }
	HDSA::Ptr<HDSA::Dense_Matrix<RealT> > x = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
	HDSA::Linear_Algebra::Upper_Tri_Solve<RealT>(*x,*b,*R);
	Std_Vector<RealT>& vec_out_std = dynamic_cast<Std_Vector<RealT>&>(*samples[i]);
	for(int k = 0; k < m_; k++)
	  {
	    vec_out_std.Replace_Element(k,(*x)(k,0));
	  }
      }

  }

};

#endif
