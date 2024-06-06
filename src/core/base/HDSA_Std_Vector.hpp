#ifndef HDSA_STDVECTOR_HPP
#define HDSA_STDVECTOR_HPP

template <class RealT>
class Std_Vector : public HDSA::Vector<RealT> {
  
private:
  int dim_;
  const HDSA::Ptr<HDSA::Random_Number_Generator<RealT> > random_number_generator_;
  HDSA::Ptr<std::vector<RealT> > vec_;

public:  
  Std_Vector(int dim):
    dim_(dim), random_number_generator_(HDSA::makePtr<HDSA::Random_Number_Generator<RealT> >())  
  {
    vec_ = HDSA::makePtr<std::vector<RealT> >(dim,0.0);
  }

  Std_Vector(int dim, int seed):
    dim_(dim), random_number_generator_(HDSA::Random_Number_Generator<RealT>(seed))
  {
    vec_ = HDSA::makePtr<std::vector<RealT> >(dim,0.0);
  }

  Std_Vector(int dim, const HDSA::Ptr<HDSA::Random_Number_Generator<RealT> > & random_number_generator)
    : dim_(dim), random_number_generator_(random_number_generator)
  {
    vec_ = HDSA::makePtr<std::vector<RealT> >(dim,0.0);
  }
  
  ~Std_Vector()
  { }

  //////////////////////////////////////////////////////////////////////////////////
  // Overloading pure virtual functions in HDSA::Vector base class
  //////////////////////////////////////////////////////////////////////////////////
  
  HDSA::Ptr<HDSA::Vector<RealT> > clone() const
  {
    HDSA::Ptr<HDSA::Vector<RealT> > vec = HDSA::makePtr<Std_Vector<RealT> >(dim_,random_number_generator_);
    return vec;
  }

  // compute the dot product of this and x
  RealT dot( const HDSA::Vector<RealT> &x ) const
  {
    RealT val = 0.0;
    const Std_Vector<RealT> x_std = dynamic_cast<const Std_Vector<RealT>&>(x);
    for(int k = 0; k < dim_; k++)
      {
	val += (*vec_)[k]*(x_std(k));
      }
    return val;
  }

  // add alpha*x to this
  void axpy( const RealT alpha, const HDSA::Vector<RealT> &x )
  {
    const Std_Vector<RealT> x_std = dynamic_cast<const Std_Vector<RealT>&>(x);
    for(int k = 0; k < dim_; k++)
      {
	(*vec_)[k] += alpha*(x_std(k));
      }
  }

  // return vector dimension
  int dimension() const
  {
    return dim_;
  }

  // set this=val elementwise
  void setScalar( const RealT val )
  {
    for(int k = 0; k < dim_; k++)
      {
	(*vec_)[k] = val;
      }
  }

  void randomize_standard_normal( ) 
  {
    for(int k = 0; k < dim_; k++)
      {
	(*vec_)[k] = random_number_generator_->Generate_Standard_Normal_Sample();
      }
  }

  void Write_to_File(std::string & name) const
  {
    std::ofstream fout;
    fout.open(name);
    for(int i = 0; i < dim_; i++)
      {
        fout << std::setprecision(16) << (*vec_)[i] << "  ";
      }
    fout.close();
  }

  //////////////////////////////////////////////////////////////////////////////////
  // Function specific to this class for convenience
  //////////////////////////////////////////////////////////////////////////////////

  // Access underlying std::vector
  const HDSA::Ptr<std::vector<RealT> > get_std_vec(void) const
  {
    return vec_;
  }
  
  // Access the (i,j) element
  RealT operator () (int k) const 
  {
    return (*vec_)[k];
  }

  // Replace the kth element of the vector by val
  void Replace_Element(int k, RealT val)
  {
    (*vec_)[k] = val;
  } 

};


#endif
