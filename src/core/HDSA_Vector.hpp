#ifndef HDSA_VECTOR_HPP
#define HDSA_VECTOR_HPP

#include <fstream>
#include <random>

namespace HDSA
{

template <class RealT>
class Vector {

protected:
  bool enforce_zeros_;
  std::vector<int> map_full_to_reduced_;
  std::vector<int> map_reduced_to_full_;
  std::vector<int> vec_zeros_;
  std::normal_distribution<RealT> distribution_;

public:

  Vector() 
  {
    enforce_zeros_ = false;
    distribution_ = std::normal_distribution<RealT>(0.0,1.0); 
  }

  virtual ~Vector()
  { }

  // Access the kth element
  virtual RealT operator () (int k) const = 0;

  // Replace the kth element of the vector by val
  virtual void Replace_Element(int k, RealT val) = 0;

  // Get the data on this processor
  virtual std::vector<RealT> Get_Data_on_Processor(void) const = 0;

  // Get the indices on this processor
  virtual std::vector<int> Get_Indices_on_Processor(void) const = 0;

  // Clone the vector
  HDSA::Ptr<HDSA::Vector<RealT> > Clone() const 
  {
    HDSA::Ptr<HDSA::Vector<RealT> > vec = Define_Clone();
    if(enforce_zeros_)
      {
	vec->enforce_zeros_ = true;
	vec->map_full_to_reduced_ = map_full_to_reduced_;
	vec->map_reduced_to_full_ = map_reduced_to_full_;
	vec->vec_zeros_ = vec_zeros_;
      }
    return vec;
  }

  // Clone the vector
  virtual HDSA::Ptr<HDSA::Vector<RealT> > Define_Clone() const = 0;

  // add x to this
  virtual void plus( const HDSA::Vector<RealT> & x ) = 0;

  // scale this by val
  virtual void scale( const RealT val ) = 0;

  // compute the dot product of this and x
  virtual RealT dot( const HDSA::Vector<RealT> &x ) const = 0;

  // compute the norm of this
  virtual RealT norm(void) const = 0;

  // add alpha*x to this
  virtual void axpy( const RealT alpha, const HDSA::Vector<RealT> &x ) = 0;
 
  // set this=0
  virtual void zero(void) = 0;

  // set this= ith canonical basis vector
  virtual void basis( const int i ) = 0;
 
  // return vector dimension
  virtual int dimension() const = 0;

  // set this=x
  virtual void set( const HDSA::Vector<RealT> &x ) = 0;

  // set this=val elementwise
  virtual void setScalar( const RealT val ) = 0;

  // set entries of this to random numbers in [l,u]
  virtual void randomize( const RealT l = 0.0, const RealT u = 1.0 ) = 0;

  virtual HDSA::Ptr<HDSA::Vector<RealT> > Generate_Gaussian_Random_Vector(void)
  {
    // Populate vectors with standard normal samples
    std::default_random_engine generator;
    HDSA::Ptr<HDSA::Vector<RealT> > vec = this->Clone();
    for(int l = 0; l < vec->dimension(); l++)
      {
	vec->Replace_Element(l,distribution_(generator));
      }
    return vec;
  }

  void Write_to_File(std::string & name) const
  {
      std::ofstream fout;
      fout.open(name);
      for(int k = 0; k < this->dimension(); k++)
	{
	  fout << std::setprecision(16) << (*this)(k) << std::endl;
	}
      fout.close();
  }

  // Access the kth element
  RealT Get_Element(int k) const 
  {
    return (*this)(k);
  }

  bool Get_enforce_zeros(void) const
  {
    return enforce_zeros_;
  }

  std::vector<int> Get_map_full_to_reduced(void) const
  {
    return map_full_to_reduced_;
  }

  std::vector<int> Get_map_reduced_to_full(void) const
  {
    return map_reduced_to_full_;
  }

  int Get_map_full_to_reduced(int i) const
  {
    int k = i;
    if(map_full_to_reduced_.size()>0)
      {
	k = map_full_to_reduced_[i];
      }
    return k;
  }

  int Get_map_reduced_to_full(int i) const
  {
    int k = i;
    if(map_reduced_to_full_.size()>0)
      {
	k = map_reduced_to_full_[i];
      }
    return k;
  }

  int Get_nonzero_dim(void) const
  {
    int k = map_reduced_to_full_.size();
    if(k == 0)
      {
	k = this->dimension();
      }
    return k;
  }

  bool Is_entry_zero(int i) const
  {
    bool b = false;
    if(vec_zeros_.size()>0)
      {
	if(vec_zeros_[i]==0)
	  {
	    b = true;
	  }
      }
    return b;
  }

  std::vector<int> Get_vec_zeros(void) const
  {
    return vec_zeros_;
  }

  void Set_Zeros(void)
  {
    for(int k = 0; k < this->dimension(); k++)
      {
	if(Is_entry_zero(k))
	  {
	    Replace_Element(k,0.0);
	  }
      }
  }

  void Enforce_Zeros(void) 
  {
    enforce_zeros_ = true;
    int dim = this->dimension();
    int nonzero_dim = dim;
      if(enforce_zeros_)
	{
	  vec_zeros_ = std::vector<int>(dim,0);
	  map_full_to_reduced_.resize(dim);   
	  nonzero_dim = 0;
	  for(int k = 0; k < dim; k++)
	    {
	      if((*this)(k) != 0.0)
		{
		  vec_zeros_[k] = 1;
		  map_full_to_reduced_[k] = nonzero_dim;
		  nonzero_dim += 1;
		}
	      else
		{
		  vec_zeros_[k] = 0;
		}
	    }
	  
	  map_reduced_to_full_.resize(nonzero_dim);
	  int counter = 0;
	  for(int k = 0; k < dim; k++)
	    {
	      if(vec_zeros_[k] == 1)
		{
		  map_reduced_to_full_[counter] = k;
		  counter += 1;
		}
	    }
	  
	}
  }

  void Check_Vector(void)
  {
    bool pass = true;
    RealT tol = 1.e-14;

    this->zero();
    if(this->norm() != 0.0)
      {
	std::cout << "Failed test 1" << std::endl;
	pass = false;
      }

    this->setScalar(2.0);
    RealT d = static_cast<RealT>(this->dimension());
    if(std::abs(this->norm() - 2.0*std::sqrt(d)) > tol)
      {
	std::cout << "Failed test 2" << std::endl;
	pass = false;
      }

    this->basis(0);
    if(std::abs(this->norm() - 1.0) > tol)
      {
	std::cout << "Failed test 3" << std::endl;
	pass = false;
      }

    this->zero();
    for(int i = 0; i < this->dimension(); i++)
      {
	this->Replace_Element(i,4.0);
      }
    if(std::abs(this->norm() - 4.0*std::sqrt(d)) > tol)
      {
	std::cout << "Failed test 4" << std::endl;
	pass = false;
      }

    this->scale(0.5);
    if(std::abs(this->norm() - 2.0*std::sqrt(d)) > tol)
      {
	std::cout << "Failed test 5" << std::endl;
	pass = false;
      }

    HDSA::Ptr<HDSA::Vector<RealT> > vec = this->Clone();
    vec->setScalar(1.0);
    this->basis(0);
    vec->axpy(2.0,*this);
    if(std::abs(vec->norm() - std::sqrt(9.0 + d-1.0)) > tol)
      {
	std::cout << "Failed test 6" << std::endl;
	pass = false;
      }

    vec->setScalar(3.0);
    this->setScalar(5.0);
    vec->plus(*this);
    if(std::abs(vec->dot(*this) - 40.0*d) > tol)
      {
	std::cout << "Failed test 7" << std::endl;
	pass = false;
      }

    if((*vec)(this->dimension()-1) != 8.0)
      {
	std::cout << "Failed test 8" << std::endl;
	pass = false;
      }

    vec->randomize(3.0,5.0);
    this->set(*vec);
    vec->axpy(-1.0,*this);
    if(vec->norm() != 0.0)
      {
	std::cout << "Failed test 9" << std::endl;
	pass = false;
      }

    if(pass == false)
      {
	std::cout << "Vector check failed" << std::endl;
      }
    else
      {
	std::cout << "Vector check passed" << std::endl;
      }

  }
 
};

}

#endif
