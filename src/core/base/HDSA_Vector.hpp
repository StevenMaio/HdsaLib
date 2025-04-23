#ifndef HDSA_VECTOR_HPP
#define HDSA_VECTOR_HPP

namespace HDSA
{

  template <class RealT>
  class Vector
  {

  public:
    Vector()
    {
    }

    virtual ~Vector()
    {
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Pure virtual functions to define when creating a vector interface
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Clone the vector
    virtual HDSA::Ptr<HDSA::Vector<RealT>> clone() const = 0;

    // compute the dot product of this and x
    virtual RealT dot(const HDSA::Vector<RealT> &x) const = 0;

    // add alpha*x to this
    virtual void axpy(const RealT alpha, const HDSA::Vector<RealT> &x) = 0;

    // return vector dimension
    virtual int dimension() const = 0;

    // set this=val elementwise
    virtual void setScalar(const RealT val) = 0;

    virtual void randomize_standard_normal() = 0;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Virtual functions available for convienence when useful, these are not called within HdsaLib but rather are for the user to call from main rather than going through casts
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    virtual void Write_to_File(const std::string &name) const
    {
      std::ofstream fout;
      fout.open(name);
      fout << "Write_to_File has not been implemented for this vector type";
      fout.close();
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Implementations using the pure virtual functions above, may be overloaded if an efficiency gain is possible
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // scale this by val
    virtual void scale(const RealT val)
    {
      this->axpy(val - 1.0, *this);
    }

    // set this=0
    virtual void zeros(void)
    {
      this->setScalar(0.0);
    }

    // compute the norm of this
    virtual RealT norm(void) const
    {
      return std::sqrt(this->dot(*this));
    }

    // add x to this
    virtual void plus(const HDSA::Vector<RealT> &x)
    {
      this->axpy(1.0, x);
    }

    // set this=x
    virtual void set(const HDSA::Vector<RealT> &x)
    {
      this->scale(0.0);
      this->axpy(1.0, x);
    }

    // Test vector implementation
    void Test_Vector(void)
    {
      bool pass = true;
      RealT tol = 1.e-14;

      this->zeros();
      if (this->norm() != 0.0)
      {
        std::cout << "Failed test 1" << std::endl;
        pass = false;
      }

      this->setScalar(2.0);
      RealT d = static_cast<RealT>(this->dimension());
      if (std::abs(this->norm() - 2.0 * std::sqrt(d)) > tol)
      {
        std::cout << "Failed test 2" << std::endl;
        pass = false;
      }

      this->scale(0.5);
      if (std::abs(this->norm() - std::sqrt(d)) > tol)
      {
        std::cout << "Failed test 3" << std::endl;
        pass = false;
      }

      HDSA::Ptr<HDSA::Vector<RealT>> vec = this->clone();
      vec->setScalar(3.0);
      this->setScalar(5.0);
      vec->plus(*this);
      if (std::abs(vec->dot(*this) - 40.0 * d) > tol)
      {
        std::cout << "Failed test 4" << std::endl;
        pass = false;
      }

      if (pass == false)
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
