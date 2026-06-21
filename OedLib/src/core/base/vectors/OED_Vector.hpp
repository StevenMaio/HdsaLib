/***********************************************************************
 HdsaLib - A library for Hyper-differential Sensitivity Analysis
 
 Questions? Contact Joseph Hart (joshart@sandia.gov)
************************************************************************/

#ifndef HDSA_VECTOR_HPP
#define HDSA_VECTOR_HPP

#include <cmath>
#include <fstream>
#include <iomanip>

namespace OED
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

    // // Clone the vector
    // virtual OED::Ptr<OED::Vector<RealT>> Clone() const = 0;

    // compute the Dot product of this and x
    virtual RealT Dot(const OED::Vector<RealT> &x) const = 0;

    // add alpha*x to this
    virtual void Scaled_Plus(const RealT alpha, const OED::Vector<RealT> &x) = 0;

    // return vector Dimension
    virtual int Dimension() const = 0;

    // Set this=val elementwise
    virtual void Set_Scalar(const RealT val) = 0;

    // virtual void Randomize_Standard_Normal() = 0;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Virtual functions available for convienence when useful, these are not called within HdsaLib but rather are for the user to call from main rather than going through casts
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // virtual void Write_to_File(const std::string &name) const
    // {
    //   HDSA_TEST_FOR_EXCEPTION(true, std::logic_error,
    //                           "Error in OED::Vector: Write_to_File has not been implemented for this vector type" << std::endl);
    // }
    //
    virtual RealT Get_Entry(int k) const
    {
      // HDSA_TEST_FOR_EXCEPTION(true, std::logic_error,
      //                         "Error in OED::Vector: Get_Entry has not been implemented for this vector type" << std::endl);
      RealT val = 0.0;
      return val;
    }

    virtual void Set_Entry(int k, RealT val)
    {
      // HDSA_TEST_FOR_EXCEPTION(true, std::logic_error,
      //                         "Error in OED::Vector: Set_Entry has not been implemented for this vector type" << std::endl);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Implementations using the pure virtual functions above, may be overloaded if an efficiency gain is possible
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Scale this by val
    virtual void Scale(const RealT val)
    {
      this->Scaled_Plus(val - 1.0, *this);
    }

    // Set this=0
    virtual void Zeros(void)
    {
      this->Set_Scalar(0.0);
    }

    // compute the Norm of this
    virtual RealT Norm(void) const
    {
      return std::sqrt(this->Dot(*this));
    }

    // add x to this
    virtual void Plus(const OED::Vector<RealT> &x)
    {
      this->Scaled_Plus(1.0, x);
    }

    // Set this=x
    virtual void Set(const OED::Vector<RealT> &x)
    {
      this->Scale(0.0);
      this->Scaled_Plus(1.0, x);
    }

    virtual void Write_To_File(const std::string& file_name)
    {
      std::ofstream fout;
      fout.open(file_name);
      for (int i = 0; i < this->Dimension(); i++)
      {
        fout << std::setprecision(16) << this->Get_Entry(i) << " ";
      }
      fout.close();
    }

    // Test vector implementation
    // void Test_Vector(void)
    // {
    //   bool pass = true;
    //   RealT tol = 1.e-14;
    //
    //   this->Zeros();
    //   if (this->Norm() != 0.0)
    //   {
    //     std::cout << "Failed test 1" << std::endl;
    //     pass = false;
    //   }
    //
    //   this->Set_Scalar(2.0);
    //   RealT d = static_cast<RealT>(this->Dimension());
    //   if (std::abs(this->Norm() - 2.0 * std::sqrt(d)) > tol)
    //   {
    //     std::cout << "Failed test 2" << std::endl;
    //     pass = false;
    //   }
    //
    //   this->Scale(0.5);
    //   if (std::abs(this->Norm() - std::sqrt(d)) > tol)
    //   {
    //     std::cout << "Failed test 3" << std::endl;
    //     pass = false;
    //   }
    //
    //   OED::Ptr<OED::Vector<RealT>> vec = this->Clone();
    //   vec->Set_Scalar(3.0);
    //   this->Set_Scalar(5.0);
    //   vec->Plus(*this);
    //   if (std::abs(vec->Dot(*this) - 40.0 * d) > tol)
    //   {
    //     std::cout << "Failed test 4" << std::endl;
    //     pass = false;
    //   }
    //
    //   if (pass == false)
    //   {
    //     std::cout << "Vector check failed" << std::endl;
    //   }
    //   else
    //   {
    //     std::cout << "Vector check passed" << std::endl;
    //   }
    // }
  };

}

#endif
