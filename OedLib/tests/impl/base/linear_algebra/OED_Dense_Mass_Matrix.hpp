#ifndef OEDLIB_TEST_DENSE_MASS_MATRIX_HPP
#define OEDLIB_TEST_DENSE_MASS_MATRIX_HPP

#include "Eigen/Core"

#include <memory>
#include "Mass_Matrix.hpp"

namespace OED_TEST
{
  template <class RealT>
  class Dense_Mass_Matrix
  {
    public:
      using Dense_Matrix = Eigen::Matrix<RealT, Eigen::Dynamic, Eigen::Dynamic>;
      using Dense_Vector = Eigen::Matrix<RealT, Eigen::Dynamic, Eigen::Dynamic>;

    private:
      std::shared_ptr<Dense_Matrix> M_; // dense mass matrix
      Eigen::SelfAdjointEigenSolver<Dense_Matrix> M_eigensolver_; // eigen solver

    public:
      Dense_Mass_Matrix(Dense_Matrix& M) : M_eigensolver_(M)
      {
        this->M_ = std::make_shared<Dense_Matrix>(M);
      }

      void Apply(Vector<RealT> &y_out,
                 Vector<RealT> &x_in)
      {
        auto &y = dynamic_cast<Test_Vector<RealT> &>(y_out);
        auto &x = dynamic_cast<Test_Vector<RealT> &>(x_in);
        // TODO: implement this
        Dense_Vector v = (*this->M_) * x.Vec();
        y.Set_Vec(v);
      }

      void Apply_Inverse(Vector<RealT> &y_out,
                         Vector<RealT> &x_in)
      {
        auto &x = dynamic_cast<Test_Vector<RealT> &>(x_in);
        auto &y = dynamic_cast<Test_Vector<RealT> &>(y_out);
        auto &V = this->M_eigensolver_.eigenvectors();
        auto &w = this->M_eigensolver_.eigenvalues();
        Dense_Vector v = V.transpose() * x.Vec();
        v = v.cwiseProduct(w.cwiseInverse());
        v = V * v;
        y.Set_Vec(v);
      }

      void Apply_Sqrt(Vector<RealT> &x_out,
                      Vector<RealT> &x_in)
      {
      }

      void Apply_Inverse_Sqrt(Vector<RealT> &x_out,
                              Vector<RealT> &x_in)
      {
      }

      Dense_Matrix &M()
      {
        return *this->M_;
      }
  };
  
}

#endif // OEDLIB_TEST_DENSE_MASS_MATRIX_HPP
