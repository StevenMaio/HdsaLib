//
// Created by Steven Maio on 5/28/26.
//

#ifndef OEDLIB_OED_DENSE_MATRIX_HPP
#define OEDLIB_OED_DENSE_MATRIX_HPP

#include <memory>
#include "Eigen/Dense"

#include "OED_Ptr.hpp"

using Eigen::MatrixXd;

namespace OED
{
  template <class RealT>
  class Dense_Matrix
  {
  private:
    int num_rows_;
    int num_cols_;
    Ptr<MatrixXd> data_;

  public:
    Dense_Matrix(int num_rows, int num_cols) : num_rows_(num_rows), num_cols_(num_cols)
    {
      this->data_ = OED::makePtr<MatrixXd>(num_rows, num_cols);
    }

    Dense_Matrix(MatrixXd &data) : num_rows_(data.rows()), num_cols_(data.cols())
    {
      this->data_ = OED::makePtr<MatrixXd>(data);
    }

    ~Dense_Matrix() {}

    double Compute_Determinant() const
    {
      // TODO: is there a better way to implement this?
      return this->data_->determinant();
    }

    void Right_Inverse_Multiply(Dense_Matrix &M_out, Dense_Matrix &M_in)
    {
      // TODO: is there a better way to implement this?
      *M_out.data_ = this->data_->inverse() * (*M_in.data_);
    }

    // TODO: fix this crime as well
    Ptr<Dense_Matrix> Select_Subsquare_Matrix(const std::vector<int> &selection)
    {
      MatrixXd submat = (*this->data_)(selection, selection);
      Ptr<Dense_Matrix> mat = OED::makePtr<Dense_Matrix>(submat);
      return mat;
    }

    // TODO: add multiply function

    int Number_of_Rows() const
    {
      return this->num_rows_;
    }

    int Number_of_Cols() const
    {
      return this->num_cols_;
    }

    void Scaled_Plus(const RealT alpha, const Dense_Matrix<RealT> &A)
    {
      *this->data_ += *A->data_;
    }

    // Access the (i,j) element
    RealT operator()(int i, int j) const
    {
      return (*this->data_)(i, j);
    }

    // Overwrite the (i,j) element
    void Set_Entry(int i, int j, RealT val)
    {
      (*this->data_)(i, j) = val;
    }

    // TODO: delete this later
    MatrixXd &Data()
    {
      return *this->data_;
    }

  };
}
#endif //OEDLIB_OED_DENSE_MATRIX_HPP
