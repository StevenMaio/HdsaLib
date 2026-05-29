//
// Created by Steven Maio on 5/28/26.
//

#ifndef OEDLIB_OED_DENSE_MATRIX_HPP
#define OEDLIB_OED_DENSE_MATRIX_HPP

#include <Eigen/Dense>

using Eigen::MatrixXd;

namespace OED
{
  template <class RealT>
  class Dense_Matrix
  {
  private:
    int num_rows_;
    int num_cols_;
    MatrixXd data_;

  public:
    Dense_Matrix(int num_rows, int num_cols) : num_rows_(num_rows), num_cols_(num_cols), data_(num_rows, num_cols) {}

    Dense_Matrix(MatrixXd &data) : num_rows_(data.rows()), num_cols_(data.cols()), data_(data) {}

    ~Dense_Matrix() {}

    double Compute_Determinant() const
    {
      // TODO: fix this crime
      return this->data_.determinant();
    }

    // TODO: add multiply by inverse
    void Right_Inverse_Multiply(Dense_Matrix &M_out, Dense_Matrix &M_in)
    {
      // TODO: fix this crime
      M_out.data_ = this->data_.inverse() * M_in.data_;
    }

    // TODO: fix this crime as well
    Dense_Matrix Select_Subsquare_Matrix(const std::vector<int> &selection)
    {
      MatrixXd submat = this->data_(selection, selection);
      Dense_Matrix Submatrix(submat);
      return Submatrix;
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

    // Access the (i,j) element
    RealT operator()(int i, int j) const
    {
      return this->data_(i, j);
    }

    // Overwrite the (i,j) element
    void Set_Entry(int i, int j, RealT val)
    {
      this->data_(i, j) = val;
    }

    // TODO: delete this later
    MatrixXd &Data()
    {
      return this->data_;
    }

  };
}
#endif //OEDLIB_OED_DENSE_MATRIX_HPP
