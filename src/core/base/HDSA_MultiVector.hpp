#ifndef HDSA_MULTIVECTOR_HPP
#define HDSA_MULTIVECTOR_HPP

#include <sys/stat.h>

namespace HDSA
{

  template <class RealT>
  class MultiVector
  {

  private:
    int num_vecs_;
    std::vector<HDSA::Ptr<HDSA::Vector<RealT>>> vecs_;

  public:
    MultiVector() : num_vecs_(0) {};

    MultiVector(const int num_vecs, const HDSA::Vector<RealT> &vec) : num_vecs_(num_vecs)
    {
      vecs_.resize(num_vecs);
      for (int k = 0; k < num_vecs; k++)
      {
        vecs_[k] = vec.clone();
      }
    }

    MultiVector(std::vector<HDSA::Ptr<HDSA::Vector<RealT>>> &vecs) : vecs_(vecs)
    {
      num_vecs_ = vecs.size();
    }

    virtual ~MultiVector()
    {
    }

    void Clear(void)
    {
      vecs_.clear();
      num_vecs_ = 0;
    }

    // push back vector
    void push_back(const HDSA::Ptr<HDSA::Vector<RealT>> &vec)
    {
      vecs_.push_back(vec);
      ++num_vecs_;
    }

    // Access the kth vector
    HDSA::Ptr<HDSA::Vector<RealT>> operator[](int k) const
    {
      return vecs_[k];
    }

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> MatVec(const HDSA::Vector<RealT> &x) const
    {
      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> Ax = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(num_vecs_, 1);
      for (int k = 0; k < num_vecs_; k++)
      {
        Ax->Replace_Element(k, 0, vecs_[k]->dot(x));
      }
      return Ax;
    }

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> MatMat(const HDSA::MultiVector<RealT> &x) const
    {
      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> C = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(x.Number_of_Vectors(), num_vecs_);
      for (int i = 0; i < x.Number_of_Vectors(); i++)
      {
        for (int j = 0; j < num_vecs_; j++)
        {
          C->Replace_Element(i, j, x[i]->dot(*vecs_[j]));
        }
      }
      return C;
    }

    int Number_of_Vectors(void) const
    {
      return num_vecs_;
    }

    void zeros(void)
    {
      for (int k = 0; k < num_vecs_; k++)
      {
        vecs_[k]->zeros();
      }
    }

    void randomize_standard_normal(void)
    {
      for (int k = 0; k < num_vecs_; k++)
      {
        vecs_[k]->randomize_standard_normal();
      }
    }

    std::vector<RealT> norms(void) const
    {
      std::vector<RealT> n = std::vector<RealT>(num_vecs_);
      for (int k = 0; k < num_vecs_; k++)
      {
        n[k] = vecs_[k]->norm();
      }
      return n;
    }

    void axpy(const RealT &alpha, const HDSA::MultiVector<RealT> &y)
    {
      if (y.Number_of_Vectors() != num_vecs_)
      {
        HDSA_TEST_FOR_EXCEPTION(true, std::logic_error,
                                "Called axpy on HDSA::MultiVector, but x and y do not have the same number of vectors" << std::endl);
      }
      for (int k = 0; k < num_vecs_; k++)
      {
        vecs_[k]->axpy(alpha, *y[k]);
      }
    }

    void axpy(const RealT &alpha, const HDSA::Vector<RealT> &y)
    {
      for (int k = 0; k < num_vecs_; k++)
      {
        vecs_[k]->axpy(alpha, y);
      }
    }

    void scale(const RealT &alpha)
    {
      for (int k = 0; k < num_vecs_; k++)
      {
        vecs_[k]->scale(alpha);
      }
    }

    void Write_to_File(const std::string &name)
    {
      mkdir(name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      for (int k = 0; k < num_vecs_; k++)
      {
        std::string name_k = name + "/Vector_" + std::to_string(k + 1) + ".txt";
        vecs_[k]->Write_to_File(name_k);
      }
    }
  };

}

#endif
