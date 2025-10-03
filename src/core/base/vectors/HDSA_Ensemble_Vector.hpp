#ifndef HDSA_ENSEMBLE_VECTOR_HPP
#define HDSA_ENSEMBLE_VECTOR_HPP

namespace HDSA
{

  template <class RealT>
  class Ensemble_Vector : public HDSA::Vector<RealT>
  {

  private:
    int num_vecs_;
    std::vector<HDSA::Ptr<HDSA::Vector<RealT>>> vecs_;

  public:
    Ensemble_Vector(const int num_vecs, const HDSA::Vector<RealT> &vec) : num_vecs_(num_vecs)
    {
      vecs_.resize(num_vecs);
      for (int k = 0; k < num_vecs; k++)
      {
        vecs_[k] = vec.Clone();
      }
    }

    Ensemble_Vector(std::vector<HDSA::Ptr<HDSA::Vector<RealT>>> &vecs) : vecs_(vecs)
    {
      num_vecs_ = vecs.size();
    }

    virtual ~Ensemble_Vector()
    {
    }

    // Access the kth vector
    HDSA::Ptr<HDSA::Vector<RealT>> operator[](int k) const
    {
      return vecs_[k];
    }

    int Number_of_Vectors(void) const
    {
      return num_vecs_;
    }

    HDSA::Ptr<HDSA::Vector<RealT>> Clone() const
    {
      HDSA::Ptr<HDSA::Vector<RealT>> ens_vec = HDSA::makePtr<HDSA::Ensemble_Vector<RealT>>(num_vecs_, *vecs_[0]);
      return ens_vec;
    }

    RealT Dot(const HDSA::Vector<RealT> &x) const
    {
      RealT val = 0.0;
      const HDSA::Ensemble_Vector<RealT> x_ens = dynamic_cast<const HDSA::Ensemble_Vector<RealT> &>(x);
      for (int k = 0; k < num_vecs_; k++)
      {
        val += x_ens[k]->Dot(*vecs_[k]);
      }
      return val;
    }

    void Scaled_Plus(const RealT alpha, const HDSA::Vector<RealT> &x)
    {
      const HDSA::Ensemble_Vector<RealT> x_ens = dynamic_cast<const HDSA::Ensemble_Vector<RealT> &>(x);
      for (int k = 0; k < num_vecs_; k++)
      {
        vecs_[k]->Scaled_Plus(alpha, *x_ens[k]);
      }
    }

    int Dimension() const
    {
      int dim = num_vecs_ * vecs_[0]->Dimension();
      return dim;
    }

    void Set_Scalar(const RealT val)
    {
      for (int k = 0; k < num_vecs_; k++)
      {
        vecs_[k]->Set_Scalar(val);
      }
    }

    void Randomize_Standard_Normal()
    {
      for (int k = 0; k < num_vecs_; k++)
      {
        vecs_[k]->Randomize_Standard_Normal();
      }
    }

    void Write_to_File(const std::string &name) const
    {
      for (int k = 0; k < num_vecs_; k++)
      {
        std::string name_k = name.substr(0, name.length() - 4) + "_ens_" + std::to_string(k + 1) + ".txt";
        vecs_[k]->Write_to_File(name_k);
      }
    }

    HDSA::Ptr<HDSA::Vector<RealT>> Ensemble_Average(void) const
    {
      HDSA::Ptr<HDSA::Vector<RealT>> avg_vec = vecs_[0]->Clone();
      avg_vec->Zeros();
      for (int k = 0; k < num_vecs_; k++)
      {
        avg_vec->Plus(*vecs_[k]);
      }
      avg_vec->Scale(1.0 / static_cast<RealT>(num_vecs_));
      return avg_vec;
    }
  };

}

#endif
