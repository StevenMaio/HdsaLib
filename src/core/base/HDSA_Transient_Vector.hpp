#ifndef HDSA_TRANSIENT_VECTOR_HPP
#define HDSA_TRANSIENT_VECTOR_HPP

template <class RealT>
class Transient_Vector : public HDSA::Vector<RealT>
{

private:
    int n_t_;
    std::vector<HDSA::Ptr<HDSA::Vector<RealT>>> vec_;

public:
    Transient_Vector(int n_t, const HDSA::Ptr<HDSA::Vector<RealT>> &spatial_vec)
    {
        n_t_ = n_t;
        vec_.resize(n_t);
        for (int k = 0; k < n_t; k++)
        {
            vec_[k] = spatial_vec->clone();
        }
    }

    ~Transient_Vector()
    {
    }

    //////////////////////////////////////////////////////////////////////////////////
    // Overloading pure virtual functions in HDSA::Vector base class
    //////////////////////////////////////////////////////////////////////////////////

    HDSA::Ptr<HDSA::Vector<RealT>> clone() const override
    {
        HDSA::Ptr<HDSA::Vector<RealT>> vec = HDSA::makePtr<Transient_Vector<RealT>>(n_t_, vec_[0]);
        return vec;
    }

    // compute the dot product of this and x
    RealT dot(const HDSA::Vector<RealT> &x) const override
    {
        RealT val = 0.0;
        const Transient_Vector<RealT> x_trans = dynamic_cast<const Transient_Vector<RealT> &>(x);
        for (int k = 0; k < n_t_; k++)
        {
            val += vec_[k]->dot(*x_trans[k]);
        }
        return val;
    }

    // add alpha*x to this
    void axpy(const RealT alpha, const HDSA::Vector<RealT> &x) override
    {
        const Transient_Vector<RealT> x_trans = dynamic_cast<const Transient_Vector<RealT> &>(x);
        for (int k = 0; k < n_t_; k++)
        {
            vec_[k]->axpy(alpha, *x_trans[k]);
        }
    }

    // return vector dimension
    int dimension() const override
    {
        return n_t_ * vec_[0]->dimension();
    }

    // set this=val elementwise
    void setScalar(const RealT val) override
    {
        for (int k = 0; k < n_t_; k++)
        {
            vec_[k]->setScalar(val);
        }
    }

    void randomize_standard_normal() override
    {
        for (int k = 0; k < n_t_; k++)
        {
            vec_[k]->randomize_standard_normal();
        }
    }

    void Write_to_File(const std::string &name) const override
    {
        int num_char = name.size();
        std::string name_tmp = name.substr(0, num_char-4);
        for(int k = 0; k < n_t_; k++)
        {
            std::string name_k = name_tmp + "_time_" + std::to_string(k+1) + ".txt";
            vec_[k]->Write_to_File(name_k);
        }
    }

    //////////////////////////////////////////////////////////////////////////////////
    // Function specific to this class for convenience
    //////////////////////////////////////////////////////////////////////////////////

    // Access the (i,j) element
    HDSA::Ptr<HDSA::Vector<RealT>> operator[](int k) const
    {
        return vec_[k];
    }
};

#endif
