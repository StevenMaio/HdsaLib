#ifndef HDSA_TRANSIENT_VECTOR_CONST_HPP
#define HDSA_TRANSIENT_VECTOR_CONST_HPP

namespace HDSA
{

    template <class RealT>
    class Transient_Vector_Const : public Transient_Vector<RealT>
    {

    private:
        int n_t_;
        std::vector<HDSA::Ptr<const HDSA::Vector<RealT>>> const_vec_;

    public:
        Transient_Vector_Const(std::vector<HDSA::Ptr<const HDSA::Vector<RealT>>> &trans_vec)
        {
            n_t_ = trans_vec.size();
            const_vec_.resize(n_t_);
            for (int k = 0; k < n_t_; k++)
            {
                const_vec_[k] = trans_vec[k];
            }
        }

        ~Transient_Vector_Const()
        {
        }

        int Get_n_t(void) const
        {
            return n_t_;
        }

        //////////////////////////////////////////////////////////////////////////////////
        // Overloading pure virtual functions in HDSA::Vector base class
        //////////////////////////////////////////////////////////////////////////////////

        HDSA::Ptr<HDSA::Vector<RealT>> clone() const override
        {
            HDSA::Ptr<HDSA::Vector<RealT>> vec = HDSA::makePtr<Transient_Vector<RealT>>(n_t_, const_vec_[0]);
            return vec;
        }

        // compute the dot product of this and x
        RealT dot(const HDSA::Vector<RealT> &x) const override
        {
            RealT val = 0.0;
            const Transient_Vector<RealT> x_trans = dynamic_cast<const Transient_Vector<RealT> &>(x);
            for (int k = 0; k < n_t_; k++)
            {
                val += const_vec_[k]->dot(*x_trans.Get_Vector_Const(k));
            }
            return val;
        }

        // compute the dot product of this and x
        RealT norm(void) const override
        {
            RealT val = 0.0;
            for (int k = 0; k < n_t_; k++)
            {
                val += const_vec_[k]->dot(*const_vec_[k]);
            }
            val = std::sqrt(val);
            return val;
        }

        // add alpha*x to this
        void axpy(const RealT alpha, const HDSA::Vector<RealT> &x) override
        {
            HDSA_TEST_FOR_EXCEPTION(true, std::logic_error,
                                    "Error in HDSA::Transient_Vector_Const: Cannot use axpy method" << std::endl);
        }

        // return vector dimension
        int dimension() const override
        {
            return n_t_ * const_vec_[0]->dimension();
        }

        // set this=val elementwise
        void setScalar(const RealT val) override
        {
            HDSA_TEST_FOR_EXCEPTION(true, std::logic_error,
                                    "Error in HDSA::Transient_Vector_Const: Cannot use setScalar method" << std::endl);
        }

        void randomize_standard_normal() override
        {
            HDSA_TEST_FOR_EXCEPTION(true, std::logic_error,
                                    "Error in HDSA::Transient_Vector_Const: Cannot use randomize_standard_normal method" << std::endl);
        }

        void Write_to_File(const std::string &name) const override
        {
            int num_char = name.size();
            std::string name_tmp = name.substr(0, num_char - 4);
            for (int k = 0; k < n_t_; k++)
            {
                std::string name_k = name_tmp + "_time_" + std::to_string(k + 1) + ".txt";
                const_vec_[k]->Write_to_File(name_k);
            }
        }

        //////////////////////////////////////////////////////////////////////////////////
        // Function specific to this class for convenience
        //////////////////////////////////////////////////////////////////////////////////

        HDSA::Ptr<const HDSA::Vector<RealT>> Get_Vector_Const(int k) const override
        {
            return const_vec_[k];
        }
    };

}

#endif
