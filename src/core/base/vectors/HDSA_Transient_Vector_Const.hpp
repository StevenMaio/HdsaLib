/***********************************************************************
 HdsaLib - A library for Hyper-differential Sensitivity Analysis
 
 Questions? Contact Joseph Hart (joshart@sandia.gov)
************************************************************************/

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

        HDSA::Ptr<HDSA::Vector<RealT>> Clone() const override
        {
            HDSA::Ptr<HDSA::Vector<RealT>> vec = HDSA::makePtr<Transient_Vector<RealT>>(n_t_, const_vec_[0]);
            return vec;
        }

        // compute the Dot product of this and x
        RealT Dot(const HDSA::Vector<RealT> &x) const override
        {
            RealT val = 0.0;
            const Transient_Vector<RealT> x_trans = dynamic_cast<const Transient_Vector<RealT> &>(x);
            for (int k = 0; k < n_t_; k++)
            {
                val += const_vec_[k]->Dot(*x_trans.Get_Vector_Const(k));
            }
            return val;
        }

        // compute the Dot product of this and x
        RealT Norm(void) const override
        {
            RealT val = 0.0;
            for (int k = 0; k < n_t_; k++)
            {
                val += const_vec_[k]->Dot(*const_vec_[k]);
            }
            val = std::sqrt(val);
            return val;
        }

        // add alpha*x to this
        void Scaled_Plus(const RealT alpha, const HDSA::Vector<RealT> &x) override
        {
            HDSA_TEST_FOR_EXCEPTION(true, std::logic_error,
                                    "Error in HDSA::Transient_Vector_Const: Cannot use Scaled_Plus method" << std::endl);
        }

        // return vector Dimension
        int Dimension() const override
        {
            return n_t_ * const_vec_[0]->Dimension();
        }

        // Set this=val elementwise
        void Set_Scalar(const RealT val) override
        {
            HDSA_TEST_FOR_EXCEPTION(true, std::logic_error,
                                    "Error in HDSA::Transient_Vector_Const: Cannot use Set_Scalar method" << std::endl);
        }

        void Randomize_Standard_Normal() override
        {
            HDSA_TEST_FOR_EXCEPTION(true, std::logic_error,
                                    "Error in HDSA::Transient_Vector_Const: Cannot use Randomize_Standard_Normal method" << std::endl);
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
