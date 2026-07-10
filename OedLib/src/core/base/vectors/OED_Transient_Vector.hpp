#pragma once

#include "OED_Vector.hpp"
#include "OED_Ptr.hpp"

#include <vector>

namespace OED
{
    template <class RealT>
    class Transient_Vector : public Vector<RealT>
    {
    private:
        int num_t_;
        int num_x_;
        std::vector<Ptr<Vector<RealT>>> vec_;

        inline RealT _Get_Entry(int t, int k) const
        {
            return this->vec_[t]->Get_Entry(k);
        }

        inline void _Set_Entry(int t, int k, RealT val)
        {
            this->vec_[t]->Set_Entry(k, val);
        }

    public:
        Transient_Vector()
        {
            // TODO: figure out why this was around in HdsaLib
        }

        Transient_Vector(int num_t, const Ptr<Vector<RealT>> &spatial_vec)
        {
            this->num_t_ = num_t;
            this->num_x_ = spatial_vec->Dimension();
            this->vec_.resize(this->num_t_);
            for (int i = 0; i < this->num_t_; i++)
            {
                vec_[i] = spatial_vec->Clone();
            }
        }

        Transient_Vector(std::vector<Ptr<Vector<RealT>>> &trans_vec)
        {
            this->num_t_ = trans_vec.size();
            this->vec_.resize(this->num_t_);
            for (int i = 0; i < this->num_t_; i++)
            {
                vec_[i] = trans_vec[i];
            }
        }

        int Get_Num_T() const
        {
            return this->num_t_;
        }

        ///////////////////////////////////////////////////
        // Virtual methods specific to transient vectors //
        ///////////////////////////////////////////////////

        Ptr<Vector<RealT>> operator[](int k)
        {
            return this->vec_[k];
        }

        virtual Ptr<Vector<RealT>> Get_Vector_Const(int k) const
        {
            return this->vec_[k];
        }

        ////////////////////////
        // Overridden Methods //
        ////////////////////////

        Ptr<Vector<RealT>> Clone() const override
        {
            auto clone = OED::makePtr<Transient_Vector<RealT>>(this->num_t_, this->vec_[0]);
            return clone;
        }

        RealT Dot(const OED::Vector<RealT> &x) const override
        {
            RealT val = 0.0;
            const Transient_Vector<RealT> &x_trans = dynamic_cast<const Transient_Vector<RealT> &>(x);
            for (int i = 0; i < this->num_t_; i++)
            {
                val += this->Get_Vector_Const(i)->Dot(*x_trans.Get_Vector_Const(i));
            }
            return val;
        }

        void Scaled_Plus(const RealT alpha, const Vector<RealT> &x) override
        {
            const Transient_Vector<RealT> &x_trans = dynamic_cast<const Transient_Vector<RealT> &>(x);
            for (int k = 0; k < this->num_t_; k++)
            {
                vec_[k]->Scaled_Plus(alpha, *x_trans.Get_Vector_Const(k));
            }
        }

        int Dimension() const override
        {
            return this->num_t_ * this->num_x_;
        }

        void Set_Scalar(const RealT val) override
        {
            for (int i = 0; i < this->num_t_; i++)
            {
                this->vec_[i]->Set_Scalar(val);
            }
        }

        void Randomize_Standard_Normal() override
        {
            // TODO: implement this eventually
        }

        RealT Get_Entry(int k) const override
        {
            int t = k / this->num_x_;
            k %= this->num_x_;
            return this->_Get_Entry(t, k);
        }

        void Set_Entry(int k, RealT val) override
        {
            int t = k / this->num_x_;
            k %= this->num_x_;
            this->_Set_Entry(t, k, val);
        }

        inline RealT Get_Entry(int t, int k) const
        {
            return this->_Get_Entry(t, k);
        }

        inline void Set_Entry(int t, int k, RealT val)
        {
            this->_Set_Entry(t, k, val);
        }

    };
}