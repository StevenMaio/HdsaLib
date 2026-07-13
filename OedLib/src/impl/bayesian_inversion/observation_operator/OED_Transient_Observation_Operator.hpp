#pragma once

#include "OED_Observation_Operator_Interface.hpp"
#include "OED_Transient_Vector.hpp"
#include "OED_Vector.hpp"
#include "OED_Ptr.hpp"

#include <iostream>

namespace OED
{

    template <class RealT>
    class Transient_Observation_Operator : public Observation_Operator_Interface<RealT>
    {
    private:
        // Transient observation settings
        OED::Ptr<Observation_Operator_Interface<RealT>> obs_;
        int num_meas_;
        std::vector<double> meas_times_;
        std::vector<int> meas_indices_;
        
        //  transient model settings
        int num_steps_;
        double initial_time_{0};
        double final_time_{0};
        double dt_;
    public:
        Transient_Observation_Operator(
            OED::Ptr<Observation_Operator_Interface<RealT>> &obs,
            double initial_time,
            double final_time,
            int num_steps,
            std::vector<double> &meas_times
        ) : obs_(obs), num_meas_(meas_times.size()), meas_times_(meas_times),
            num_steps_(num_steps), initial_time_(initial_time), final_time_(final_time)
        {
            this->dt_ = (final_time - initial_time) / num_steps;

            double t = initial_time;
            int idx = 0;
            for (double tau : meas_times)
            {
                while (t < tau - this->dt_ / 2)
                {
                    t += this->dt_;
                    idx++;
                }
                this->meas_indices_.push_back(idx);
                std::cout << "Transient_Observation_Operator::Constructor meas_idx=" << idx << std::endl;
            }
            this->meas_indices_.resize(this->num_meas_);
        }

        int Data_Dimension() override
        {
            return this->num_meas_ * this->obs_->Data_Dimension();
        }

        int State_Dimension() override
        {
            return this->num_steps_ * this->obs_->State_Dimension();
        }

        void Observation_Operator_Apply(Vector<RealT> &d_out, Vector<RealT> &u_in) override
        {
            auto &d_trans = dynamic_cast<Transient_Vector<RealT> &>(d_out);
            auto &u_trans = dynamic_cast<Transient_Vector<RealT> &>(u_in);
            for (int i = 0; i < this->meas_indices_.size(); i++)
            {
                int t = this->meas_indices_[i];
                this->obs_->Observation_Operator_Apply(*d_trans[i], *u_trans[t]);
            }
        }

        void Observation_Operator_Transpose_Apply(Vector<RealT> &u_out, Vector<RealT> &d_in) override
        {
            auto &d_trans = dynamic_cast<Transient_Vector<RealT> &>(d_in);
            auto &u_trans = dynamic_cast<Transient_Vector<RealT> &>(u_out);
            for (int i = 0; i < this->meas_indices_.size(); i++)
            {
                int t = this->meas_indices_[i];
                this->obs_->Observation_Operator_Transpose_Apply(*u_trans[t], *d_trans[i]);
            }
        }

        Ptr<Vector<RealT>> Get_Empty_Data_Vector() override {
            Ptr<Transient_Vector<RealT>> data_vector = OED::makePtr<Transient_Vector<RealT>>(this->num_meas_, this->obs_->Get_Empty_Data_Vector());
            return data_vector;
        }

    };
}