// 

#ifndef OEDLIB_POISSON_OBS_HPP
#define OEDLIB_POISSON_OBS_HPP

using namespace Eigen;
using namespace OED;

namespace OED_TEST
{

    class Poisson_Forward_Precision
    {
    private:
        std::shared_ptr<Poisson_Constraint> constraint_;
        std::shared_ptr<Poisson_Error_Model> likelihood_;
        std::shared_ptr<Poisson_Prior> prior_;
        int param_dim_;
        int state_dim_;
        int data_dim_;
    public:
        using Scalar = double;
        using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
        using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
        using MapConstVec = Eigen::Map<const Vector>;
        using MapVec = Eigen::Map<Vector>;
        int rows() const { return this->param_dim_; }
        int cols() const { returh this->param_dim_; }

        // TODO: create a constructor

        void perform_op(const Scalar *x_in, Scalar *y_out) const
        {
            // create vectors to wrap input and output
            MapConstVec x(x_in);
            MapVec y(y_out);
            Test_Vector m(this->param_dim_);
            Test_Vector u(this->state_dim_);
            Test_Vector d(this->data_dim_);
            m.Set_Vec(x);

            this->constraint_->State_Solve(u, m);
            this->likelihood_->Observation_Operator_Apply(d, u);
            this->likelihood_->Noise_Covariance_Apply(d, d);

            y.noalias() = d.Vec();
        }
    }
}

#endif // OEDLIB_POISSON_OBS_HPP
