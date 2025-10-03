#ifndef HDSA_HESSIAN_INVERSION_HPP
#define HDSA_HESSIAN_INVERSION_HPP

namespace HDSA
{

  template <class RealT>
  class Hessian_Inversion
  {

  private:
    RealT tol_;
    std::string solver_;
    bool verbose_;

  public:
    Hessian_Inversion(RealT tol = 1.e-8, std::string solver = "CG", bool verbose = false)
    {
      tol_ = tol;
      solver_ = solver;
      verbose_ = verbose;
    }

    virtual ~Hessian_Inversion()
    {
    }

    virtual void Apply_RS_Hessian(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in, const HDSA::Vector<RealT> &z) const = 0;

    void Apply_RS_Hessian_Inverse(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in, const HDSA::Vector<RealT> &z) const
    {
      HDSA::Ptr<Hessian_Operator<RealT>> hess_op = HDSA::makePtr<Hessian_Operator<RealT>>(this, &z);
      HDSA::Linear_Algebra::Iterative_Linear_Solve<RealT>(z_out, z_in, *hess_op, tol_, solver_, verbose_);
    }

    template <class ScalarType>
    class Hessian_Operator : public HDSA::Linear_Operator<ScalarType>
    {
    private:
      const HDSA::Hessian_Inversion<ScalarType> *hess_invert_;
      const HDSA::Vector<ScalarType> *z_;

    public:
      Hessian_Operator(const HDSA::Hessian_Inversion<ScalarType> *hess_invert, const HDSA::Vector<RealT> *z) : hess_invert_(hess_invert), z_(z)
      {
      }

      ~Hessian_Operator()
      {
      }

      void Apply(HDSA::Vector<ScalarType> &y, const HDSA::Vector<ScalarType> &x) const
      {
        hess_invert_->Apply_RS_Hessian(y, x, *z_);
      }
    };
  };

}

#endif
