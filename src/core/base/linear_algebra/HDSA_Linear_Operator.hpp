#ifndef HDSA_LINEAR_OPERATOR_HPP
#define HDSA_LINEAR_OPERATOR_HPP

namespace HDSA
{

  template <class RealT>
  class Linear_Operator
  {

  public:
    Linear_Operator()
    {
    }

    virtual ~Linear_Operator()
    {
    }

    // evaluate matvec y=A*x
    virtual void matvec(HDSA::Vector<RealT> &y, const HDSA::Vector<RealT> &x) const = 0;
  };

}

#endif
