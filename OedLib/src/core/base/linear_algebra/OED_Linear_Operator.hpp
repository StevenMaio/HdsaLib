/***********************************************************************
 HdsaLib - A library for Hyper-differential Sensitivity Analysis

 Questions? Contact Joseph Hart (joshart@sandia.gov)
************************************************************************/

#ifndef OED_LINEAR_OPERATOR_HPP
#define OED_LINEAR_OPERATOR_HPP

namespace OED
{

  template <class RealT>
  class Linear_Operator
  {

  private:
    bool is_symmetric_;

  public:
    Linear_Operator()
    {
      is_symmetric_ = false;
    }

    virtual ~Linear_Operator()
    {
    }

    void Set_Symmetric(void)
    {
      is_symmetric_ = true;
    }

    bool Is_Symmetric(void) const
    {
      return is_symmetric_;
    }

    // evaluate Apply y=A*x
    virtual void Apply(OED::Vector<RealT> &y, const OED::Vector<RealT> &x) const = 0;
  };

}

#endif
