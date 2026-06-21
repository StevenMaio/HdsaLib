#ifndef OED_MASS_MATRIX_HPP
#define OED_MASS_MATRIX_HPP

#include <memory>

#include "OED_Vector.hpp"
#include "OED_Ptr.hpp"

namespace OED
{
  template <class RealT>
  class Mass_Matrix
  {
  public:

    // Required
    virtual void Apply(Ptr<Vector<RealT>> &x_out,
                       Ptr<Vector<RealT>> &x_in) = 0;

    // Not necessarily required. May depend on application
    virtual void Apply_Inverse(Ptr<Vector<RealT>> &x_out,
                               Ptr<Vector<RealT>> &x_in) {}

    virtual void Apply_Sqrt(Ptr<Vector<RealT>> &x_out,
                            Ptr<Vector<RealT>> &x_in) {}

    virtual void Apply_Inverse_Sqrt(Ptr<Vector<RealT>> &x_out,
                                    Ptr<Vector<RealT>> &x_in) {}

  };
}

#endif // OED_MASS_MATRIX_HPP
