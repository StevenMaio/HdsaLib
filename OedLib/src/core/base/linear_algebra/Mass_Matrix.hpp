#ifndef OED_MASS_MATRIX_HPP
#define OED_MASS_MATRIX_HPP

#include <memory>
#include "OED_Vector.hpp"

namespace OED
{
  template <class RealT>
  class Mass_Matrix
  {
  public:

    // Required
    virtual void Apply(std::shared_ptr<Vector<RealT>> &x_out,
                       std::shared_ptr<Vector<RealT>> &x_in) = 0;

    // Not necessarily required. May depend on application
    virtual void Apply_Inverse(std::shared_ptr<Vector<RealT>> &x_out,
                               std::shared_ptr<Vector<RealT>> &x_in) {}

    virtual void Apply_Sqrt(std::shared_ptr<Vector<RealT>> &x_out,
                            std::shared_ptr<Vector<RealT>> &x_in) {}

    virtual void Apply_Inverse_Sqrt(std::shared_ptr<Vector<RealT>> &x_out,
                                    std::shared_ptr<Vector<RealT>> &x_in) {}

  };
}

#endif // OED_MASS_MATRIX_HPP
