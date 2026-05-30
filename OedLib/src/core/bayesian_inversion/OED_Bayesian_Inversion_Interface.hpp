//
// Created by Steven Maio on 5/27/26.
//

#ifndef OEDLIB_BAYESIAN_INVERSION_HPP
#define OEDLIB_BAYESIAN_INVERSION_HPP

#include <memory>

#include "OED_Likelihood_Interface.hpp"
#include "OED_Prior_Model.hpp"
#include "OED_Constraint_Interface.hpp"

namespace OED
{
  template <class RealT>
  class Bayesian_Inversion_Interface
  {
  private:
    std::shared_ptr<Likelihood_Interface<RealT>> likelihood_;
    std::shared_ptr<Prior_Interface<RealT>> prior_;
    std::shared_ptr<Constraint_Interface<RealT>> constraint_;

  public:

    Bayesian_Inversion_Interface(
        std::shared_ptr<Likelihood_Interface<RealT>> &likelihood,
        std::shared_ptr<Prior_Interface<RealT>> &prior,
        std::shared_ptr<Constraint_Interface<RealT>> &constraint
    )
      : likelihood_(likelihood), prior_(prior), constraint_(constraint) {}

    std::shared_ptr<Likelihood_Interface<RealT>> &Likelihood()
    {
      return likelihood_;
    }

    std::shared_ptr<Prior_Interface<RealT>> &Prior()
    {
      return prior_;
    }

    std::shared_ptr<Constraint_Interface<RealT>> &Constraint()
    {
      return constraint_;
    }

    virtual ~Bayesian_Inversion_Interface() {}

    // Need these to build the relevant things
    // TODO: change these to use smart pointers later
    virtual std::shared_ptr<Vector<RealT>> Get_Empty_Parameter_Vector() = 0;
    virtual std::shared_ptr<Vector<RealT>> Get_Empty_State_Vector() = 0;
    virtual std::shared_ptr<Vector<RealT>> Get_Empty_Data_Vector() = 0;    // TODO: need to eventually think about how to deal with reducing size of data vector
  };
}

#endif //OEDLIB_BAYESIAN_INVERSION_HPP
