//
// Created by Steven Maio on 5/27/26.
//

#ifndef OEDLIB_BAYESIAN_INVERSION_HPP
#define OEDLIB_BAYESIAN_INVERSION_HPP
#include "OED_Likelihood_Model.hpp"
#include "OED_Prior_Model.hpp"
#include "../optimization/OED_Constraint.hpp"

namespace OED
{
  template <class RealT>
  class Bayesian_Inversion_Interface
  {
  private:
    // TODO: clean this up -- don't use references
    // Components of Bayesian inverse problem
    Likelihood_Model<RealT> &likelihood_;
    Prior_Model<RealT> &prior_;
    Constraint<RealT> &constraint_;

  public:

    // TODO: do I need this? It's a virtual class, why would I need to create a cosntructor?
    Bayesian_Inversion_Interface(Likelihood_Model<RealT> &likelihood, Prior_Model<RealT> &prior, Constraint<RealT> &constraint)
      : likelihood_(likelihood), prior_(prior), constraint_(constraint) {}

    Likelihood_Model<RealT> &Likelihood()
    {
      return likelihood_;
    }

    Prior_Model<RealT> &Prior()
    {
      return prior_;
    }

    Constraint<RealT> &Constraint()
    {
      return constraint_;
    }

    virtual ~Bayesian_Inversion_Interface() {}

    // Need these to build the relevant things
    // TODO: change these to use smart pointers later
    virtual Vector<RealT> *Get_Empty_Parameter_Vector() = 0;
    virtual Vector<RealT> *Get_Empty_State_Vector() = 0;
    virtual Vector<RealT> *Get_Empty_Data_Vector() = 0;    // TODO: need to eventually think about how to deal with reducing size of data vector
  };
}

#endif //OEDLIB_BAYESIAN_INVERSION_HPP
