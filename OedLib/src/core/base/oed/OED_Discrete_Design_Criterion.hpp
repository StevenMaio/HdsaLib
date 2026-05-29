//
// Created by Steven Maio on 5/22/26.
//

#ifndef OEDLIB_DISCRETE_DESIGN_CRITERION_HPP
#define OEDLIB_DISCRETE_DESIGN_CRITERION_HPP

#include "OED_Active_Sensors.hpp"

namespace OED
{

    class Discrete_Design_Criterion
    {
    public:
        virtual double Evaluate(Active_Sensors &sensors)
        {
          return 0;
        }

        virtual double Compute_Marginal_Gain(Active_Sensors &sensors, int v)
        {
          return 0;
        }

        virtual ~Discrete_Design_Criterion() {}
    };

}

#endif //OEDLIB_DISCRETE_DESIGN_CRITERION_HPP
