//
// Created by Steven Maio on 5/22/26.
//

#ifndef OEDLIB_DESIGN_CRITERION_HPP
#define OEDLIB_DESIGN_CRITERION_HPP
#include "OED_Design.hpp"

namespace OED
{

    class DesignCriterion {
    public:
        virtual float Evaluate(Design &design) {
            return 0;
        };

    };

}

#endif //OEDLIB_DESIGN_CRITERION_HPP
