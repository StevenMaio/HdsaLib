cmake_minimum_required(VERSION 4.2)
project(OedLib)

set(CMAKE_CXX_STANDARD 20)

add_library(OedLib INTERFACE
        src/core/optimization/OED_Constraint.hpp
        src/core/bayesian_inversion/OED_Prior_Model.hpp
        src/core/bayesian_inversion/OED_Likelihood_Model.hpp
        src/core/optimal_experimental_design/OED_Linear_OED_D_Opt.hpp
        src/core/base/oed/OED_Design.hpp
        src/core/base/oed/OED_Active_Sensors.hpp
        src/core/algorithms/OED_Lazy_Greedy.hpp
        src/core/algorithms/Lazy_Eval_Data.hpp
        src/core/base/oed/OED_Design_Criterion.hpp
        src/core/bayesian_inversion/OED_Bayesian_Inversion_Interface.hpp
        src/core/base/linear_algebra/OED_Dense_Matrix.hpp)
