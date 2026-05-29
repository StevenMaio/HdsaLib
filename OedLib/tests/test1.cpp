//
// Created by Steven Maio on 5/22/26.
//
#include <iostream>
#include <queue>
#include <math.h>

#include "../src/core/algorithms/OED_Lazy_Greedy.hpp"
#include "../src/core/base/oed/OED_Active_Sensors.hpp"

using namespace OED;

class TestCriterion : public Discrete_Design_Criterion
{
public:
    double Evaluate(Active_Sensors &sensors) override
    {
        double sum = 0;
        for (int v: sensors.selected_sensors_)
        {
            sum += v;
        }
        return sum;
    }

    double Compute_Marginal_Gain(Active_Sensors &sensors, int v) override
    {
        return v;
    }

    ~TestCriterion() override {}
};

int main()
{
    Active_Sensors sensors(5);
    sensors.Add_Sensor(3);
    sensors.Add_Sensor(2);

    TestCriterion f;
    std::cout << f.Evaluate(sensors) << std::endl;

    Active_Sensors optimal_design = Lazy_Greedy_Solve(f, 5, 2);
    std::cout << f.Evaluate(optimal_design) << std::endl;
    return 0;
}
