//
// Created by Steven Maio on 5/22/26.
//

#ifndef OEDLIB_LAZY_GREEDY_HPP
#define OEDLIB_LAZY_GREEDY_HPP

#include <functional>
#include <queue>
#include <math.h>

#include "../base/oed/OED_Active_Sensors.hpp"
#include "Lazy_Eval_Data.hpp"
#include "../base/oed/OED_Discrete_Design_Criterion.hpp"

namespace OED
{
  inline Active_Sensors Lazy_Greedy_Solve(
    Discrete_Design_Criterion &f,
    int total_sensors,
    int budget)
  {
    Active_Sensors sensors(total_sensors);
    std::priority_queue<Lazy_Eval_Data> marginal_gains;
    double objective_value{0};
    for (int v = 0; v < total_sensors; v++)
    {
      Lazy_Eval_Data data(INFINITY, v, -1);
      marginal_gains.push(data);
    }

    // determine optimal sensors
    for (int i = 0; i < budget; i++)
    {
      while (true)
      {
        Lazy_Eval_Data data = marginal_gains.top();
        int v = data.sensor;
        int t = data.time;
        double gain = data.gain;
        marginal_gains.pop();
        sensors.Add_Sensor(v);
        if (t == i)
        {
          objective_value += gain;
          break;
        }
        gain = f.Evaluate(sensors) - objective_value;
        data = Lazy_Eval_Data(gain, v, i);
        marginal_gains.push(data);
        sensors.Remove_Sensor(v);
      }
    }
    return sensors;
  }
}

#endif //OEDLIB_LAZY_GREEDY_HPP
