//
// Created by Steven Maio on 5/22/26.
//

#ifndef OEDLIB_SUBSET_DESIGN_HPP
#define OEDLIB_SUBSET_DESIGN_HPP

#include "OED_Design.hpp"
#include <vector>
#include <iostream>

namespace OED
{
  class Active_Sensors : public Design
  {
  private:
    int total_sensors_;
    std::vector<int> selected_sensors_;

  public:
    Active_Sensors(int total_sensors)
    {
      this->total_sensors_ = total_sensors;
    }

    void Add_Sensor(int v)
    {
      this->selected_sensors_.push_back(v);
    }

    void Remove_Sensor(int v)
    {
      auto ne = std::remove(this->selected_sensors_.begin(), this->selected_sensors_.end(), v);
      this->selected_sensors_.erase(ne);
    }

    void Print_Sensors()
    {
      for (int v: this->selected_sensors_)
      {
        std::cout << v << " ";
      }
      std::cout << std::endl;
    }

    std::vector<int> &Selection()
    {
      return this->selected_sensors_;
    }
  };
}

#endif //OEDLIB_SUBSET_DESIGN_HPP
