#ifndef OBJECTIVE_SIMOPT_TS_HPP
#define OBJECTIVE_SIMOPT_TS_HPP

template <class Real>
class Objective_SimOpt_TS : public ROL::Objective_SimOpt<Real>
{

protected:
  std::vector<ROL::TimeStamp<Real>> timeStamp_;
  int current_TS_;

public:
  Objective_SimOpt_TS(std::vector<ROL::TimeStamp<Real>> &timeStamp) : ROL::Objective_SimOpt<Real>(), timeStamp_(timeStamp)
  {
    current_TS_ = 0;
  }

  virtual void update(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &t, bool flag = true, int iter = -1)
  {
    Update_current_TS(t);
  }

  void Update_current_TS(Real &t)
  {
    bool success = false;
    for (unsigned int k = 0; k < timeStamp_.size(); k++)
    {
      if (timeStamp_[k].t.at(0) == t)
      {
        current_TS_ = k;
        success = true;
        break;
      }
    }
    if (success == false)
    {
      // std::cout << "The time stamp did not update correctly" << std::endl;
    }
  }
};

#endif
