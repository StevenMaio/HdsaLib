/***********************************************************************
 HdsaLib - A library for Hyper-differential Sensitivity Analysis
 
 Questions? Contact Joseph Hart (joshart@sandia.gov)
************************************************************************/

#ifndef OED_TIMER_HPP
#define OED_TIMER_HPP

#include "Teuchos_Time.hpp"
#include "OED_Ptr.hpp"

namespace OED
{

  namespace Trilinos_Adapter
  {

  template <class RealT>
  class Timer
  {
  private:
    std::ostream& out_stream_;
    OED::Ptr<Teuchos::Time> timer_;

  public:

    Timer(std::ostream& out_stream = std::cout): out_stream_(out_stream)
    {
      timer_ = Teuchos::TimeMonitor::getNewCounter("My timer");
    }

    ~Timer()
    {}

    void Start_Timer(void)
    {
      timer_->reset();
      timer_->start();
    }

    RealT End_Timer(void)
    {
      timer_->stop();
      RealT elapsed_time = timer_->totalElapsedTime();
      return elapsed_time;
    }

   
  };

  }

}

#endif
