#ifndef HDSA_TIMER_HPP
#define HDSA_TIMER_HPP

#include "Teuchos_Time.hpp"
#include "HDSA_Ptr.hpp"

namespace HDSA
{

  template <class RealT>
  class Timer
  {
  private:
    std::ostream& out_stream_;
    HDSA::Ptr<Teuchos::Time> timer_;

  public:

    Timer(std::ostream& out_stream = std::cout): out_stream_(out_stream)
    {
      timer_ = Teuchos::TimeMonitor::getNewCounter("My timer");
    }

    ~Timer()
    {}

    void Start_Timer(void)
    {
      timer_->start();
    }

    void End_Timer(const std::string & code_description)
    {
      timer_->stop();
      out_stream_ << code_description << " took " << timer_->totalElapsedTime() << " seconds to execute" << std::endl;
    }

   
  };

}

#endif
