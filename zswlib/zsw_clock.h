#ifndef ZSW_CLOCK_H
#define ZSW_CLOCK_H

#include <boost/chrono.hpp>
#include <boost/lexical_cast.hpp>

namespace zsw{
  namespace common {
    class Clock {
      typedef boost::chrono::high_resolution_clock::steady_clock::time_point TimePoint;
    public:
      Clock() {
        cur_time_ = start_time_ = boost::chrono::high_resolution_clock::now();
      }
      std::string time() {
        TimePoint pre_time = cur_time_;
        cur_time_ = boost::chrono::high_resolution_clock::now();
        return boost::lexical_cast<std::string>(boost::chrono::duration_cast<boost::chrono::milliseconds>(cur_time_ - pre_time).count());
      }
      std::string totalTime() {
        cur_time_ = boost::chrono::high_resolution_clock::now();
        return boost::lexical_cast<std::string>(boost::chrono::duration_cast<boost::chrono::milliseconds>(cur_time_ - start_time_).count());
      }
    private:
      TimePoint start_time_;
      TimePoint cur_time_;
    };
  }
}

#endif /* ZSW_CLOCK_H */
