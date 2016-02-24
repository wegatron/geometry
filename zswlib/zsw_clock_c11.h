#ifndef ZSW_CLOCK_C11_H
#define ZSW_CLOCK_C11_H

#include <string>
#include <chrono>

namespace zsw {
  namespace common {
    class ClockC11
    {
      typedef  std::chrono::high_resolution_clock::time_point TimePoint;
    public:
      ClockC11() { cur_time_ = start_time_ = std::chrono::high_resolution_clock::now(); }
      void clearAll() { cur_time_ = start_time_ = std::chrono::high_resolution_clock::now(); }
      void clearCur() { cur_time_ = std::chrono::high_resolution_clock::now(); }
      std::string time() {
        TimePoint pre_time = cur_time_;
        cur_time_ =  std::chrono::high_resolution_clock::now();
        return std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(cur_time_ - pre_time).count());
      }
      std::string totalTime() {
        cur_time_ = std::chrono::high_resolution_clock::now();
        return std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(cur_time_ - start_time_).count());
      }
    private:
      TimePoint start_time_;
      TimePoint cur_time_;
    };
  }
}

#endif /* ZSW_CLOCK_C11_H */
