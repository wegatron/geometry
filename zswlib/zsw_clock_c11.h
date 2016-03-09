#ifndef ZSW_CLOCK_C11_H
#define ZSW_CLOCK_C11_H

#include <string>
#include <chrono>
#include <map>

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

      size_t timeCount() {
        TimePoint pre_time = cur_time_;
        cur_time_ =  std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::milliseconds>(cur_time_ - pre_time).count();
      }

      size_t totalTimeCount() {
        cur_time_ = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::milliseconds>(cur_time_ - start_time_).count();
      }

      std::string totalTime() {
        cur_time_ = std::chrono::high_resolution_clock::now();
        return std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(cur_time_ - start_time_).count());
      }
    private:
      ClockC11(const ClockC11 &clock) = delete;
      TimePoint start_time_;
      TimePoint cur_time_;
    };

    class TimeAnalysis
    {
    public:
      TimeAnalysis() {}
      size_t getTime(const std::string &name) const {
        auto it = time_cost_map_.find(name);
        if(it == time_cost_map_.end()) { return 0; }
        return it->second;
      }

      void setTime(const std::string &name, size_t ms) {
        time_cost_map_[name] = ms;
      }
    private:
      TimeAnalysis(const TimeAnalysis &ta) = delete;
      std::map<std::string, size_t> time_cost_map_; // cost function - millisecond
    };

    TimeAnalysis &getTA();
  }
}

#define GET_TA(name) zsw::common::getTA().getTime(name)
#define SET_TA(name, ms) zsw::common::getTA().setTime(name, ms)

#endif /* ZSW_CLOCK_C11_H */
