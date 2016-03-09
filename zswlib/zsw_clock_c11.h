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

    class TimeCostMap
    {
    public:
      TimeCostMap() {}
      size_t getTime(const std::string &name) const {
        auto it = time_cost_map_.find(name);
        if(it == time_cost_map_.end()) { return 0; }
        return it->second;
      }

      void setTime(const std::string &name, size_t ms) {
        time_cost_map_[name] = ms;
      }

      void print() const;
    private:
      TimeCostMap(const TimeCostMap &tcm) = delete;
      std::map<std::string, size_t> time_cost_map_; // cost function - millisecond
    };

    class TimeAnalysis
    {
    public:
      static const TimeCostMap &getTCM();
      TimeAnalysis(const std::string &name);
      ~TimeAnalysis();
    private:
      TimeAnalysis(const TimeAnalysis &ta) = delete;
      static TimeCostMap tcm_;
      std::string name_;
      size_t ms_count_;
      ClockC11 clock_;
    };

    TimeCostMap &getTCM();
  }
}

#define FUNCTION_TIME_ANALYSIS() zsw::common::TimeAnalysis tm(__FUNCTION__)
#define PRINT_FUNCTION_COST() zsw::common::TimeAnalysis::getTCM().print();

#endif /* ZSW_CLOCK_C11_H */
