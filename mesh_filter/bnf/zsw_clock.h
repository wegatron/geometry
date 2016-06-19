#ifndef ZSW_CLOCK_H
#define ZSW_CLOCK_H

#include <boost/chrono.hpp>
#include <boost/lexical_cast.hpp>

namespace zsw{
  namespace common {
    /***
     * \brief 封装了chrono::high_resolution_clock用来计算程序运行的时间。
     * */
    class Clock {
      typedef boost::chrono::high_resolution_clock::steady_clock::time_point TimePoint;
    public:
      /***
       * \brief 构造函数新建一个时钟记录当前时间点。
       * */
      Clock() {
        cur_time_ = start_time_ = boost::chrono::high_resolution_clock::now();
      }

      /***
       * \brief 计算上时间点到现在的时间，并更新当前时间点.
       * */
      std::string time() {
        TimePoint pre_time = cur_time_;
        cur_time_ = boost::chrono::high_resolution_clock::now();
        return boost::lexical_cast<std::string>(boost::chrono::duration_cast<boost::chrono::milliseconds>(cur_time_ - pre_time).count());
      }

      /***
       * \brief 计算时钟从构件到现在的时间。
       * */
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
