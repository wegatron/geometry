#ifndef ZSW_LOG_H
#define ZSW_LOG_H

#include <fstream>
#include <iostream>
#include <sstream>

#include <memory>
#include <map>
#include <vector>
#include <time.h>
#include <boost/foreach.hpp>

#define ZSW_LOG_ACTIVE

namespace zsw{
  enum STANDARD_OUTPUT{
    COUT = 1,
    CERR = 2
  };

  class CompositOstream {
  public:
  CompositOstream(int standard_output=0) :standard_output_(standard_output) {}

    void setStandardOutput(int standard_output) {
      standard_output_ = standard_output;
    }

    void clear() {
      oss_.clear();
    }

    void addOstream(const std::shared_ptr<std::ostream> &os) {
      oss_.push_back(os);
    }

    template <typename T>
      CompositOstream& operator<<(const T& t) {
      BOOST_FOREACH(std::shared_ptr<std::ostream>& os, oss_) {
        (*os) << t;
      }
      if(standard_output_ & STANDARD_OUTPUT::COUT) {
        std::cout << t;
      }
      if(standard_output_ & STANDARD_OUTPUT::CERR) {
        std::cerr << t;
      }
      return *this;
    }

    CompositOstream& operator<<(std::ostream& (*__pf)(std::ostream&))
      {
        BOOST_FOREACH(std::shared_ptr<std::ostream>& os, oss_) {
          __pf(*os);
        }
        if(standard_output_ & STANDARD_OUTPUT::COUT) {
          __pf(std::cout);
        }
        if(standard_output_ & STANDARD_OUTPUT::CERR) {
          __pf(std::cerr);
        }
        return *this;
      }
    ~CompositOstream() {}
  private:
    int standard_output_;
    CompositOstream(const CompositOstream& cls) {}
    std::vector<std::shared_ptr<std::ostream>> oss_;
  };

  inline std::string getTime()
  {
    time_t now = time(0);
    tm *localtm = localtime(&now);
    std::stringstream ss;
    ss << "[" << localtm->tm_year+1900 << "-" << localtm->tm_mon+1 << "-" << localtm->tm_mday
       << "-" << localtm->tm_hour
       << ":" << localtm->tm_min
       << ":" << localtm->tm_sec << "] ";
    return ss.str();
  }

  class ZswLog{
  public:
    static void resetLogFile(const std::string &log_file) { p_instance = std::shared_ptr<ZswLog>(new ZswLog(log_file)); }
    static std::shared_ptr<ZswLog> getInstance(const std::string& log_type){
      if(p_instance == NULL){ p_instance = std::shared_ptr<ZswLog>(new ZswLog()); }
      p_instance->setCLogtype(log_type);
      return p_instance;
    }

    void setCLogtype(const std::string &clog_type) { clog_type_ = clog_type;  }
    void log(const std::string& info);
    template<typename T>
      CompositOstream& operator << (const T& t)
      {
        os << "# [" << clog_type_ << "] "<< getTime() << t;
        os << t;
        return os;
      }
    CompositOstream& operator<<(std::ostream& (*__pf)(std::ostream&)) {
      os << __pf;
      return os;
    }
  private:
#ifdef WIN32
  ZswLog(const std::string log_file="C:\\users\\wegatron\\AppData\\Local\\Temp\\zswlog.log"):
    isos_ready_(false), os(1),clog_type_("") { init(log_file); }
#else
  ZswLog(const std::string log_file="/tmp/zswlog.log"):
    isos_ready_(false), os(1),clog_type_("") { init(log_file); }
#endif
    void init(const std::string &log_file);
    std::string clog_type_;
    bool isos_ready_;
    CompositOstream os;
    std::map<std::string,int> log_type_map;
    static std::shared_ptr<ZswLog> p_instance;
  };
  typedef std::shared_ptr<ZswLog> pZswLog;

}//end of namespace

#ifdef ZSW_LOG_ACTIVE
#define RESETLOG_FILE(log_file) do{             \
    zsw::ZswLog::resetLogFile(log_file);        \
  }while(0)
#define ZSWLOG(log_type, info) do{                      \
    zsw::pZswLog log = zsw::ZswLog::getInstance(log_type);      \
    log->log(info);                           \
  }while(0)

#define NZSWLOG(log_type) (*zsw::ZswLog::getInstance(log_type))
#else
#define ZSWLOG(log_type, info)
#endif /* ZSW_LOG_ACTIVE */

#endif /*ZSW_LOG_H*/
