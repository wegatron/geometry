#ifndef ERROR_CTRL_H
#define ERROR_CTRL_H

#include <iostream>
#include <fstream>
#include <sstream>

#ifndef WIN32
#include <execinfo.h>
#include <unistd.h>
#endif

#include <zswlib/zsw_log.h>

#define OK 0

#ifdef PRINT_BACKTRACE_ACTIVE
#define PRINT_BACKTRACE() do{                                   \
    void *buffer[300];                                          \
    int nbtrace = backtrace(buffer, 100);                       \
    backtrace_symbols_fd(buffer, nbtrace, STDOUT_FILENO);       \
  }while(0)
#else
#define PRINT_BACKTRACE()
#endif

// err_hnd throw exception, or return __LINE__;

#define CHECK_RETCODE(ret_code, err_hnd) do{            \
    if (ret_code != OK) {                               \
      std::stringstream ss;                             \
      ss << "Error occured in " << __FILE__ << __LINE__ \
         << " as ret_code is abnormal!" << std::endl;   \
      ZSWLOG("zsw_err", ss.str());                      \
      PRINT_BACKTRACE();                                \
      err_hnd;                                          \
    }                                                   \
  }while(0)

#define CALL_FUNC(pn, err_hnd) do{                                      \
    ZSWLOG("zsw_info", std::string("Run: ") + #pn + std::string("\n")); \
    if(pn != OK){ ZSWLOG("zsw_err", std::string(#pn)+std::string(" Failed!\n")); PRINT_BACKTRACE(); err_hnd; } \
  }while(0)


#define CALL_FUNC2(ptr, pn, err_hnd) do{                                \
    ptr = pn;                                                           \
    ZSWLOG("zsw_info", std::string("Run: ") + #pn + std::string("\n")); \
    if (ptr == nullptr) {                                               \
      std::stringstream ss;                                             \
      ss <<  __FILE__ << __LINE__                                       \
         << #pn << " return null." << std::endl;                        \
        ZSWLOG("zsw_err", ss.str());                                    \
        PRINT_BACKTRACE();                                              \
        err_hnd;                                                        \
    }                                                                   \
  }while(0)


#define OPEN_STREAM(path, fs, mode, err_hnd) do{                        \
    fs.open(path, mode);                                                \
    if(!fs || !fs.good()) {                                             \
      std::stringstream ss;                                             \
      ss<< __FILE__ << __LINE__ << " could not open file " << path <<  "for " << #mode << std::endl; \
        ZSWLOG("zsw_err", ss.str());                                    \
        PRINT_BACKTRACE();                                              \
        err_hnd;                                                        \
    }                                                                   \
  }while(0)

#if ASSURE_ACTIVE
#define ASSURE(cond, err_hnd) do{                                       \
    if(!(cond)) {                                                       \
      std::stringstream ss;                                             \
      ss << "condition " << #cond << " not satisfied!" << std::endl;    \
        ZSW_LOG("zsw_err", ss.str());                                   \
        PRINT_BACKTRACE();                                              \
        err_hnd;                                                        \
    }                                                                   \
  }while(0)

#define ASSURE_MSG(cond, msg, err_hnd) do{      \
    if(!(cond)) {                               \
      ZSWLOG("zsw_err", msg);                   \
      err_hnd;                                  \
    }                                           \
  }while(0)
#else
#define ASSURE(cond, err_hnd) ((void)(cond))
#define ASSURE_MSG(cond, msg, err_hnd) ((void)(cond))
#endif

#endif /* ERROR_CTRL_H */
