#ifndef COMMAND_LINE_PRASER_H
#define COMMAND_LINE_PRASER_H

#include <string>
#include <map>
#include <boost/lexical_cast.hpp>
#include <iostream>

namespace zsw {
  namespace common {
    class CommandLinePraser
    {
    public:
      void prase(int argc, char *argv[]) {
        for(size_t i=1; i<argc; i+=2) {
          vm[std::string(argv[i]+1)] = argv[i+1];
        }
      }
      template<typename T>
        bool get(const std::string &name, T &t) const {
        auto it = vm.find(name);
        if(it == vm.end()) return false;
        try{
          t = boost::lexical_cast<T>(it->second);
        }
        catch(const boost::bad_lexical_cast &) {
          return false;
        }
        return true;
      }
    private:
      std::map<std::string, std::string> vm;
    };
  }
}


#endif /* COMMAND_LINE_PRASER_H */
