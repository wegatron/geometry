#ifndef _CONFIG_PARSER_H_
#define _CONFIG_PARSER_H_
#include <cstring>
#include <iostream>
#include <fstream>
#include <map>

#include "error.h"

namespace CONFIG_PARSER {

#define MAX_BUFFER_SIZE 50000
  
  int parserFile(const std::string& path, std::map<std::string, std::string>& cfg_map)
  {
  
    std::ifstream fs(path.c_str(), std::ifstream::in);
    if (!fs.is_open())
    {
      std::cerr << __FILE__ << ":" << __LINE__ << ":"
                << "can not open file:"<< path << std::endl;
      return FILEOPEN_ERROR;
    }
  
    char buffer[MAX_BUFFER_SIZE];
    char * str0=buffer;
    char * str1=buffer;
    int line_num = 1;
    while(fs.getline(buffer, MAX_BUFFER_SIZE))
    {
      if(buffer[0] == '#')
        continue;
      str1=std::strchr(buffer, '=');
      if (str1 == NULL)
      {
        std::cerr << __FILE__ << ":" << __LINE__ << ":"
                  << "parse file " << path << "error on line:" << line_num
                  << std::string(buffer) << std::endl;
        return FILEPARSE_ERROR;
      }
      *str1 = '\0';
      cfg_map.insert(
          std::pair<std::string,std::string>(
              std::string(str0),std::string(++str1)
                                             )
                     );
    }
    return 0;
  }
}

#endif
