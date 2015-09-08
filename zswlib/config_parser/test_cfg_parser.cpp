#include <cstdlib>
#include "config_parser.h"

using namespace std;

int main(int argc, char *argv[])
{
  map<string, string> cfg_map;
  int error_code = 0;
  if(error_code = CONFIG_PARSER::parserFile("/home/wegatron/workspace/assure_poj/example.conf", cfg_map))
  {
    exit(error_code);
  }

  for (map<string, string>::iterator it = cfg_map.begin();
       it != cfg_map.end(); ++it)
  {
    cout << it->first << ":" << it->second << endl;
  }
  return 0;
}
