#include "zsw_log.h"

#include <iostream>
#include <string>
#include <time.h>
#include <utility>

using namespace std;
using namespace zsw;

#ifdef WIN32
static const string log_path_prefix="C:\\users\\wegatron\\AppData\\Local\\Temp\\zswlog";
#else
static const string log_path_prefix="/tmp/zswlog";
#endif /* WIN32 */

static pair<string, int> log_types[] = {
  pair<string, int>("zsw_info",1),
  pair<string, int>("zsw_err",1)
};

pZswLog ZswLog::p_instance;

void zsw::ZswLog::init()
{
  int log_type_size = sizeof(log_types)/sizeof(pair<string,int>);
  for (int i=0; i<log_type_size; ++i)
  {
	  log_type_map.insert(log_types[i]);
  }
  string file_path = log_path_prefix + ".log";
  std::shared_ptr<std::ofstream> ofs_ptr(new ofstream(file_path.c_str(), std::fstream::app));
  if(!ofs_ptr || !(*ofs_ptr)|| !ofs_ptr->good()) {
    cerr << "[WARNING] can not open log file: " << file_path << "to write!" << endl;
    return;
  }
  os.addOstream(ofs_ptr);
  isos_ready_ = true;
}

void zsw::ZswLog::log(const string& info)
{
  map<string,int>::const_iterator it = log_type_map.find(clog_type_);
  if (!isos_ready_ || it == log_type_map.end() || it->second==0) {
    return;
  }
  os << "# [" << clog_type_ << "] "<< getTime() << info << endl;
}
