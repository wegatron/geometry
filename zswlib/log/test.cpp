#include "zsw_log.h"

using namespace zsw;

int main(int argc, char *argv[])
{
  ZSWLOG("zsw_info", "good time");
  NZSWLOG("zsw_info") << "next good time" << std::endl;
  return 0;
}
