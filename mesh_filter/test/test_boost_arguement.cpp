#include <iostream>

#include "../common/command_line_praser.h"

using namespace std;

int main(int ac, char *av[])
{
  zsw::common::CommandLinePraser clp;
  clp.prase(ac, av);
  int val = 0;
  if(clp.get<int>("val", val)) {
    std::cerr << val << std::endl;
  }
  return 0;
}
