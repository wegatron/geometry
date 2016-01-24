#include <iostream>
#include "isotopic_debug.h"
using namespace std;

namespace zsw{
  bool isZeroTetExist(const TTds &tds)
  {
    for(auto cit=tds.cells_begin(); cit!=tds.cells_end(); ++cit) {
      if(cit->vertex(0)->info().pt_type_==zsw::ZERO_POINT &&
         cit->vertex(1)->info().pt_type_==zsw::ZERO_POINT &&
         cit->vertex(2)->info().pt_type_==zsw::ZERO_POINT &&
         cit->vertex(3)->info().pt_type_==zsw::ZERO_POINT) {
        return true;
      }
    }
    return false;
  }
}
