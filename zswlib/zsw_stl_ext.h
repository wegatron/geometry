#ifndef ZSW_STL_EXT_H
#define ZSW_STL_EXT_H

#include <vector>

namespace zsw
{
  template<typename T>
    size_t bsearch(const std::vector<T> vals, T value)
    {
      int l=0;
      int r=vals.size()-1;
      int middle=-1;
      while(l<r) {
        middle=(l+r)>>1;
        if(vals[middle] < value) { l=middle+1; }
        else if(vals[middle] > value) { r=middle-1; }
        else { break; }
      }
      if(middle!=-1 && vals[middle]!=value) { middle=-1; }
      return middle;
    }
}

#endif /* ZSW_STL_EXT_H */
