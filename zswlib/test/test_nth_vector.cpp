#include <iostream>
#include <zswlib/config.h>
#include <zswlib/data_type.h>
#include <boost/foreach.hpp>

int main(int argc, char *argv[])
{
  zsw::NthVector<zsw::Scalar> nthv(5, [](zsw::Scalar a, zsw::Scalar b)->bool{ return a<b; });

  for(size_t i=100; i>30; --i) {
    nthv.insert(i);
  }
  const std::list<zsw::Scalar> &data = nthv.getData();
  BOOST_FOREACH(zsw::Scalar tmp, data) {
    std::cerr << tmp << " " << std::endl;
  }

  return 0;
}
