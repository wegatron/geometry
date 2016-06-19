#ifndef data_TYPE_H
#define data_TYPE_H

#include <list>
#include <set>
#include <boost/function.hpp>
#include <Eigen/Dense>

namespace zsw
{
  /***
   * \brief 自定义的set 用来存储每个面片的r环邻域的面片id等.因为std::set耗用内存过大.
   **/
  template <typename Scalar>
    class FakeSet
    {
    public:
      typedef typename std::list<Scalar>::const_iterator const_iterator;
      typedef typename std::list<Scalar>::iterator iterator;
      FakeSet() {}
    template<typename UNIQUE_CONTAINER>
      void initFromSet(const UNIQUE_CONTAINER &in_set) {
        data_.resize(in_set.size());
        std::copy(in_set.cbegin(), in_set.cend(), data_.begin());
      }
      void insert(Scalar val)
      {
        std::list<Scalar>::iterator it = data_.begin();
        while(it!=data_.end()) {
          if(val == *it) { return; }
          ++it;
        }
        data_.insert(it,val);
      }

      void erase(Scalar val) {
        for(std::list<Scalar>::iterator it = data_.begin();
            it!=data_.end(); ++it) {
          if(*it == val) {
            data_.erase(it);
            break;
          }
        }
      }

      const_iterator cbegin() const { return data_.cbegin(); }
      const_iterator cend() const { return data_.cend(); }
      iterator begin() { return data_.begin(); }
      iterator end() { return data_.end(); }
      size_t size() const { return data_.size(); }
    private:
      std::list<Scalar> data_;
    };
}

#endif /* data_TYPE_H */
