#ifndef data_TYPE_H
#define data_TYPE_H

#include <list>
#include <set>
#include <iostream>
#include <boost/function.hpp>
#include <Eigen/Dense>

namespace zsw
{
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
        typename std::list<Scalar>::iterator it = data_.begin();
        while(it!=data_.end()) {
          if(val == *it) { return; }
          ++it;
        }
        data_.insert(it,val);
      }

      void erase(Scalar val) {
        for(typename std::list<Scalar>::iterator it = data_.begin();
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

  template <typename Scalar>
    class NthVector
    {
    public:
      NthVector(size_t n, boost::function<bool (Scalar, Scalar)> pre_func) {
        pre_func_ = pre_func;
        n_ = n;
      }
      void insert(Scalar v) {
        typename std::list<Scalar>::reverse_iterator it;
        for(it=data_.rbegin(); it!=data_.rend();
            ++it) {
          if(!pre_func_(v, *it)) { break; }
        }
        data_.insert(it.base(), v);
        if(data_.size()>n_) { data_.pop_back(); }
      }

      size_t size() const { return data_.size(); }
      zsw::Scalar getIthVal(size_t ith) {
        if(ith >=n_ ) {
          std::cerr << "can't get " << ith << "th elemen, exceed " << n_ << std::endl;
        }
        typename std::list<Scalar>::iterator it = data_.begin();
        for(size_t i=0; i<ith; ++i) { ++it; }
        return *it;
      }
      const std::list<Scalar> & getData() const { return data_; }
      void clear() { data_.clear(); }
      Scalar back() {
        if(data_.size()>0) { return data_.back(); }
        return 0;
      }
    private:
      size_t n_;
      boost::function<bool (Scalar, Scalar)> pre_func_;
      std::list<Scalar> data_;
    };
}

#endif /* data_TYPE_H */
