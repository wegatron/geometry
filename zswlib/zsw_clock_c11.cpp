#include "zsw_clock_c11.h"
#include <zswlib/zsw_log.h>

zsw::common::TimeCostMap zsw::common::TimeAnalysis::tcm_;

void zsw::common::TimeCostMap::print() const
{
  for(auto it = time_cost_map_.begin(); it != time_cost_map_.end(); ++it) {
    NZSWLOG("func_time_cost") << it->first << " : " << it->second <<std::endl;
  }
}

zsw::common::TimeAnalysis::TimeAnalysis(const std::string &name) : name_(name)
{
  ms_count_ = tcm_.getTime(name);
}

zsw::common::TimeAnalysis::~TimeAnalysis()
{
  tcm_.setTime(name_, ms_count_ + clock_.totalTimeCount());
}

const zsw::common::TimeCostMap &zsw::common::TimeAnalysis::getTCM()
{
  return tcm_;
}
