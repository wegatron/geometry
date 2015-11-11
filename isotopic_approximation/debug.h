#ifndef DEBUG_H
#define DEBUG_H

#include <fstream>
#include "basic_data_structure.h"

/// \file debug.h
/// \brief debug functioins
///
///  fuctions for debug
///
/// \author wegatron
/// \email wegatron@hotmail.com
/// \version 1.0
/// \date 2015-11-10


namespace zsw
{
  template< template<typename, typename...> class Container>
    void writeJudgePoints(const std::string &filepath, const Container<JudgePoint> &jpts)
    {
      std::ofstream out_ofs(filepath+"_out.obj");
      std::ofstream in_ofs(filepath+"_in.obj");
      std::ofstream all_ofs(filepath+"_all.obj");
      for(const JudgePoint &jpt : jpts) {
        all_ofs << "v " << jpt.pt_[0] << " " << jpt.pt_[1] << " " << jpt.pt_[2] << std::endl;
        if(jpt.val_exp_>0.5) {
          out_ofs << "v " << jpt.pt_[0] << " " << jpt.pt_[1] << " " << jpt.pt_[2] << std::endl;
        } else {
          in_ofs<< "v " << jpt.pt_[0] << " " << jpt.pt_[1] << " " << jpt.pt_[2] << std::endl;
        }
      }
    }
}

#endif /* DEBUG_H */
