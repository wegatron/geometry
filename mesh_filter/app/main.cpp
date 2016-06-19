#include <iostream>
#include <boost/algorithm/string/predicate.hpp>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#include "../bnf/command_line_praser.h"
#include "../bnf/zsw_clock.h"

#include "../bnf/bilateral_normal_filter.h"

using namespace std;
using namespace zsw;

int initParam(std::shared_ptr<zsw::BNFParam> param, const zsw::common::CommandLinePraser &clp)
{
  param->st_ = 4;
  param->ut_ = 15;
  param->bs_ = 0.1;
  clp.get<size_t>("st", param->st_); clp.get<size_t>("ut",param->ut_);
  clp.get<zsw::Scalar>("bs", param->bs_);
  clp.get<bool>("need_flip_edge", param->need_flipedge_);
  param->need_smooth_ = false;
  clp.get<bool>("need_smooth", param->need_smooth_);

  if(param->need_smooth_ && !clp.get<zsw::Scalar>("smooth_threshold", param->smooth_threshold_)) {
    std::cerr << "error, smooth_threshold unset!" << std::endl;
    return __LINE__;
  }

  std::string filter_type;
  if(!clp.get<std::string>("filter_type",filter_type)) {
    std::cerr << "filter_type haven't specified!" << std::endl;
    return __LINE__;
  }

  if(filter_type=="bnf_fe") {
    param->bnf_type_=zsw::BNFParam::BNF_TYPE::BASIC_FE_RING;
    if(!clp.get<size_t>("rn",param->rn_)) {
      std::cerr << "number of r ring not specified!" << std::endl;
      return __LINE__;
    }
  } else if(filter_type=="bnf_fv") {
    param->bnf_type_=zsw::BNFParam::BNF_TYPE::BASIC_FV_RING;
    if(!clp.get<size_t>("rn",param->rn_)) {
      std::cerr << "number of r ring not specified!" << std::endl;
      return __LINE__;
    }
  } else if(filter_type=="bnf_fre") {
    // ring adj face, expand withing radius r_max, r_min, normal_threshold;
    param->bnf_type_=zsw::BNFParam::BNF_TYPE::EXTENDED_RADIUS_ERING;
    bool suc=clp.get<zsw::Scalar>("normal_threshold", param->normal_threshold_);
    suc=suc&&(clp.get<size_t>("r_min",param->r_min_));
    suc=suc&&(clp.get<size_t>("r_max",param->r_max_));
    if(!suc) {
      std::cerr << "error, some param unset!" << std::endl;
      return __LINE__;
    }
  } else if(filter_type == "bnf_frv") {
    // ring adj face, expand withing radius r_max, r_min, normal_threshold;
    param->bnf_type_=zsw::BNFParam::BNF_TYPE::EXTENDED_RADIUS_VRING;
    bool suc=clp.get<zsw::Scalar>("normal_threshold", param->normal_threshold_);
    suc=suc&&(clp.get<size_t>("r_min",param->r_min_));
    suc=suc&&(clp.get<size_t>("r_max",param->r_max_));
    if(!suc) {
      std::cerr << "error, some param unset!" << std::endl;
      return __LINE__;
    }
  }
  return 0;
}

void run(const zsw::common::CommandLinePraser &clp)
{
  std::string input_file;
  if(!clp.get<std::string>("input_file",input_file)) {
    std::cerr << "no input file specified!" << std::endl;
    return;
  }
  zsw::common::Clock clock;
  std::shared_ptr<zsw::mesh::TriMesh> tm(new zsw::mesh::TriMesh);
  if(!OpenMesh::IO::read_mesh(*tm,input_file)) {
    std::cerr << "Failed to read mesh: " << input_file << std::endl;
    return;
  }
  std::cerr << "read mesh time: " << clock.time() << std::endl;
  std::cerr << "mesh points:" << tm->n_vertices() << std::endl;
  std::cerr << "mesh faces:" << tm->n_faces() << std::endl;

  std::shared_ptr<zsw::BNFParam> param(new zsw::BNFParam);
  if(!initParam(param, clp)) {
    std::cerr << "init param suc!" << std::endl;
    zsw::common::Clock tmp_clock;
    BilateralNormalFilter bnf(param);
    bnf.run(tm);
    std::cerr << "run_total:" << tmp_clock.time() << std::endl;
    std::string output_file;
    if(!clp.get<std::string>("output_file",output_file)) {
      std::cerr << "no output file specified!" << std::endl;
      return;
    }
    if(OpenMesh::IO::write_mesh(*tm,output_file)) {
      std::cerr << "write "+ output_file + " suc!" << std::endl;
    } else {
      std::cerr << "Failed to write "+ output_file + " suc!" << std::endl;
    }
  } else {
    std::cerr << "initParam failed!" << std::endl;
  }
}

int main(int argc, char *argv[])
{
  zsw::common::CommandLinePraser clp;
  clp.prase(argc, argv);
  zsw::common::Clock clock;
  run(clp);
  std::cout << "Total used: " << clock.totalTime() << std::endl;
  return 0;
}
