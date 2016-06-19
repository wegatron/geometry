#include <iostream>
#include <boost/algorithm/string/predicate.hpp>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#include "../bilateral_normal_filter.h"
#include "../data.h"
#include "../common/command_line_praser.h"
#include "../common/zsw_clock.h"

using namespace std;
using namespace zsw;

void filterTrimesh(std::shared_ptr<zsw::data::BNFDataInterFace> df)
{
  BilateralNormalFilter bnf(df);
  bnf.run();

  if(OpenMesh::IO::write_mesh(df->trimesh_,df->output_file_)) {
    std::cout << "write "+ df->output_file_ + " suc!" << std::endl;
  }
}

int initBNDf(std::shared_ptr<zsw::data::BNFDataInterFace> df, const zsw::common::CommandLinePraser &clp)
{
  std::string filter_type;
  clp.get<std::string>("filter_type", filter_type);

  if(!clp.get<std::string>("input_file", df->input_file_)) {
    std::cerr << "haven't specify input_file, -input_file" << std::endl;
    return __LINE__;
  }
  if(!OpenMesh::IO::read_mesh(df->trimesh_, df->input_file_)) {
    std::cerr << "read mesh:" << df->input_file_ << " error!" << std::endl;
    return __LINE__;
  }
  // common param
  df->st_ = 4;
  df->ut_ = 15;
  df->b_s_ = 0.1;
  clp.get<size_t>("st", df->st_); clp.get<size_t>("ut", df->ut_);
  clp.get<zsw::Scalar>("b_s", df->b_s_);
  df->output_file_ = "e:/tmp/out.vtk";
  clp.get<std::string>("output_file", df->output_file_);
  clp.get<bool>("need_flip_edge", df->need_flip_edge_);
  df->need_smooth_ = false;
  clp.get<bool>("need_smooth", df->need_smooth_);

  if(filter_type == "bnf_fe1") {
    // NI: one ring adj face, expand withing edge
    df->adj_ring_->ring_type_ = zsw::mesh::AdjRing::ONE_FE_RING;
  } else if(filter_type == "bnf_fv1") {
    // NII: one ring adj face, expand withing vertex
    df->adj_ring_->ring_type_ = zsw::mesh::AdjRing::ONE_FV_RING;
  } else if(filter_type == "bnf_fe2") {
    // NI: two ring adj face, expand withing edge
    df->adj_ring_->ring_type_ = zsw::mesh::AdjRing::TWO_FE_RING;
  } else if(filter_type == "bnf_fv2") {
    // NII: two ring adj face, expand withing vertex
    df->adj_ring_->ring_type_ = zsw::mesh::AdjRing::TWO_FV_RING;
  } else if(filter_type == "bnf_fre") {
    // ring adj face, expand withing radius r_max, r_min, normal_threshold, smooth_threshold_;
    bool suc = clp.get<zsw::Scalar>("r_min", df->r_min_);
    suc = suc && clp.get<zsw::Scalar>("r_max", df->r_max_);
    suc = suc && clp.get<zsw::Scalar>("normal_threshold", df->normal_threshold_);
    suc = suc && (df->need_smooth_ == clp.get<zsw::Scalar>("smooth_threshold", df->smooth_threshold_));
    if(!suc) {
      std::cerr << "error, some param unset!" << std::endl;
      return __LINE__;
    }
    df->adj_ring_->ring_type_ = zsw::mesh::AdjRing::R_FE_RING;
  } else if(filter_type == "bnf_frv") {
    // ring adj face, expand withing radius r_max, r_min, normal_threshold, smooth_threshold_;
    bool suc = clp.get<zsw::Scalar>("r_min", df->r_min_);
    suc = suc && clp.get<zsw::Scalar>("r_max", df->r_max_);
    suc = suc && clp.get<zsw::Scalar>("normal_threshold", df->normal_threshold_);
    if(!suc) {
      std::cerr << "error, some param unset!" << std::endl;
      return __LINE__;
    }
    df->adj_ring_->ring_type_ = zsw::mesh::AdjRing::R_FV_RING;
  }
  return 0;
}

void run(const zsw::common::CommandLinePraser &clp)
{
  std::string filter_type;
  if(!clp.get<std::string>("filter_type", filter_type)) {
    std::cerr << "have't specify filter_type, -input_file" << std::endl;
    return;
  }
  if(boost::starts_with(filter_type, "bnf")) {
    std::shared_ptr<zsw::data::BNFDataInterFace> df(new zsw::data::BNFDataInterFace());
    if(!initBNDf(df, clp)) {
      filterTrimesh(df);
    }
  } else if(filter_type == "zsw") {
    std::cerr << __FILE__ << " " << __LINE__ << __FUNCTION__ << " haven't implemented!" << std::endl;
    return;
  }
}

int main(int argc, char *argv[])
{
  zsw::common::CommandLinePraser clp;
  clp.prase(argc, argv);
  std::string debug_output_dir = "e:/tmp/";
  if(!clp.get<std::string>("debug_output_dir", debug_output_dir)) {
    std::cerr << "using default debug out put dir:" << debug_output_dir << std::endl;
  }
  zsw::data::setDebugOutputDir(debug_output_dir);
  zsw::common::Clock clock;
  run(clp);
  std::cout << "Total used: " << clock.totalTime() << std::endl;
  return 0;
}
