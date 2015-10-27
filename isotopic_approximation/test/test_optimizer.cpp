#define BOOST_AUTO_TEST_MAIN


#include <boost/test/included/unit_test.hpp>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <zswlib/mesh/mesh_type.h>
#include "../optimizer.h"

using namespace std;
using namespace Ipopt;

BOOST_AUTO_TEST_SUITE(Inequality_equations_optimizer)

// BOOST_AUTO_TEST_CASE(optimizer_basic)
// {
//   Eigen::Matrix<Ipopt::Number,3,1> cx; cx<<1,1,1;
//   SmartPtr<zsw::Optimizer> mynlp = new zsw::Optimizer(cx);

//   mynlp->addConstraint(1,0,0,0);
//   mynlp->addConstraint(0,1,0,0);
//   mynlp->addConstraint(0,0,1,0);
//   mynlp->addConstraint(-1,-1,-1,-1);

//   SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

//   app->Options()->SetNumericValue("tol", 1e-3);
//   app->Options()->SetStringValue("mu_strategy", "adaptive");
//   app->Options()->SetStringValue("output_file", "ipopt.out");

//   ApplicationReturnStatus status;
//   status = app->Initialize();
//   if (status != Solve_Succeeded) {
//     printf("\n\n*** Error during initialization!\n");
//     return;
//   }

//   // Ask Ipopt to solve the problem
//   status = app->OptimizeTNLP(mynlp);
//   if (status == Solve_Succeeded) {
//     printf("\n\n*** The problem solved!\n");
//   }
//   else {
//     printf("\n\n*** The problem FAILED!\n");
//   }
//   BOOST_REQUIRE(status == Solve_Succeeded);
//   BOOST_REQUIRE(mynlp->verify()==true);
// }

BOOST_AUTO_TEST_CASE(optimizer_adv0)
{
  const string file_path="/home/wegatron/workspace/geometry/data/sphere.obj";
  zsw::mesh::TriMesh input_mesh;
  if(!OpenMesh::IO::read_mesh(input_mesh, file_path)) {
    std::cerr << "[ERROR] can't read mesh!" << std::endl;
  }
  Eigen::Matrix<zsw::Scalar,3,1> center=Eigen::Matrix<zsw::Scalar,3,1>::Zero();
  for(zsw::mesh::TriMesh::VertexIter vit=input_mesh.vertices_begin(); vit!=input_mesh.vertices_end(); ++vit) {
    center += input_mesh.point(*vit);
  }
  center/=input_mesh.n_vertices();
  std::cerr << "center " << center.transpose() << std::endl;
  Eigen::Matrix<Ipopt::Number,3,1> cx=Eigen::Matrix<Ipopt::Number,3,1>::Random();
  SmartPtr<zsw::Optimizer> mynlp = new zsw::Optimizer(cx);
  for(zsw::mesh::TriMesh::FaceIter fit=input_mesh.faces_begin();
      fit!=input_mesh.faces_end(); ++fit) {
    Eigen::Matrix<zsw::Scalar,3,1> v[3];
    int i=0;
    for(zsw::mesh::TriMesh::FaceVertexIter fvit=input_mesh.fv_iter(*fit); fvit.is_valid(); ++fvit) {
      v[i]=input_mesh.point(*fvit); ++i;
    }
    Eigen::Matrix<zsw::Scalar,3,1> n=(v[1]-v[0]).cross(v[2]-v[0]);
    Eigen::Matrix<zsw::Scalar,3,1> m=center-v[0];
    zsw::Scalar k=n.dot(m);
    Eigen::Matrix<zsw::Scalar,3,1> kn=k*n;
    mynlp->addConstraint(kn[0],kn[1],kn[2], kn.dot(v[0]));
  }

  SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
  app->Options()->SetNumericValue("tol", 1e-3);
  app->Options()->SetStringValue("mu_strategy", "adaptive");
  app->Options()->SetStringValue("output_file", "ipopt.out");

  ApplicationReturnStatus status;
  status = app->Initialize();
  if (status != Solve_Succeeded) {
    printf("\n\n*** Error during initialization!\n");
    return;
  }

  // Ask Ipopt to solve the problem
  status = app->OptimizeTNLP(mynlp);
  if (status == Solve_Succeeded) {
    printf("\n\n*** The problem solved!\n");
  }
  else {
    printf("\n\n*** The problem FAILED!\n");
  }
  Eigen::Matrix<zsw::Scalar,3,1> res;
  mynlp->getResult(res);
  std::cerr << "----------- res:" << res.transpose() << std::endl;
  BOOST_REQUIRE(status == Solve_Succeeded);
  BOOST_REQUIRE(mynlp->verify(center) == true);
  BOOST_REQUIRE(mynlp->verify()==true);
}

BOOST_AUTO_TEST_SUITE_END()
