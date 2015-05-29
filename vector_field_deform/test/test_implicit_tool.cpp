#include <iostream>

#include "../implicit_tools.h"

using namespace std;
using namespace zsw;

void testPlane()
{
  // create integrator
  std::shared_ptr<VectorFieldIntegrator> vf_integrator(new AdVectorIntegrator());
  // create deformer
  std::shared_ptr<VfDeformer> vf_deformer(new VfDeformer());
  vf_deformer->setVectorFieldIntegrator(vf_integrator);
  vf_deformer->loadModel("/home/wegatron/tmp/plane_input.obj");
  // creatr implicit_tool
  double center[3] = {0,0,0};
  SphereDeformTool spdf_tool(center, 0.01, 0.8);
  spdf_tool.setDeformer(vf_deformer);
  spdf_tool.setTimeSlice(10);
  double trans_vec[3][3] = {{0,0.1,0}, {0,0.1,0}, {0,0.1,0}};
  // doing deform
  {
    spdf_tool.translateAndDeform(trans_vec[0]);
    spdf_tool.translateAndDeform(trans_vec[1]);
    spdf_tool.translateAndDeform(trans_vec[2]);
  }
  vf_deformer->saveModel("/home/wegatron/tmp/plane_final.obj");
}

void testSphere()
{
  // create integrator
  std::shared_ptr<VectorFieldIntegrator> vf_integrator(new AdVectorIntegrator());
  // create deformer
  std::shared_ptr<VfDeformer> vf_deformer(new VfDeformer());
  vf_deformer->setVectorFieldIntegrator(vf_integrator);
  vf_deformer->loadModel("/home/wegatron/tmp/sphere_input.obj");
  // creatr implicit_tool
  double center[3] = {0,1,0};
  SphereDeformTool spdf_tool(center, 0.01, 1.0);
  spdf_tool.setDeformer(vf_deformer);
  spdf_tool.setTimeSlice(10);
  double trans_vec[2][3] = {{0,0.1,0}, {0.1, 0,0}};
  size_t time = 20;
  // doing deform
  do {
    spdf_tool.translateAndDeform(trans_vec[0]);
    // spdf_tool.translateAndDeform(trans_vec[1]);
    // spdf_tool.translateAndDeform(trans_vec[2]);
  } while(--time);

  time = 20;
  do {
    spdf_tool.translateAndDeform(trans_vec[1]);
    // spdf_tool.translateAndDeform(trans_vec[1]);
    // spdf_tool.translateAndDeform(trans_vec[2]);
  } while(--time);

  // vf_deformer->saveModel("/home/wegatron/tmp/sphere_final.obj");
}

void testBendBeam()
{
  std::shared_ptr<VectorFieldIntegrator> vf_integrator(new AdVectorIntegrator());
  std::shared_ptr<VfDeformer> vf_deformer(new VfDeformer());
  vf_deformer->setVectorFieldIntegrator(vf_integrator);
  vf_deformer->loadModel("/home/wegatron/tmp/beam_input.obj");
  double b[3] = {1,0,0};
  double a[3] = {0, 1, 0};
  double c[3] = {0,0,0};
  double r[2] = {-2.0, 2.0};
  zsw::BendDeformTool bdf_tool(b,a,c, r[0], r[1]);
  bdf_tool.setDeformer(vf_deformer);
  bdf_tool.rotateAndDeform(3.1415926/2.0);
  vf_deformer->saveModel("/home/wegatron/tmp/beam_bend_output.obj");
}

void testTwistBeam()
{
  std::shared_ptr<VectorFieldIntegrator> vf_integrator(new AdVectorIntegrator());
  std::shared_ptr<VfDeformer> vf_deformer(new VfDeformer());
  vf_deformer->setVectorFieldIntegrator(vf_integrator);
  vf_deformer->loadModel("/home/wegatron/tmp/beam_input.obj");
  double a[3] = {1,0,0};
  double c[3] = {-5.5, 0,0};
  double r[2] = {6.0, 200};// all in inner region
  zsw::TwistDeformTool tdf_tool(a, c, r[0], r[1]);
  tdf_tool.setDeformer(vf_deformer);
  tdf_tool.twistAndDeform(3.1415926);
  vf_deformer->saveModel("/home/wegatron/tmp/beam_twist_output.obj");
}

int main(int argc, char *argv[])
{
  // testPlane();
  // testSphere();
  testBendBeam();
  // testTwistBeam();
  return 0;
}
