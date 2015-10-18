#define BOOST_AUTO_TEST_MAIN

#include "../optimizer.h"

#include <boost/test/included/unit_test.hpp>

using namespace std;

BOOST_AUTO_TEST_SUITE(Inequality_equations_optimizer)

BOOST_AUTO_TEST_CASE(optimizer_0)
{
  Eigen::Matrix<Number,3,1> cx; cx<<1,1,1;
  SmartPtr<zsw::Optimizer> mynlp = new zsw::Optimizer(cx);

  mynlp->addConstraint(1,0,0,0);
  mynlp->addConstraint(0,1,0,0);
  mynlp->addConstraint(0,0,1,0);
  mynlp->addConstraint(-1,-1,-1,-1);

  SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

  app->Options()->SetNumericValue("tol", 1e-3);
  app->Options()->SetStringValue("mu_strategy", "adaptive");
  app->Options()->SetStringValue("output_file", "ipopt.out");

  ApplicationReturnStatus status;
  status = app->Initialize();
  if (status != Solve_Succeeded) {
    printf("\n\n*** Error during initialization!\n");
    return (int) status;
  }

  // Ask Ipopt to solve the problem
  status = app->OptimizeTNLP(mynlp);
  if (status == Solve_Succeeded) {
    printf("\n\n*** The problem solved!\n");
  }
  else {
    printf("\n\n*** The problem FAILED!\n");
  }
  BOOST_REQUIRE(status == Solve_Succeeded);
  BOOST_REQUIRE(mynlp.verify()==true);
}

BOOST_AUTO_TEST_SUITE_END()
