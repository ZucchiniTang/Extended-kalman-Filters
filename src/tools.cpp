#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;
  // check the validity of the following inputs:
  // * the estimation vector size should not be zero
  // * the estimation vector size should equal ground truth vector size
  if (estimations.size() != ground_truth.size() || estimations.size() == 0) {
    std::cout << "Invalid estimation or ground size " << std::endl;
    return rmse;
  }
  //accumulate squared residuals
  for(unsigned int i = 0; i < estimations.size(); ++i){
    VectorXd residual = estimations[i] - ground_truth[i];
    // coefficient-wise multiplication
    residual = residual.array() * residual.array();
    rmse += residual;
  }
  // Calculate the mean
  rmse = rmse / estimations.size();

  // calculate the squared roor
  rmse = rmse.array().sqrt();
  std::cout << "rmse = " << rmse << std::endl;
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3,4);
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // pre-compute a set of terms to aviod repeated calculation
  float h1 = px*px + py*py;
  float h2 = sqrt(h1);
  float h3 = (h1 * h2);

  // check division by zero
  if(fabs(h1)<0.0001){
    std::cout << "CalculateJacobian() - Error - division by Zero" << std::endl;
    return Hj;
  }

  // calculate Jacobian matrix
  Hj << px/h2, py/h2, 0, 0,
        -px/h1, px/h1, 0 ,0,
        py*(vx*py-vy*px)/h3, px*(vy*px - vx*py)/h3, px/h2, py/h2;

  return Hj;
}
