#include <iostream>
#include "tools.h"
#define RMSE_SIZE 4

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
  VectorXd rmse(RMSE_SIZE);
  VectorXd error(RMSE_SIZE);
  rmse.setZero();

  for (int i = 0; i < estimations.size(); ++i){
    // transform the estimations state to the ground truth state (px, py, vx, vy)
    error(0) = estimations[i](0);
    error(1) = estimations[i](1);
    error(2) = estimations[i](2) * cos(estimations[i](3);
    error(3) = estimations[i](2) * sin(estimations[i](3));

    //compute error
    error -= ground_truth[i];
    rmse += error.array()*error.array();
  }
  //compute rmse
  rmse /= estimations.size();
  rmse = rmse.array().sqrt();
  return rmse;
}
