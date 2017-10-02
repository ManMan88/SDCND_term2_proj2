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
  VectorXd rmse = VectorXd::Zero(RMSE_SIZE);
  VectorXd error(RMSE_SIZE);

  for (int i = 0; i < estimations.size(); ++i){
    //compute error
    error = estimations[i] - ground_truth[i];
    error = error.array()*error.array();
    rmse += error;
  }
  //compute rmse
  rmse /= (double)estimations.size();
  rmse = rmse.array().sqrt();
  return rmse;
}
