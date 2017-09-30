#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  // the filter starts uninitialized
  is_initialized_= false;

  // set nx_x and n_aug
  n_x_ = 5;
  n_aug_ = 7;
  n_sig_ = 2*n_aug_ + 1;

  // initialize sigma point spreading parameter
  lambda_ = 3 - n_x_;

  // set weights
  weights_.fill(1/(2*(lambda_ + n_sig_)));
  weights_(0) = lambda_/(lambda_+n_sig_);
}

UKF::~UKF() {}

/**
 *  First measurement initializer
 */
void UKF::FirstMeasurementInitializer(MeasurementPackage meas_package){
  P_.setZero();
  x_.setZero();
  if (meas_package.sensor_type_ == LASER){
    x_(0) = meas_package.raw_measurements_(0);
    x_(1) = meas_package.raw_measurements_(1);
    P_(0,0) = std_laspx_*2;
    P_(1,1) = std_laspy_*2;
  }
  else {
    double r = meas_package.raw_measurements_(0);
    double phi = meas_package.raw_measurements_(1);
    x_(0) = r*cos(phi);
    x_(1) = r*sin(phi);
    P_(0,0) = sqrt(std_radr_*std_radr_ + std_radphi_*std_radphi_)*2;
    P_(1,1) = P_(0,0);
  }
  P_(2,2) = 1000;
  P_(3,3) = 1000;
  P_(4,4) = 1000;

  time_us_ = meas_package.timestamp_;
  is_initialized_ = true;
  return;
}


/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_){
    FirstMeasurementInitializer(meas_package);




  }

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}
