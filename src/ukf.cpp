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

  // set nx_x and n_aug
  n_x_ = 5;
  n_aug_ = 7;
  n_sig_ = 2*n_aug_ + 1;

  // initial state vector
  x_ = VectorXd::Zero(n_x_);

  // initial covariance matrix
  P_ = MatrixXd::Zero(n_x_,n_x_);

  // initial lidar measurement matrix
  H_ = MatrixXd::Zero(2,n_x);
  H_(0,0) = 1;
  H_(1,1) = 1;
  Ht_ = H_.transpose();

  // initial covariance matrix
  Xsig_pred_ = MatrixXd(n_x_, n_sig_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 3;

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

  // initial lidar measurement noise matrix
  Rl_= MatrixXd::Zero(2,2);
  Rl_(0,0) = std_laspx_*std_laspx_;
  Rl_(1,1) = std_laspy_*std_laspy_;

  // initial radar measurement noise matrix
  Rr_= MatrixXd::Zero(3,3);
  Rr_(0,0) = std_radr_*std_radr_;
  Rr_(1,1) = std_radphi_*std_radphi_;
  Rr_(2,2) = std_radrd_*std_radrd_;

  // the filter starts uninitialized
  is_initialized_= false;

  // initialize sigma point spreading parameter
  lambda_ = 3 - n_aug_;
  s_lam_n_x_ = sqrt(lambda_ + n_aug_);

  // set weights
  weights_.fill(1/(2*(lambda_ + n_sig_)));
  weights_(0) = lambda_/(lambda_+n_sig_);

  // set process covariance matrix
  Q_ = MatrixXd::Zero(2,2);
  Q_(0,0) = std_a_*std_a_;
  Q_(1,1) = std_yawdd_*std_yawdd_;
}

UKF::~UKF() {}

/**
 *  First measurement initializer
 */
void UKF::FirstMeasurementInitializer(MeasurementPackage meas_package){

  if (meas_package.sensor_type_ == meas_package.LASER){
    // initialize state and process uncertainty for laser
    x_(0) = meas_package.raw_measurements_(0);
    x_(1) = meas_package.raw_measurements_(1);
    P_(0,0) = std_laspx_*2;
    P_(1,1) = std_laspy_*2;
  }
  else {
    // initialize state and process uncertainty for radar
    double r = meas_package.raw_measurements_(0);
    double phi = meas_package.raw_measurements_(1);
    x_(0) = r*cos(phi);
    x_(1) = r*sin(phi);
    P_(0,0) = sqrt(std_radr_*std_radr_ + std_radphi_*std_radphi_)*2;
    P_(1,1) = P_(0,0);
  }
  // set the rest of the process uncertainty
  P_(2,2) = 1000;
  P_(3,3) = 1000;
  P_(4,4) = 1000;

  //update time step
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
    return;
  }

  // update time and dt
  double dt = (meas_package.timestamp_ - time_us_)/1000000;
  time_us_ = meas_package.timestamp_;

  ///* Prediction
  // if very small dt, skip prediction
  if (dt > 0.001)
    Prediction(dt);

  ///* update measurement
  if (meas_package.sensor_type_== meas_package.LASER){
    UpdateLidar(meas_package);
  }
  else {
    UpdateRadar(meas_package);
  }

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double dt) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  VectorXd x_aug_ = VectorXd::Zero(n_aug_);
  MatrixXd P_aug_ = MatrixXd::Zero(n_aug_,n_aug_);

  // set augmented state vector
  x_aug_.head(n_x_) = x_;

  // set augmented uncertainty matrix
  P_aug_.topLeftCorner(n_x_,n_x_) = P_;
  P_aug_.bottomRightCorner(n_aug_-n_x_,n_aug_-n_x_) = Q_;

  ///* generate sigma points
  MatrixXd Xsig_aug_(n_aug_,n_sig_);
  MatrixXd A_ = sqrt P_aug_.llt().matrixL();
  Xsig_aug_.col(0) = x_aug_;
  for (int i = 0; i < n_aug_; ++i){
    Xsig_aug_.col(i+1) = x_aug_+ s_lam_n_x_*A.col(i);
    Xsig_aug_.col(i+n_aug_) = x_aug_- s_lam_n_x_*A.col(i);
  }

  ///* predict sigma points
  double dt2 = dt*dt;
  double px,py,v,xi,xi_d,nu_a,nu_xi;

  for (int i = 0; i < n_sig_; ++i){
      px = Xsig_aug_(0,i);
      py = Xsig_aug_(1,i);
      v = Xsig_aug_(2,i);
      xi = Xsig_aug_(3,i);
      xi_d = Xsig_aug_(4,i);
      nu_a = Xsig_aug_(5,i);
      nu_xid = Xsig_aug_(6,i);

      if (fabs(xi_d) < 0.001) {
          Xsig_pred_(0,i) = px + v*cos(xi)*dt;
          Xsig_pred_(1,i) = py + v*sin(xi)*dt;
      }
      else {
          Xsig_pred_(0,i) = px + (v/xi_d)*(sin(xi+xi_d*dt)-sin(xi));
          Xsig_pred_(1,i) = py + (v/xi_d)*(-cos(xi+xi_d*dt)+cos(xi));
      }

      Xsig_pred_(0,i) += 0.5*dt2*cos(xi)*nu_a;
      Xsig_pred_(1,i) += 0.5*dt2*sin(xi)*nu_a;
      Xsig_pred_(2,i) = v + dt*nu_a;
      Xsig_pred_(3,i) = xi + xi_d*dt + 0.5*dt2*nu_xid;
      Xsig_pred_(4,i) = xi_d + dt*nu_xid;
  }

  ///* derive mean
  x_.setZero();

  vectorXd err(n_x);
  for (int i = 0; i < n_sig_; ++i){
    x_ += weights_(i)*Xsig_pred_.col(i);
  }

  ///* derive covariance
  P_.setZero();
  vectorXd err(n_x);
  for (int i = 0; i < n_sig_; ++i){
    err = Xsig_pred_.col(i) - x_;
    P_ += weights_* err * err.transpose();
  }
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
  VectorXd y_ = meas_package.raw_measurements_- H_*x_;
  MatrixXd S_ = H_*P_*Ht_ + Rl_;
  MatrixXd K_ = P_*Ht_*S_.inverse();

  x_ += K_*y_;
  P_ = (MatrixXd::Identity(2,2) - K_*H_)*P_;

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
