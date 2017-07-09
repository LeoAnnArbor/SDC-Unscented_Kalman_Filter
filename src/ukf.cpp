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
  std_a_ = 0.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.0;

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

  */
  
  // initialized and timestamp
  is_initialized_ = false;
  previous_timestamp_ = 0;

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = n_x_ + 2;

  // state vector: [p_x p_x vel_abs yaw_angle yaw_rate] in SI units and rad
  x_ = VectorXd(n_x_);

  // // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate nu psi] in SI units and rad
  // x_aug_ = VectorXd(n_aug_);

  // state covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  
  // // aug state covariance matrix
  // P_aug_ = MatrixXd(n_aug_, n_aug_);

  // // Sigma Points
  // Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Initialize the other Parameters from ukf.h
  time_us_ = 0.0;
  
  // Number of sigma points
  n_sig_ = 2 * n_aug_ + 1;
  
  // Weights of sigma points
  weights_ = VectorXd(n_sig_);

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  Xsig_pred_ = MatrixXd(n_x_, n_sig_);

  // Weights of sigma points

  // Measurement noise covariance matrices initialization
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << std_radr_*std_radr_, 0, 0,
              0, std_radphi_*std_radphi_, 0,
              0, 0,std_radrd_*std_radrd_;
  R_lidar_ = MatrixXd(2, 2);
  R_lidar_ << std_laspx_*std_laspx_,0,
              0,std_laspy_*std_laspy_;
}

UKF::~UKF() {}

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
  if (!is_initialized_)
  {
    // first measurement
    cout << "UKF: " << endl;
    
    // initialize state vector
    x_ << 0, 0, 0, 0, 0;

    
    // intialized state covariance matrix P 
    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;
    
    //Initialize weights
    weights_(0) = lambda_/(lambda_ + n_aug_);

    for (int i=1; i<2*n_aug_+1; i++) 
    {  
      weights_(i) = 0.5/(n_aug_+lambda_);
    }

    double p_x = 0;
    double p_y = 0;
    double v = 0;

    if (meas_package.sensor_type_== MeasurementPackage::RADAR) 
    {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double rho_dot = meas_package.raw_measurements_[2];

      // coordinate conversion
      p_x = rho * cos(phi);
      p_y = rho * sin(phi);

      double v_x = rho_dot * cos(phi);
      double v_y = rho_dot * sin(phi);

      v = sqrt(v_x * v_x +  v_y * v_y);

      x_ << p_x, p_y, v, 0, 0;

    }
    else if (meas_package.sensor_type_== MeasurementPackage::LASER) 
    {
      /**
      Initialize state.
      */
      p_x = meas_package.raw_measurements_[0];
      p_y = meas_package.raw_measurements_[1];

      x_ << p_x, p_y, 0, 0, 0;
    }

    // Check for zeros in the initial values of px and py
    if (fabs(x_(0)) < 0.001 && fabs(x_(1)) < 0.001) {
      x_(0) = 0.001;
      x_(1) = 0.001;
    }

    previous_timestamp_ = meas_package.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  } 

  /*****************************************************************************
  *  Prediction
  ****************************************************************************/
  float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = meas_package.timestamp_;

  Prediction(dt);

  /*****************************************************************************
  *  Update
  ****************************************************************************/
  
  // Use the sensor type to perform the update step.
  // Update the state and covariance matrices.
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    
    //Measurement update RADAR
    UpdateRadar(meas_package);
    
  } else {
    
    //Measurement update LIDAR
    UpdateLidar(meas_package);
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

  //create augmented mean state
  VectorXd x_aug_ = VectorXd(n_aug_);
  // Augmented state covarience matrix
  MatrixXd P_aug_ = MatrixXd(n_aug_, n_aug_);
  
  x_aug_.head(5) = x_;
  x_aug_(5) = 0;
  x_aug_(6) = 0;
  
  //create augmented covariance matrix
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(n_x_,n_x_) = P_;
  P_aug_(5,5) = std_a_*std_a_;
  P_aug_(6,6) = std_yawdd_*std_yawdd_;
  
  //create square root matrix
  MatrixXd L = P_aug_.llt().matrixL();
  
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  
  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug_;

  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i + 1) = x_aug_ + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug_ - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  //predict sigma points
  for (int i = 0; i< 2 * n_aug_ + 1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd * delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw + yawd * delta_t) );
    }
    else {
        px_p = p_x + v * delta_t * cos(yaw);
        py_p = p_y + v * delta_t * sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5 * nu_a*delta_t * delta_t * cos(yaw);
    py_p = py_p + 0.5 * nu_a*delta_t * delta_t * sin(yaw);
    v_p = v_p + nu_a * delta_t;

    yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p = yawd_p + nu_yawdd * delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) 
  {  //iterate over sigma points
    x_ = x_+ weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) 
  {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_= P_ + weights_(i) * x_diff * x_diff.transpose() ;
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

   //set measurement dimension, px, py
  int n_z = 2;
  
  //set the measurements
  VectorXd z = meas_package.raw_measurements_;
  
  //initialize Zsig
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  Zsig.fill(0.0);
  
  //initialize S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  
  //initialize measurement prediction vector
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //sigma point predictions in process space
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);

    //sigma point predictions in measurement space
    Zsig(0,i) = px;                   
    Zsig(1,i) = py;                                 
  }

  //mean predicted measurement
  z_pred = Zsig * weights_; // Matrix Calc (2X15) * (15X1) = (2X1)

  //measurement covariance matrix S
  for (int i = 0; i < 2 * n_aug_ + 1; i++) 
  { 
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_laspx_*std_laspx_, 0,
          0, std_laspy_*std_laspy_;
  S = S + R;
  
  //Now Update the State x_ and state covariance P_
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  
  //calculate cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
  
  //NIS Lidar Update
  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff; 


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

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;
  
  //set the measurements
  VectorXd z = meas_package.raw_measurements_;
  
  //initialize Zsig
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  Zsig.fill(0.0);
  
  //initialize S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  
  //initialize measurement prediction vector
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  
  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //extract values for better readibility 
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;
    
    //check for zeros
    if (fabs(p_x) < 0.001) {
     p_x = 0.001;
    }
    if (fabs(p_y) < 0.001) {
     p_y = 0.001;
    }

    //measurement model (sigma point predictions in measurement space)
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  z_pred = Zsig * weights_; //Matrix Calc (3X15) * (15X1) = (3X1)

  //measurement covariance matrix S
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0, std_radrd_*std_radrd_;
  S = S + R;
  
  //Now Update the State x_ and state covariance P_
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  
  //calculate cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
  
  //NIS Update
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
}
