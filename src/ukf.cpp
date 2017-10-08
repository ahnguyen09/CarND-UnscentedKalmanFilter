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

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 4;

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

  //set state dimension
  n_x_ = 5;

  //set augmented state dimension
  n_aug_ = 7;

  //set measurement space dimension
  n_z_lidar_ = 2;
  n_z_radar_ = 3;

  //define spreading parameter
  lambda_ = 3 - n_aug_;

  num_sigma_pts_ = 2 * n_aug_ + 1;

  // initial state vector
  x_ = VectorXd(n_x_);
  
  // initial covariance matrix
  P_ = MatrixXd(n_x_,n_x_);

  //set size of sigma point matrix
  Xsig_ = MatrixXd(n_aug_, num_sigma_pts_);

  //set size of matrix with predicted sigma points
  Xsig_pred_ = MatrixXd(n_x_, num_sigma_pts_);

  //set size of weights vector
  weights_ = VectorXd(num_sigma_pts_);

  //set size of matrix for sigma points in measurement space
  Zsig_lidar_ = MatrixXd(n_z_lidar_, num_sigma_pts_);
  Zsig_radar_ = MatrixXd(n_z_radar_, num_sigma_pts_);

  //set size of measurement noise covariance matrix
  R_lidar_ = MatrixXd(n_z_lidar_,n_z_lidar_);
  //initialize measurement noise
  R_lidar_ <<   std_laspx_*std_laspx_,0,
                0,std_laspy_*std_laspy_;

  //set size of measurement noise covariance matrix
  R_radar_ = MatrixXd(n_z_radar_,n_z_radar_);
  //initialize measurement noise
  R_radar_ <<   std_radr_*std_radr_, 0, 0,
                0, std_radphi_*std_radphi_, 0,
                0, 0,std_radrd_*std_radrd_;

  // x,y,v,phi,phi dot
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  // TODO: Initialize P
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

    /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
      * Initialize the state x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
      */

    // first measurement
    // px,py,v,yaw,yaw_dot
    x_ << 0, 0, 0, 0, 0;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
        * Convert radar from polar to cartesian coordinates and initialize state.
        */
      float rho = meas_package.raw_measurements_(0);
      float theta = meas_package.raw_measurements_(1);
      float rho_dot = meas_package.raw_measurements_(2);

      float vx = rho_dot * cos(theta);
			float vy = rho_dot * sin(theta);
      x_(0) = rho * cos(theta);
      x_(1) = rho * sin(theta);
      if (rho_dot == 0.0) {
        x_(2) = 0;
      } else {
        x_(2) = sqrt(vx * vx + vy * vy)*(rho_dot/fabs(rho_dot))/2;
      }
      cout << "radar first" << endl;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
        *Initialize state.
        */
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);
      cout << "x,y first: " << x_(0) << " , " << x_(1);
      if (x_(0) > 0) {
        x_(2) = 0.05;
      } else if (x_(0) < 0) {
        x_(2) = -0.05;
      } else {
        x_(0) = 0.002;
        x_(1) = 0.002;
        x_(2) = 0.0;
      }
      cout << " , " << x_(2) << endl;
      cout << "lidar first" << endl;
    }

    //set the weight vector once. weight vector is only dependent on state dimension which doesn't change
    SetWeights(); 

    // done initializing, no need to predict or update
    previous_timestamp_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  /**
    * Update the state transition matrix F according to the new elapsed time.
    - Time is measured in seconds.
    * Update the process noise covariance matrix.
    * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
    */
  double dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  Prediction(dt);

  /*****************************************************************************
   *  Update
   ****************************************************************************/
  /**
    * Use the sensor type to perform the update step.
    * Update the state and covariance matrices.
    */

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    if (use_radar_) {UpdateRadar(meas_package);}
  } else {
    // Laser updates
    if (use_laser_) {UpdateLidar(meas_package);}
  }

  previous_timestamp_ = meas_package.timestamp_;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double dt) {
  /**
  * 1) generate sigma points
  * 2) use generated sigma points to predict sigma points after time dt
  * 3) calculate mean and covariance from predicted sigma points
  */
  GenerateSigmaPoints(); 
  cout << "dt: " << dt << endl;
  SigmaPointPrediction(dt);
  PredictMeanAndCovariance(); 
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  //set size of predict state vector in radar space 
  VectorXd z_pred_lidar = VectorXd(n_z_lidar_);
  //fill z_pred_radar_ with correct values by calling predict radar measurement
  PredictLidarMeasurement(n_z_lidar_, z_pred_lidar, R_lidar_, Zsig_lidar_); 

  VectorXd z = VectorXd(n_z_lidar_);
  z << meas_package.raw_measurements_;

  UpdateState(n_z_lidar_, z_pred_lidar, z, Zsig_lidar_); 
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  //set size of predict state vector in radar space 
  VectorXd z_pred_radar_ = VectorXd(n_z_radar_);
  //fill z_pred_radar_ with correct values by calling predict radar measurement
  PredictRadarMeasurement(n_z_radar_, z_pred_radar_,R_radar_,Zsig_radar_);

  VectorXd z = VectorXd(n_z_radar_);
  z << meas_package.raw_measurements_;

  UpdateState(n_z_radar_, z_pred_radar_, z, Zsig_radar_);
}

/**
  * Generate sigma points for the prediction step 
  * Augment sigma point matrix with process moise
  */
void UKF::GenerateSigmaPoints() {
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  VectorXd x_aug = VectorXd(n_aug_);

  //init Xsig_ to 0; zero mean error for index 5 and 6
  x_aug.fill(0.0);
  //set mean state to current state
  x_aug.head(n_x_) = x_;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.block(0,0,n_x_,n_x_) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;
  
  //calculate square root of P_aug
  MatrixXd sqrt_P_aug = P_aug.llt().matrixL();

  //caclulate pos and neg of sqrt_P_aug: sqrt(lambda_+n_x)*sqrt_P_aug and -(sqrt(lambda_+n_x))*sqrt_P_aug
  MatrixXd pos_sqrt_P_aug = sqrt(lambda_+n_aug_)*sqrt_P_aug;
  MatrixXd neg_sqrt_P_aug = (-sqrt(lambda_+n_aug_))*sqrt_P_aug;

  //perform colwise add and store back into pos and neg of sqrt_P_aug
  //now pos_sqrt_P_aug hold all x + sqrt(lambda+n_x_) * sqrt_P_aug and 
  //neg_sqrt_P_aug hold all x - sqrt(lambda+n_x_) * sqrt_P_aug and 
  pos_sqrt_P_aug.colwise() += x_aug;
  neg_sqrt_P_aug.colwise() += x_aug;

  //fill correct block of Xsig_
  Xsig_.fill(0.0);
  Xsig_.block(0,0,n_aug_,1) = x_aug;
  Xsig_.block(0,1,n_aug_,n_aug_) = pos_sqrt_P_aug;
  Xsig_.block(0,n_aug_+1,n_aug_,n_aug_) = neg_sqrt_P_aug;
}

void UKF::SigmaPointPrediction(double dt) {
  //predict sigma points
  for (int i = 0; i< num_sigma_pts_; i++)
  {
    //extract values for better readability
    double p_x = Xsig_(0,i);
    double p_y = Xsig_(1,i);
    double v = Xsig_(2,i);
    double yaw = Xsig_(3,i);
    double yawd = Xsig_(4,i);
    double nu_a = Xsig_(5,i);
    double nu_yawdd = Xsig_(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*dt) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*dt) );
    }
    else {
        px_p = p_x + v*dt*cos(yaw);
        py_p = p_y + v*dt*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*dt;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*dt*dt * cos(yaw);
    py_p = py_p + 0.5*nu_a*dt*dt * sin(yaw);
    v_p = v_p + nu_a*dt;

    yaw_p = yaw_p + 0.5*nu_yawdd*dt*dt;
    yawd_p = yawd_p + nu_yawdd*dt;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
}

void UKF::SetWeights() {
  weights_(0) = lambda_/(lambda_+n_aug_);
  for (int i=1; i<num_sigma_pts_; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }
}

/**
* Predict mean and covariance from predicted sigma point matrix
*/
void UKF::PredictMeanAndCovariance() {
  //predicted state mean
  x_.fill(0.0);

  for (int i = 0; i < num_sigma_pts_; i++) {  //iterate over sigma points
    x_ = x_ + weights_(i)*Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < num_sigma_pts_; i++) {  //iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    tools.NormalizeAngle(x_diff(3));

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
}

void UKF::PredictLidarMeasurement(int n_z, VectorXd& z_pred_out, MatrixXd R, MatrixXd& Zsig_out) {
  Zsig_out.fill(0.0);
  //transform sigma points into measurement space
  for (int i = 0; i < num_sigma_pts_; i++) {  //2n+1 simga points

    // extract values for better readibility
    Zsig_out(0,i) = Xsig_pred_(0,i);
    Zsig_out(1,i) = Xsig_pred_(1,i);
  }

  //mean predicted measurement
  z_pred_out.fill(0.0);
  for (int i=0; i < num_sigma_pts_; i++) {
      z_pred_out = z_pred_out + weights_(i) * Zsig_out.col(i);
  }

  //set size of measurement covariance matrix S
  S_ = MatrixXd(n_z_lidar_,n_z_lidar_);
  //measurement covariance matrix S
  S_.fill(0.0);
  for (int i = 0; i < num_sigma_pts_; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig_out.col(i) - z_pred_out;

    S_ = S_ + weights_(i) * z_diff * z_diff.transpose();
  }
  //add measurement noise covariance matrix
  S_ = S_ + R;
}

void UKF::PredictRadarMeasurement(int n_z, VectorXd& z_pred_out, MatrixXd R, MatrixXd& Zsig_out) {
  //transform sigma points into measurement space
  for (int i = 0; i < num_sigma_pts_; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    double p_xx = p_x*p_x;
    double p_yy = p_y*p_y;

    Zsig_out(0,i) = sqrt(p_xx + p_yy);                        //r
    Zsig_out(1,i) = atan2(p_y,p_x);                                 //phi
    if (p_xx + p_yy != 0)
    {
      Zsig_out(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
    } else {
      Zsig_out(2,i) = 0;
    }
  }

  //mean predicted measurement
  z_pred_out.fill(0.0);
  for (int i=0; i < num_sigma_pts_; i++) {
      z_pred_out = z_pred_out + weights_(i) * Zsig_out.col(i);
  }

  S_ = MatrixXd(n_z_radar_,n_z_radar_);
  //measurement covariance matrix S
  S_.fill(0.0);
  for (int i = 0; i < num_sigma_pts_; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig_out.col(i) - z_pred_out;

    //angle normalization
    tools.NormalizeAngle(z_diff(1));

    S_ = S_ + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  S_ = S_ + R;
}

void UKF::UpdateState(int n_z, VectorXd z_pred, VectorXd z, MatrixXd Zsig) {
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < num_sigma_pts_; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization only if it is radar measurement
    if (n_z == n_z_radar_) {tools.NormalizeAngle(z_diff(1));}
    

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization only if it is radar measurement
    if (n_z == n_z_radar_) {tools.NormalizeAngle(x_diff(3));}

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S_.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  tools.NormalizeAngle(z_diff(1));

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S_*K.transpose();
  //calculate NIS and probably write to reference address
}

