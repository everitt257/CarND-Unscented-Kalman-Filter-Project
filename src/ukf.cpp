#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
//#include <math.h>

//#define PI 3.14159265

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
  std_a_ = 0.8; //this could be wildy off, need to reconsider

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.6; //this is could be wildly off, need to reconsider

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
  // condition to determine whether to initialize the data
  is_initialized_ = false;

  // dimension of internal state
  n_x_ = 5;

  // dimension off augmented internal state
  n_aug_ = 7;

  // factor for the weights and scales
  lambda_ = 3 - n_aug_;

  // initialize weights_
  weights_ = VectorXd(15);
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = weight_0;
  for (int i = 1; i < 2 * n_aug_ + 1; i++) 
  {  //2n+1 weights
    double weight = 0.5/(n_aug_ + lambda_);
    weights_(i) = weight;
  }

  // X state sigma points
  Xsig_pred_ = MatrixXd(5,15); 

  // NIS parameters
  NIS_radar_ = 0;
  NIS_laser_ = 0;
  //cout << "Constructed" << endl;
  
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
  if(!is_initialized_){
    cout << "UKF: " << endl;
    //initialize state convariance matrix
    //set P to be the same as in the excercise
    P_ <<  0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
          -0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
           0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
          -0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
          -0.0020,    0.0060,    0.0008,    0.0100,    0.0123;

    time_us_ = meas_package.timestamp_;
    //cout << "Initial time stamp: " << time_us_ << endl;


    if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
      /**
      Convert radar from polar to cartesian coordinates and initialize the state
      */
      double ro = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double ro_dot = meas_package.raw_measurements_[2];
      /**
      need to convert to px, py, v, yaw_angle, yaw_rate 
      */
      //normalize the angle
      while(phi > M_PI) phi -= 2. * M_PI;
      while(phi < -M_PI) phi += 2. * M_PI;
      
      // initialize 5 states
      double px = ro*cos(phi);
      double py = ro*sin(phi);
      double v = ro_dot;
      double yaw_angle = atan2(py,px);
      double yaw_rate = 0.001;
      
      //normalize the angle
      while(yaw_angle > M_PI) yaw_angle -= 2. * M_PI;
      while(yaw_angle < -M_PI) yaw_angle += 2. * M_PI;
      
      cout << "Radar measured" << endl;
      cout << "---------------------------------------------------" << endl;

      x_ << px,py,v,yaw_angle,yaw_rate; 
    }
    else if(meas_package.sensor_type_ == MeasurementPackage::LASER){
      /**
      Initialize state for laser
      */

      //double vx = 0.001;
      //double vy = 0.001;
      // initialize 5 states
      double px = meas_package.raw_measurements_[0];
      double py = meas_package.raw_measurements_[1];
      double v = 0.01;
      double yaw_angle = atan2(py,px);
      double yaw_rate = 0.001;
      //normalize the angle
      while(yaw_angle > M_PI) yaw_angle -= 2. * M_PI;
      while(yaw_angle < -M_PI) yaw_angle += 2. * M_PI;
      cout << "Laser measured" << endl;

      x_ << px,py,v,yaw_angle,yaw_rate;
      //cout << "x_ is: " << endl << x_ << endl;
    }

    is_initialized_ = true;
    //cout << "initialized" << endl;
    return;

  }

  double dt = (meas_package.timestamp_ - time_us_)/1000000.0;
  time_us_ = meas_package.timestamp_;


  cout << "dt: " << dt << endl;
  //cout << "time_us_: " << time_us_ << endl; 

  Prediction(dt);
  //cout << "Predicted" << endl;
  //cout << "x_ is: " << endl << x_ << endl;
  if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
    if(fabs(x_(0) * x_(0) + x_(1) * x_(1)) > 0.0001){
      n_z_ = 3;
      cout << "Radar measured" << endl;
      UpdateRadar(meas_package);
    }
    else
      cout << "Skipped measurement" << endl;
  }
  else{
    n_z_ = 2;
    cout << "Laser measured" << endl;
    UpdateLidar(meas_package);
  }

  // print the output
  cout << "x_ = " << x_(0) << " || " << x_(1) << " || " << x_(2) << " || " << x_(3) << " || " << x_(4) << endl;
  //cout << "P_ = " << "\n" << P_ << endl;

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
  /*
  MatrixXd Xsig = MatrixXd(11,5);
  GenerateSigmaPoints(&Xsig);*/

  MatrixXd Xsig_aug = MatrixXd(7,15);
  AugmentedSigmaPoints(&Xsig_aug);
  //cout << "Augmented" << endl;

  //MatrixXd Xsig_pred = MatrixXd(15,5);
  SigmaPointPrediction(Xsig_aug, delta_t);
  //update internal Xsig_pred
  //Xsig_pred_ = Xsig_pred;
  //cout << "Sigma Predicted" << endl;
  //VectorXd x_pred = VectorXd(5);
  //MatrixXd P_pred = MatrixXd(5,5);
  PredictMeanAndCovariance();
  //cout << "Mean and Covariance Predicted" << endl;
  //update
  //x_ = x_pred;
  //P_ = P_pred;
}

void UKF::GenerateSigmaPoints(MatrixXd* Xsig_out){
  //create sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);
  //calculate square root of P_
  MatrixXd A = P_.llt().matrixL();
  //set first column of sigma point matrix
  Xsig.col(0) = x_;
  //set lambda
  lambda_ = 3 - n_x_;

  //set remaining sigma points
  for (int i = 0; i < n_x_; i++)
  {
    Xsig.col(i+1)     = x_ + sqrt(lambda_ + n_x_) * A.col(i);
    Xsig.col(i+1+n_x_) = x_ - sqrt(lambda_ + n_x_) * A.col(i);
  }
  //write result
  *Xsig_out = Xsig;
}

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out){
  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;
  //cout << "x_ before x_aug: " << endl << x_ << endl;
  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();
  
  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  
  //set lambda
  lambda_ = 3 - n_aug_;

  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }
  *Xsig_out = Xsig_aug;
  //cout << "Xsig_aug: " << endl << Xsig_aug << endl;
}

void UKF::SigmaPointPrediction(MatrixXd & Xsig_aug, double delta_t){
  //create matrix with predicted sigma pints as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
  //cout << "Xsig_aug: " << endl << Xsig_aug << endl;
  //predict sigma points
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
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
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }
    //cout << "v: " << v << endl;
    //cout << "yawd: " << endl;
    //cout << "v/yawd: " << v/yawd << endl;
    //cout << "aprt from py_p: " << sin (yaw + yawd*delta_t) - sin(yaw) << endl;
    //cout << "delta_t: " << delta_t << endl;
    //cout << "loops :" << endl;

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
  }
  //cout << "loops ended" << endl;
  //write result
  Xsig_pred_ = Xsig_pred;
  //cout << "Xsig_pred_: " << endl << Xsig_pred_ << endl;
}

void UKF::PredictMeanAndCovariance(){
  //define spreading parameter
  lambda_ = 3 - n_aug_;

  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);

  //predicted state mean
  x.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x = x + weights_(i) * Xsig_pred_.col(i);
  }
  //cout << "x predicted" << endl << x << endl;
  //cout << "old x: " << endl << x_ << endl;
  //cout << "Xsig_aug: " << endl << Xsig_aug << endl;
  //cout << "Xsig Predi: " << endl << Xsig_pred_ << endl;
  //predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

    P = P + weights_(i) * x_diff * x_diff.transpose() ;
    //cout << "P: " << endl << x_diff * x_diff.transpose() << endl;
  }
  //write result
  //cout << "P predicted" << endl;
  x_ = x;
  P_ = P;
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
  VectorXd z_pred = VectorXd(n_z_);
  MatrixXd S_out = MatrixXd(n_z_, n_z_);
  MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);
  PredictLidarMeasurement(&z_pred, &S_out, &Zsig);
  // get sample measurement
  //cout << "PredictLidarMeasurement Complete" << endl << endl;
  VectorXd z = VectorXd(n_z_);
  z(0) = meas_package.raw_measurements_(0);
  z(1) = meas_package.raw_measurements_(1);
  //cout << "z is: " << endl << z << endl;
  UpdateState(Zsig, z_pred, S_out, z);

}

void UKF::PredictLidarMeasurement(VectorXd* z_out, MatrixXd* S_out, MatrixXd* Zsig_out){

  lambda_ = 3 - n_aug_;
  
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_, 2*n_aug_ + 1);
  
  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);

    // measurement model
    Zsig(0,i) = p_x;
    Zsig(1,i) = p_y;
  }
  //cout << "Zsig is: " << endl << Zsig << endl;
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);
  
  z_pred.fill(0.0);

  for (int i=0; i < 2 * n_aug_ + 1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  //cout << "z_pred is: " << endl << z_pred << endl;
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_, n_z_);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    //while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
    //while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  //cout << "S is: " << endl << S << endl;
  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z_,n_z_);
  
  R << std_laspx_ * std_laspx_, 0, 
       0, std_laspy_ * std_laspy_;

  S = S + R;
  //cout << "S is after R: " << endl << S << endl;
  //write result
  *z_out = z_pred;
  *S_out = S;
  *Zsig_out = Zsig;
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
  VectorXd z_pred = VectorXd(n_z_);
  MatrixXd S_out = MatrixXd(n_z_, n_z_);
  MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);
  PredictRadarMeasurement(&z_pred, &S_out, &Zsig);

  // get sample measurement
  VectorXd z = VectorXd(n_z_);
  z(0) = meas_package.raw_measurements_(0);
  z(1) = meas_package.raw_measurements_(1);
  z(2) = meas_package.raw_measurements_(2);
  UpdateState(Zsig, z_pred, S_out, z);

}

void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out, MatrixXd* Zsig_out) {
  lambda_ = 3 - n_aug_;
  
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_, 2*n_aug_ + 1);
  
  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / Zsig(0,i);                 //r_dot
    
    //cout << "v1 is: " << v1 << endl;
    //cout << "v2 is: " << v2 << endl;
    /*cout << "px is: " << p_x << endl;
    cout << "py is: " << p_y << endl;*/
    //cout << "p_x*v1 + p_y*v2: " << p_x*v1 + p_y*v2 << endl;
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_ + 1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_,n_z_);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  //cout << "S: " << endl << S << endl;
  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z_,n_z_);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;
  
  // write result
  // normalize ?
  //while (z_pred(1) > M_PI) z_pred(1) -= 2. * M_PI;
  //while (z_pred(1) < -M_PI) z_pred(1) += 2. * M_PI;
  *z_out = z_pred;
  *S_out = S;
  *Zsig_out = Zsig;

}

void UKF::UpdateState(MatrixXd & Zsig, VectorXd & z_pred, MatrixXd & S, VectorXd & z) {

  lambda_ = 3 - n_aug_;
  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  // calculate cross corre  Niion matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //cout << "z_diff: " << endl << z_diff << endl;

    if (n_z_ == 3) {
      // angle normalization
      while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
      while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;
    }
    
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
   
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;
    //cout << "x_diff(3): " << x_diff(3) << endl;
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = z - z_pred;

  if (n_z_ == 3) {
    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;   
  }

  // update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  //cout << "K*Z_diff: " << K * z_diff << endl;
  P_ = P_ - K*S*K.transpose();

  // update NIS
  double NIS = z_diff.transpose() * S.inverse() * z_diff;
  if(n_z_ == 3)
    NIS_radar_ =  NIS;
  else
    NIS_laser_ =  NIS;
}
