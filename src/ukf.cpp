#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

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
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
   is_initialized_ = false;
   previous_timestamp_ = 0;
   
   n_x_ = 5;
   lambda_= 3 - n_x_;
   n_aug_ = 7; 

}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
 
   //Initialization 
   if (!is_initialized_){
	
	if (meas_package.sensor_type_ == MeasurementPackage::LASER){
		x_ << meas_package.raw_measurements_[0],
		      meas_package.raw_measurements_[1],
		      0,
		      0,
		      0;
		
	}
	else if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
                
		double position = meas_package.raw_measurements_[0];
		double angle = meas_package.raw_measurements_[1];
		x_ << position*cos(angle),
		      position*sin(angle),
		      0,
		      0,
		      0;
		
	}

	P_ << 1,0,0,0,0,
	      0,1,0,0,0,
	      0,0,1000,0,0,
	      0,0,0,1000,0,
	      0,0,0,0,1000;

	previous_timestamp_ = meas_package.timestamp_;
	is_initialized_ = true;
	return;
   }

  if (meas_package.sensor_type_ == MeasurementPackage::LASER && !use_laser_){
	return;	
   }
   else if (meas_package.sensor_type_ == MeasurementPackage::RADAR && !use_radar_){
	return;
   }

   //Time between two sensor readings. 
   float delta_t = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
   previous_timestamp_ = meas_package.timestamp_;
   
   Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);
  
   //Weights vetor which will be used to calculate means state and measurement vector and state and measurement covariance matrix
   weights_ = VectorXd(2*n_aug_+1);
   double weight_0 = lambda_/(lambda_+n_aug_);
   weights_(0) = weight_0; 
   for(int i=1; i<2*n_aug_+1;i++){
	weights_(i) = 0.5/(n_aug_+lambda_);
   }

   //Prediction of State at k+1 
   Prediction(delta_t);
   
  
   //Updating the State with sensor readings
   if (meas_package.sensor_type_ == MeasurementPackage::LASER){
	UpdateLidar(meas_package);	
   }
   
   if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
	UpdateRadar(meas_package);
   }
   

}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
	
   //Augumenting State vector and State Covariance matrix 
   VectorXd x_aug_ = VectorXd(n_aug_);
   MatrixXd P_aug_ = MatrixXd(n_aug_, n_aug_);
   

   x_aug_.head(5) = x_;
   x_aug_(5) = 0;
   x_aug_(6) = 0;
  
   P_aug_.fill(0.0);
   P_aug_.topLeftCorner(5,5) = P_;
   P_aug_(5,5) = std_a_*std_a_;
   P_aug_(6,6) = std_yawdd_*std_yawdd_;
   
   //Augumented Sigma Matrix
   MatrixXd xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
   xsig_aug_.col(0) = x_aug_;

   MatrixXd L = P_aug_.llt().matrixL();

   for(int i = 0; i < n_aug_; i++){
	
	xsig_aug_.col(i+1) = x_aug_ + sqrt(lambda_ + n_aug_)*L.col(i);
	xsig_aug_.col(i+n_aug_+1) = x_aug_ - sqrt(lambda_ + n_aug_)*L.col(i);
   }
 
   //Predicted Sigma Matrix (from k to k+1)
   for(int i = 0; i < 2*n_aug_+1; i++){

	VectorXd col = VectorXd(n_x_);
        VectorXd c = VectorXd(n_aug_);
	c = xsig_aug_.col(i);

	double p_x = c(0);
	double p_y = c(1);
	double vel = c(2);
	double yaw = c(3);
	double yaw_rate = c(4);
	double acc = c(5);
	double yaw_acc = c(6);

	if (yaw_rate != 0){
		
		col(0) = p_x + (sin(yaw + yaw_rate*delta_t)-sin(yaw))*vel/yaw_rate + delta_t*delta_t*cos(yaw)*acc/2;
		col(1) = p_y + (-cos(yaw + yaw_rate*delta_t)+cos(yaw))*vel/yaw_rate + delta_t*delta_t*sin(yaw)*acc/2;	
	}
	else{
		col(0) = p_x + vel*cos(yaw)*delta_t + delta_t*delta_t*cos(yaw)*acc/2;
		col(1) = p_y + vel*sin(yaw)*delta_t + delta_t*delta_t*sin(yaw)*acc/2;
	}

	col(2) = vel + acc*delta_t;
	col(3) = yaw + yaw_rate*delta_t + delta_t*delta_t*yaw_acc/2;
	col(4) = yaw_rate + delta_t*yaw_acc;
        
        Xsig_pred_.col(i) = col;
   }

   //Calculating Mean State Vector and State Covariance Matrix 
   VectorXd x_pred = VectorXd(n_x_);
   MatrixXd P_pred = MatrixXd(n_x_,n_x_);

   x_pred.fill(0.0);
   for(int i=0; i<2*n_aug_+1; i++){
	
	x_pred = x_pred + weights_(i)*Xsig_pred_.col(i);
   }   

   P_pred.fill(0.0);
   for(int i=0; i<2*n_aug_+1; i++){
	
	VectorXd x_diff = Xsig_pred_.col(i) - x_pred;

	while(x_diff(3) > M_PI) x_diff(3)-=2.*M_PI;
	while(x_diff(3) < -M_PI) x_diff(3)+=2.*M_PI;

	P_pred = P_pred + weights_(i) * x_diff * x_diff.transpose();
   }
   
   x_ = x_pred;
   P_ = P_pred;
} 

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */

   //Predicting Measurement Sigma matrix from Predicted State Sigma Matrix
   int n_z_ = 2;
   MatrixXd Zsig = MatrixXd(n_z_,2*n_aug_+1);
   for(int i=0; i<2*n_aug_+1;i++){
	
	VectorXd col = 	VectorXd(n_z_);
	VectorXd c = VectorXd(n_x_);
	c = Xsig_pred_.col(i);
	
	double p_x = c(0);
	double p_y = c(1);

	col(0) = p_x;
	col(1) = p_y;
	

	Zsig.col(i) = col;
   }

   //Calculating Mean Measurement Vector 
   VectorXd z_pred = VectorXd(n_z_);
   z_pred.fill(0.0);
   for(int i = 0; i < 2*n_aug_+1;i++){
	
	z_pred = z_pred + weights_(i)*Zsig.col(i);
   }  

   //Calculating Measurement Covariance Matrix
   MatrixXd S = MatrixXd(n_z_,n_z_);
   S.fill(0);
   for(int i = 0; i < 2*n_aug_+1;i++){
	
	VectorXd vec_diff = Zsig.col(i) - z_pred;

	S = S + weights_(i)*vec_diff*vec_diff.transpose();
   }

   MatrixXd R = MatrixXd(n_z_,n_z_);
   R << std_laspx_*std_laspx_, 0,
        0, std_laspy_*std_laspy_;
	
   
   S = S + R;


   MatrixXd Tc = MatrixXd(n_x_,n_z_);


   Tc.fill(0.0);
   for (int i = 0; i < 2*n_aug_+1; i++){

	VectorXd x_diff = Xsig_pred_.col(i) - x_;
	VectorXd z_diff = Zsig.col(i) - z_pred;

	while(x_diff(3)>M_PI) x_diff(3)-=2.*M_PI;
	while(x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

	Tc = Tc + weights_(i)*x_diff*z_diff.transpose();
   }

   //Calculating Kalman Gain
   MatrixXd K = Tc*S.inverse();

   //Data from Sensor
   VectorXd Z = VectorXd(n_z_); 
   Z << meas_package.raw_measurements_[0],
        meas_package.raw_measurements_[1];

   //Updating State Vector and State Covariance Matrix
   x_ = x_ + K*(Z-z_pred);
   P_ = P_ - K * S * K.transpose();
   
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

   //Predicting Measurement Sigma matrix from Predicted State Sigma Matrix
   int n_z_ = 3;
   MatrixXd Zsig = MatrixXd(n_z_,2*n_aug_+1);
   for(int i=0; i<2*n_aug_+1;i++){
	
	VectorXd col = 	VectorXd(n_z_);
	VectorXd c = VectorXd(n_x_);
	c = Xsig_pred_.col(i);
	
	double p_x = c(0);
	double p_y = c(1);
	double vel = c(2);
	double yaw = c(3);
	double yaw_rate = c(4);

	col(0) = sqrt(p_x*p_x + p_y*p_y);
	col(1) = atan2(p_y,p_x);
	col(2) = (p_x*cos(yaw)*vel+p_y*sin(yaw)*vel)/sqrt(p_x*p_x + p_y*p_y);

	Zsig.col(i) = col;
   }
  
   //Calculating Mean Measurement Vector
   VectorXd z_pred = VectorXd(n_z_);
   z_pred.fill(0.0);
   for(int i = 0; i < 2*n_aug_+1;i++){
	
	z_pred = z_pred + weights_(i)*Zsig.col(i);
   }  

   //Calculating Measurement Covariance Matrix
   MatrixXd S = MatrixXd(n_z_,n_z_);
   S.fill(0);
   for(int i = 0; i < 2*n_aug_+1;i++){
	
	VectorXd vec_diff = Zsig.col(i) - z_pred;
	
	while (vec_diff(1)>M_PI) vec_diff(1)-=2.*M_PI;
	while (vec_diff(1)<-M_PI) vec_diff(1)+=2.*M_PI;

	S = S + weights_(i)*vec_diff*vec_diff.transpose();
   }
   

   MatrixXd R = MatrixXd(n_z_,n_z_);
   R << std_radr_*std_radr_, 0, 0,
        0, std_radphi_*std_radphi_, 0,
	0, 0, std_radrd_*std_radrd_;
   
   S = S + R;


   MatrixXd Tc = MatrixXd(n_x_,n_z_);


   Tc.fill(0.0);
   for (int i = 0; i < 2*n_aug_+1; i++){

	VectorXd x_diff = Xsig_pred_.col(i) - x_;
	VectorXd z_diff = Zsig.col(i) - z_pred;

	while(x_diff(3)>M_PI) x_diff(3)-=2.*M_PI;
	while(x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

	while(z_diff(1)>M_PI) z_diff(1)-=2.*M_PI;
	while(z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

	Tc = Tc + weights_(i)*x_diff*z_diff.transpose();
   }
   
   //Calculating Kalman Gain
   MatrixXd K = Tc*S.inverse();

   //Data from Sensor
   VectorXd Z = VectorXd(n_z_); 
   Z << meas_package.raw_measurements_[0],
        meas_package.raw_measurements_[1],
        meas_package.raw_measurements_[2];

   //Updating State Vector and State Covariance Matrix
   x_ = x_ + K*(Z-z_pred);
   P_ = P_ - K * S * K.transpose();
}

























