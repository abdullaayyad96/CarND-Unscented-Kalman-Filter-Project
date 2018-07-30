#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
* Initializes Unscented Kalman filter
* This is scaffolding, do not modify
*/
UKF::UKF() {
	// if this is false, laser measurements will be ignored (except during init)
	use_laser_ = true;

	// if this is false, radar measurements will be ignored (except during init)
	use_radar_ = true;

	// initially set to false, set to true in first call of ProcessMeasurement
	is_initialized_ = false;

	//states and augmented states dimensions 
	n_x_ = 5;
	n_aug_ = 7;

	// initial state vector
	x_ = VectorXd(n_x_);

	// initial covariance matrix
	P_ = MatrixXd(n_x_, n_x_);

	// predicted sigma point matrix
	Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

	// Process noise standard deviation longitudinal acceleration in m/s^2
	std_a_ = 1;

	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = 1;

	//DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
	//DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

	//Sigma point spreading parameter
	lambda_ = 3 - n_aug_;

	//inistialize weights vector
	weights_ = VectorXd(2 * n_aug_ + 1);

	for (int i = 0; i < (2*n_aug_+1) ; i++)
	{
		if (i == 0)
			weights_(i) = lambda_ / (lambda_ + n_aug_);
		else
			weights_(i) = 0.5 / (lambda_ + n_aug_);
	}

	// open files to store NIS values
	NIS_radar_f.open("NIS_radar.csv", std::ofstream::out | std::ofstream::trunc);
	NIS_laser_f.open("NIS_laser.csv", std::ofstream::out | std::ofstream::trunc);

	/**
	TODO:

	Complete the initialization. See ukf.h for other member properties.

	Hint: one or more values initialized above might be wildly off...
	*/
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
	if (!is_initialized_) {
		/**
		* Initialize the state ekf_.x_ with the first measurement.
		* Create the covariance matrix.
		*/

		// first measurement

		if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
			/**
			Convert radar from polar to cartesian coordinates and initialize state.
			*/
			// polar coordinates angle
			float phi = meas_package.raw_measurements_[1];

			// check that polar angle is between +-pi
			while (phi > PI) {
				phi -= 2 * PI;
			}
			while (phi  < -PI) {
				phi += 2 * PI;
			}

			// calculate initial states from radar measurments
			x_ << meas_package.raw_measurements_[0] * cos(phi),
				meas_package.raw_measurements_[0] * sin(phi),
				0,
				0,
				0;

			// initialize state covariance matrix
			P_ << 1, 0, 0, 0, 0,
				0, 1, 0, 0, 0,
				0, 0, 1, 0, 0,
				0, 0, 0, 1, 0,
				0, 0, 0, 0, 1;

		}
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
			/**
			Initialize state.
			*/

			// calculate initial states from radar measurments
			x_ << meas_package.raw_measurements_[0],
				meas_package.raw_measurements_[1],
				0,
				0,
				0;

			// initialize state covariance matrix
			P_ << pow(std_laspx_,2), 0, 0, 0, 0,
				0, pow(std_laspy_,2), 0, 0, 0,
				0, 0, 1, 0, 0,
				0, 0, 0, 1, 0,
				0, 0, 0, 0, 1;
		}

		//update timestamp
		time_us_ = meas_package.timestamp_;
		// done initializing, no need to predict or update
		is_initialized_ = true;
		return;
	}

	double dt = (meas_package.timestamp_ - time_us_) / 1000000.0; //dt - expressed in seconds
	time_us_ = meas_package.timestamp_;

	//apply prediction algorithim
	Prediction(dt);

	//apply measurment update algorithim
	if ((meas_package.sensor_type_ == MeasurementPackage::RADAR) && (use_radar_)) {
		UpdateRadar(meas_package);
	}
	else if ((meas_package.sensor_type_ == MeasurementPackage::LASER) && (use_laser_)) {
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
	// create augmented sigma point matrix
	MatrixXd Xsig_aug = MatrixXd::Zero(n_aug_, 2 * n_aug_ + 1);

	//create augmented mean state
	VectorXd x_aug = VectorXd::Zero(n_aug_);
	x_aug.head(n_x_) = x_;

	//create augmented state covariance
	MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);

	//update augmented covariance matrix
	P_aug.topLeftCorner(n_x_, n_x_) = P_;
	P_aug(n_x_, n_x_) = std_a_ * std_a_;
	P_aug(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;

	//update square root matrix
	MatrixXd A = P_aug.llt().matrixL();

	//update augmented sigma points
	Xsig_aug.block(0, 0, n_aug_, 1) = x_aug;
	Xsig_aug.block(0, 1, n_aug_, n_aug_) = x_aug.replicate(1, n_aug_) + std::sqrt(lambda_ + n_aug_) *A;
	Xsig_aug.block(0, n_aug_ + 1, n_aug_, n_aug_) = x_aug.replicate(1, n_aug_) - std::sqrt(lambda_ + n_aug_) *A;

	//predict the sigma points using CTRV model
	for (int i = 0; i<Xsig_aug.cols(); i++)
	{
		float x = Xsig_aug(0, i); //x position
		float y = Xsig_aug(1, i); //y position
		float v = Xsig_aug(2, i); //velocity
		float yaw = Xsig_aug(3, i); //heading angle
		float yawd = Xsig_aug(4, i); //turn rate

		float noise_a = Xsig_aug(5, i); //linear acceleration noise
		float noise_yawdd = Xsig_aug(6, i); //turn rate acceleration noise

		if (yawd != 0)
		{
			Xsig_pred_(0, i) = x + (v / yawd)*(sin(yaw + yawd * delta_t) - sin(yaw)) + 0.5*pow(delta_t, 2)*cos(yaw)*noise_a; //update x position
			Xsig_pred_(1, i) = y + (v / yawd)*(-cos(yaw + yawd * delta_t) + cos(yaw)) + 0.5*pow(delta_t, 2)*sin(yaw)*noise_a; //update y position
		}
		else
		{
			//avoid division by zero
			Xsig_pred_(0, i) = x + delta_t * v*cos(yaw) + 0.5*pow(delta_t, 2)*cos(yaw)*noise_a; //update x position
			Xsig_pred_(1, i) = y + delta_t * v*sin(yaw) + 0.5*pow(delta_t, 2)*sin(yaw)*noise_a; //update y position
		}
		Xsig_pred_(2, i) = v + delta_t * noise_a; //update velocity
		Xsig_pred_(3, i) = yaw + delta_t * yawd + 0.5*pow(delta_t, 2)*noise_yawdd; //update heading position
		Xsig_pred_(4, i) = yawd + delta_t * noise_yawdd; //update turn rate
	}

	//update state mean
	x_.fill(0.0);
	for (int i = 0; i<Xsig_pred_.cols(); i++)
	{
		x_ += weights_(i) * Xsig_pred_.col(i);
	}

	//update state covariance matrix
	P_.fill(0.0);
	for (int i = 0; i<Xsig_pred_.cols(); i++)
	{
		P_ += weights_(i) * (Xsig_pred_.col(i) - x_) * (Xsig_pred_.col(i) - x_).transpose();
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
	//number of measurment variables
	int n_z = 2;

	//create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

	//mean predicted measurement
	VectorXd z_pred = VectorXd::Zero(n_z);

	//innovation covariance matrix S
	MatrixXd S = MatrixXd::Zero(n_z, n_z);

	//measurment covariance matrix
	VectorXd R_vec(n_z);
	R_vec << pow(std_laspx_, 2), pow(std_laspy_, 2);
	MatrixXd R = R_vec.asDiagonal();

	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);

	//transform sigma points into measurement space
	for (int i = 0; i<Zsig.cols(); i++) {
		float x = Xsig_pred_(0, i); //x position
		float y = Xsig_pred_(1, i); //y position

		Zsig(0, i) = x; //predicted x position measurment 
		Zsig(1, i) = y; //predicted x position measurment 
	}

	//calculate mean predicted measurement
	for (int i = 0; i < Zsig.cols(); i++) {  //iterate over sigma points
		z_pred += weights_(i) * Zsig.col(i);
	}

	//calculate innovation covariance matrix S
	for (int i = 0; i < Zsig.cols(); i++) {  //iterate over sigma points
		VectorXd meas_diff = Zsig.col(i) - z_pred;
		S = S + weights_(i) * meas_diff*meas_diff.transpose();
	}
	S += R;

	//create matrix for cross correlation Tc
	for (int i = 0; i < Zsig.cols(); i++)
	{
		Tc += weights_(i)*(Xsig_pred_.col(i) - x_)*(Zsig.col(i) - z_pred).transpose();
	}

	//calculate Kalman gain K;
	MatrixXd K = Tc * S.inverse();

	//update state mean and covariance matrix
	x_ += K * (meas_package.raw_measurements_ - z_pred);
	P_ -= K * S * K.transpose();

	//calculate NIS
	NIS_laser = (meas_package.raw_measurements_ - z_pred).transpose() * S.inverse() * (meas_package.raw_measurements_ - z_pred);
	//store NIS value
	NIS_laser_f << NIS_laser << endl;
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

	//number of measurment variables
	int n_z = 3;

	//create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

	//mean predicted measurement
	VectorXd z_pred = VectorXd::Zero(n_z);

	//innovation covariance matrix S
	MatrixXd S = MatrixXd::Zero(n_z, n_z);

	//measurment covariance matrix
	VectorXd R_vec(n_z);
	R_vec << pow(std_radr_, 2), pow(std_radphi_, 2), pow(std_radrd_, 2);
	MatrixXd R = R_vec.asDiagonal();

	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);


	//transform sigma points into measurement space
	for (int i = 0; i<Zsig.cols(); i++) {
		float x = Xsig_pred_(0, i); //x position
		float y = Xsig_pred_(1, i); //y position
		float v = Xsig_pred_(2, i); //velocity
		float yaw = Xsig_pred_(3, i); //heading angle

		Zsig(0, i) = sqrt(x*x + y * y); //predicted distance measurment 

		Zsig(1, i) = atan2(y, x); //predicted polar angle
								  //normalize angle
		while (Zsig(1, i)> M_PI) Zsig(1, i) -= 2.*PI;
		while (Zsig(1, i)<-M_PI) Zsig(1, i) += 2.*PI;

		Zsig(2, i) = (x*cos(yaw)*v + y * sin(yaw)*v) / (Zsig(0, i)); //predicted velocity parallel to position vector
	}

	//calculate mean predicted measurement
	for (int i = 0; i < Zsig.cols(); i++) {  //iterate over sigma points
		z_pred += weights_(i) * Zsig.col(i);
	}

	//calculate innovation covariance matrix S
	for (int i = 0; i < Zsig.cols(); i++) {  //iterate over sigma points
		VectorXd meas_diff = Zsig.col(i) - z_pred;
		S = S + weights_(i) * meas_diff*meas_diff.transpose();
	}
	S += R;

	//create matrix for cross correlation Tc
	for (int i = 0; i < Zsig.cols(); i++)
	{
		Tc += weights_(i)*(Xsig_pred_.col(i) - x_)*(Zsig.col(i) - z_pred).transpose();
	}

	//calculate Kalman gain K;
	MatrixXd K = Tc * S.inverse();

	//update state mean and covariance matrix
	x_ += K * (meas_package.raw_measurements_ - z_pred);
	P_ -= K * S * K.transpose();

	//calculate NIS
	NIS_radar = (meas_package.raw_measurements_ - z_pred).transpose() * S.inverse() * (meas_package.raw_measurements_ - z_pred);
	//store NIS value
	NIS_radar_f << NIS_radar << endl;
}
