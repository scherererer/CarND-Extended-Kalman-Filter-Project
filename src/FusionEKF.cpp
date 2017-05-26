#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF()
{
	is_initialized_ = false;

	previous_timestamp_ = 0;

	ekf_.x_ = VectorXd(4);

	// initializing matrices
	R_laser_ = MatrixXd(2, 2);
	R_radar_ = MatrixXd(3, 3);
	H_laser_ = MatrixXd(2, 4);
	Hj_ = MatrixXd(3, 4);

	//measurement covariance matrix - laser
	R_laser_ << 0.0225, 0,
	            0,      0.0225;

	//measurement covariance matrix - radar
	R_radar_ << 0.09, 0,      0,
	            0,    0.0009, 0,
	            0,    0,      0.09;

	/**
	* Set the process and measurement noises
	*/

	//measurement matrix
	H_laser_ << 1, 0, 0, 0,
	            0, 1, 0, 0;

	/// \todo I think I might need to initialize this differently if we get radar first
	ekf_.P_ = MatrixXd(4, 4);
	ekf_.P_ << 1, 0, 0, 0,
	           0, 1, 0, 0,
	           0, 0, 1000, 0,
	           0, 0, 0, 1000;

	ekf_.F_ = MatrixXd(4, 4);
	ekf_.F_ << 1, 0, 1, 0,
	           0, 1, 0, 1,
	           0, 0, 1, 0,
	           0, 0, 0, 1;

	ekf_.Q_ = MatrixXd(4, 4);
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF()
{
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack)
{
	/*****************************************************************************
	*  Initialization
	****************************************************************************/
	if (! is_initialized_)
	{
		/**
		* Initialize the state ekf_.x_ with the first measurement.
		* Create the covariance matrix.
		*/

		// first measurement
		ekf_.x_ = VectorXd(4);
		ekf_.x_ << 1, 1, 0, 0;

		previous_timestamp_ = measurement_pack.timestamp_;

		if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
		{
			/**
			Convert radar from polar to cartesian coordinates and initialize state.
			*/
			float const p = measurement_pack.raw_measurements_[0];
			float const phi = measurement_pack.raw_measurements_[1];

			ekf_.x_ << p * cos (phi),  // x
			           p * sin (phi), // y
					   0,
					   0;
		}
		else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
		{
			/**
			Initialize state.
			*/
			ekf_.x_ << measurement_pack.raw_measurements_[0],
			           measurement_pack.raw_measurements_[1],
			           0,
			           0;
		}

		// done initializing, no need to predict or update
		is_initialized_ = true;
		return;
	}

	/*****************************************************************************
	*  Prediction
	****************************************************************************/

	/**
	TODO:
	* Update the state transition matrix F according to the new elapsed time.
	- Time is measured in seconds.
	* Update the process noise covariance matrix.
	* Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
	*/

	double const dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;

	/// \todo Maybe would be interesting to buffer up measurements, and then roll back in time
	/// and calculate forward again if time goes backwards like this -- really depends on if we
	/// have a delayed source though.
	if (dt < 0.0)
	{
		cerr << "Time moved backwards at " << measurement_pack.timestamp_ << "\n";
		return;
	}

	previous_timestamp_ = measurement_pack.timestamp_;

	ekf_.F_(0, 2) = dt;
	ekf_.F_(1, 3) = dt;

	float const dt2 = dt * dt;
	float const dt3 = dt2 * dt;
	float const dt4 = dt3 * dt;

	float const dt3_2 = dt3 / 2.0f;
	float const dt4_4 = dt4 / 4.0f;

	float const noise_ax = 9.0f;
	float const noise_ay = 9.0f;

	ekf_.Q_ << (dt4_4 * noise_ax), 0, (dt3_2 * noise_ax), 0,
	           0, (dt4_4 * noise_ay), 0, (dt3_2 * noise_ay),
	           (dt3_2 * noise_ax), 0, (dt2 * noise_ax), 0,
	           0, (dt3_2 * noise_ay), 0, (dt2 * noise_ay);

	ekf_.Predict();

	/*****************************************************************************
	*  Update
	****************************************************************************/

	/**
	* Use the sensor type to perform the update step.
	* Update the state and covariance matrices.
	*/

	switch (measurement_pack.sensor_type_)
	{
	case MeasurementPackage::RADAR:
		// Radar updates
		ekf_.H_ = tools.CalculateJacobian (ekf_.x_);
		ekf_.UpdateEKF (measurement_pack.raw_measurements_);
		break;
	case MeasurementPackage::LASER:
		// Laser updates
		ekf_.R_ = R_laser_;
		ekf_.H_ = H_laser_;
		ekf_.Update (measurement_pack.raw_measurements_);
		break;
	}

	// print the output
	cout << "x_ = " << ekf_.x_ << endl;
	cout << "P_ = " << ekf_.P_ << endl;
}
