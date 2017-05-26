#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;


///////////////////////////////////////////////////////////////////////////
namespace
{

inline Eigen::VectorXd hRadar (Eigen::VectorXd const &x)
{
	Eigen::VectorXd y = VectorXd(3); // p, phi, pd

	float const x2 = x[0] * x[0];
	float const y2 = x[1] * x[1];
	float const p = sqrt (x2 + y2);

	y << p,     // p
	     atan2 (x[1], x[0]), // phi
	     (x[0] * x[2] + x[1] * x[3]) / p;

	return y;
}

}


///////////////////////////////////////////////////////////////////////////
KalmanFilter::KalmanFilter()
{
}

KalmanFilter::~KalmanFilter()
{
}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in)
{
	x_ = x_in;
	P_ = P_in;
	F_ = F_in;
	H_ = H_in;
	R_ = R_in;
	Q_ = Q_in;
}

void KalmanFilter::Predict()
{
	/**
	* predict the state
	*/

	x_ = F_ * x_;
	MatrixXd const Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z)
{
	/**
	* update the state by using Kalman Filter equations
	*/

	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long const x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z)
{
	/**
	* update the state by using Extended Kalman Filter equations
	*/

	VectorXd z_pred = hRadar (x_);
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long const x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}
