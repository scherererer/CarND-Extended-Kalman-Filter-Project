#pragma once

#include "Eigen/Dense"

#include <vector>

class Tools
{
public:
	/**
	* Constructor.
	*/
	Tools();

	/**
	* Destructor.
	*/
	virtual ~Tools();

	/**
	* A helper method to calculate RMSE.
	*/
	Eigen::VectorXd CalculateRMSE(std::vector<Eigen::VectorXd> const &estimations,
	                              std::vector<Eigen::VectorXd> const &ground_truth);

	/**
	* A helper method to calculate Jacobians.
	*/
	Eigen::MatrixXd CalculateJacobian(Eigen::VectorXd const &x_state);
};
