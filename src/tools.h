#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

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
	VectorXd CalculateRMSE(vector<VectorXd> const &estimations,
	                       vector<VectorXd> const &ground_truth);

	/**
	* A helper method to calculate Jacobians.
	*/
	MatrixXd CalculateJacobian(VectorXd const &x_state);
};

#endif /* TOOLS_H_ */
