#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools()
{
}

Tools::~Tools()
{
}

VectorXd Tools::CalculateRMSE(vector<VectorXd> const &estimations,
                              vector<VectorXd> const &ground_truth)
{
	/**
	TODO:
	* Calculate the RMSE here.
	*/

	VectorXd rmse(4);
	rmse << 0,0,0,0;

	if (estimations.empty () || ground_truth.empty ())
    {
        cerr << "Empty estimation or ground truth vector\n";
        return rmse;
    }

    if (estimations.size () != ground_truth.size ())
    {
        cerr << "Estimation and ground truth vector must be same length\n";
        return rmse;
    }

	// accumulate squared residuals
	for(int i=0; i < estimations.size(); ++i)
	{
		VectorXd const diff = estimations[i] - ground_truth[i];
		VectorXd const squared = diff.array () * diff.array ();
		rmse = rmse + squared;
	}

	rmse = rmse / estimations.size ();

	return rmse.array ().sqrt ();
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state)
{
	/**
	TODO:
	* Calculate a Jacobian here.
	*/
	MatrixXd Hj(3,4);

	//recover state parameters
	float const px = x_state(0);
	float const py = x_state(1);
	float const vx = x_state(2);
	float const vy = x_state(3);

	//pre-compute a set of terms to avoid repeated calculation
	float const c1 = px*px+py*py;
	float const c2 = sqrt(c1);
	float const c3 = (c1*c2);

	//check division by zero
	if(fabs(c1) < 0.0001)
	{
		cerr << "CalculateJacobian () - Error - Division by Zero\n";
		return Hj;
	}

	//compute the Jacobian matrix
	Hj << (px/c2), (py/c2), 0, 0,
		  -(py/c1), (px/c1), 0, 0,
		  py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

	return Hj;
}
