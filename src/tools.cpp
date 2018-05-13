#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
		const vector<VectorXd> &ground_truth){

	VectorXd rmse(4);
	rmse << 0,0,0,0;

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	// ... your code here
    if(estimations.size() == 0 || estimations.size() != ground_truth.size()) {
        std::cout << "CalculateRMSE() - Error - Vector size is zero";
        return rmse;
    }

	//accumulate squared residuals
	for(int i=0; i < estimations.size(); ++i){
        // ... your code here
		VectorXd residual = estimations[i] - ground_truth[i];

		//coefficient-wise multiplication
		residual = residual.array()*residual.array();
		rmse += residual;
	}

	//calculate the mean
	rmse /= estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
}


MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

	MatrixXd Hj(3,4);

	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

    float c1 = px*px + py*py;
    float c2 = sqrt(c1);
    float c3 = py*vx - px*vy;

	//check division by zero
	if(fabs(c1) < 0.0001) {
	    std::cout << "CalculateJabian() - Error - Division by Zero";
	    return Hj;
	}

	//compute the Jacobian matrix
    Hj << px/c2,  py/c2, 0, 0,
          -py/c1, px/c1, 0, 0,
          py*c3/pow(c1, 1.5), -px*c3/pow(c1, 1.5), px/c2, py/c2;

	return Hj;
}
