#include "kalman_filter.h"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Predict() {
	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {

  Tools tools;
  VectorXd hx = VectorXd(3);

  hx(0) = sqrt(x_(0)*x_(0) + x_(1)*x_(1));
  hx(1) = atan2(x_(1),x_(0));

  if(sqrt(x_(0)*x_(0) + x_(1)*x_(1)) < 0.0001) {
      hx(2) = 0;
  }
  else {
      hx(2) = (x_(0)*x_(2) + x_(1)*x_(3)) / hx(0);
  }

  //hx(2) = (x_(0)*x_(2) + x_(1)*x_(3)) / hx(0);

  VectorXd y = z - hx;

  //normalizing angles
  while(y(1) >  M_PI) y(1) -= 2*M_PI;
  while(y(1) < -M_PI) y(1) += 2*M_PI;

  MatrixXd Hj = tools.CalculateJacobian(x_);
  MatrixXd Hjt = Hj.transpose();
  MatrixXd S = Hj * P_ * Hjt + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHjt = P_ * Hjt;
  MatrixXd K = PHjt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * Hj) * P_;

}
