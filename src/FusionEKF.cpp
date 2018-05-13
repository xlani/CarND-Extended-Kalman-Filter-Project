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
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  //Initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //Measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //Measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /** Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

  //Measurement matrix for laser measurements
  ekf_.H_ = MatrixXd(2, 4);
  ekf_.H_ << 1, 0, 0, 0,
             0, 1, 0, 0;

  //The initial transition matrix F_
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0,
             0, 0, 0, 1;

  //Initializing process covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_.setZero();

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {

    /** Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    //First measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      float rhod = measurement_pack.raw_measurements_[2];
      ekf_.x_ << sin(phi)*rho, cos(phi)*rho, 1, 1;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 1, 1;
    }

    //Create and initialize covariance matrix
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 0, 0,
               0, 1, 0, 0,
               0, 0, 1000, 0,
               0, 0, 0, 1000;

    //Save timestamp of first measurement
    previous_timestamp_ = measurement_pack.timestamp_;

    //Done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
  *  Prediction
  ****************************************************************************/

  /*** Update the state transition matrix F according to the new elapsed time.
       - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  //compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  //update state transition matrix F
  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;

  //set the acceleration noise components
  float noise_ax = 9;
  float noise_ay = 9;

  //update process covariance matrix Q
  ekf_.Q_(0,0) = dt_4 / 4 * noise_ax;
  ekf_.Q_(0,2) = dt_3 / 2 * noise_ax;
  ekf_.Q_(1,1) = dt_4 / 4 * noise_ay;
  ekf_.Q_(1,3) = dt_3 / 2 * noise_ay;
  ekf_.Q_(2,0) = dt_3 / 2 * noise_ax;
  ekf_.Q_(2,2) = dt_2 / 1 * noise_ax;
  ekf_.Q_(3,1) = dt_3 / 2 * noise_ay;
  ekf_.Q_(3,3) = dt_2 / 1 * noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /*** Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  }

  else {
    // Laser updates
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  //cout << endl << "dt = " << dt << endl;
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
