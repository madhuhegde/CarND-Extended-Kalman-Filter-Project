#include "kalman_filter.h"
#define PI 3.14159265

#include <iostream>
using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(MatrixXd &H_in,
                        MatrixXd &R_in)
{
  H_ = H_in;
  R_ = R_in;
  
}

void KalmanFilter::Predict() 
{
  // predict the state
  x_ = F_*x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_*P_*Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) 
{
  // Update the state using Kalman Filter equations
    
  VectorXd y = z - H_*x_;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_*P_*Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_*Ht*Si;

  // New state
 
  x_ = x_ + ( K*y );
  P_ = ( I_ - K*H_ )*P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) 
{
  float px = x_[0];
  float py = x_[1];
  float vx = x_[2];
  float vy = x_[3];

  // If rho == 0, skip the update step to avoid dividing by zero.
  if( px == 0. && py == 0. )
    return;

  H_ = tools.CalculateJacobian( x_ );
    
  VectorXd hofx(3);
  float rho = sqrt( px*px + py*py );
    
  hofx << rho, atan2( py, px ), ( px*vx + py*vy )/rho;
    
  /* rho_dot check */
  if(fabs(rho) < 0.0001)
  {
      hofx[2] = 0;
  }
  // Update the state using Extended Kalman Filter equations
  VectorXd y = z - hofx;
    
  /* keep theta within -pi to + pi */
  if( y[1] > PI )
    y[1] -= 2.f*PI;
  if( y[1] < -PI )
    y[1] += 2.f*PI;
    
  
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_*P_*Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_*Ht*Si;

  // Compute new state
  x_ = x_ + ( K*y );
  P_ = ( I_ - K*H_ )*P_;
}
