#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() 
{
  mse  =   VectorXd(4);
  alpha =  VectorXd(4);
  beta  =  VectorXd(4);
  rmse  =   VectorXd(4);
  
  mse << 0, 0, 0, 0;
  alpha << 0.02, 0.02, 0.01, 0.01;
  
  /* beta = (1-alpha) */
  beta << 0.98, 0.98, 0.99, 0.99;
  
}

Tools::~Tools() {}



// Calculate the RMSE
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) 

{
  int index = estimations.size();    
  

  VectorXd res = estimations[index-1] - ground_truth[index-1];
  
  res = res.array() * res.array();
  res = res.array() * alpha.array();
  
  mse = mse.array() * beta.array();
  mse = res + mse;
     

  // Calculate the RMSE
  rmse = mse.array().sqrt();

  if( rmse(0) > .12 ||
      rmse(1) > .12 ||
      rmse(2) > .52 ||
      rmse(3) > .52 )
    cout << "Warning  " << ":  rmse = " 
         << rmse(0) << "  " << rmse(1) << "  " 
         << rmse(2) << "  " << rmse(3) << endl;
    
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) 
{
  // Calculate a Jacobian for the transformation from the state vector 
  // px, py, vx, vy to the radar measurement space
  // rho, phi, rhodot.
  

  MatrixXd Hj(3,4);

  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);


  //pre-compute a set of terms to avoid repeated calculation
   float c1 = px*px+py*py;
   float c2 = sqrt(c1);
   float c3 = (c1*c2);
   
   //check division by zero
	if(fabs(c1) < 0.0001){
		cout << "CalculateJacobian () - Error - Division by Zero" << endl;
		return Hj;
	}
	
  //compute the Jacobian matrix
	Hj << (px/c2),               (py/c2),         0,              0,
		  -(py/c1),              (px/c1),         0,              0,
		   py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;
  return Hj;
}
