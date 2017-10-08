#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}


/**
* Takes in ground truth and estimation state (px,py,vx,vy) and calculates root mean squred error
*/
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
                    
  VectorXd rmse(4);
	VectorXd err(4);
	rmse << 0,0,0,0;

  /**
   * check the validity of the following inputs:
    * the estimation vector size should not be zero
    * the estimation vector size should equal ground truth vector size
    */
	if (estimations.size() != 0 and estimations.size() == ground_truth.size()) {
    	for(int i=0; i < estimations.size(); ++i){
            err = estimations[i] - ground_truth[i];
            err = err.array()*err.array();
            //accumulate squared residuals
            rmse = rmse + err; 
    	}
    
    	//calculate the mean
    	rmse = rmse/estimations.size();
    
    	//calculate the squared root
    	rmse = rmse.array().sqrt();	
	}
	else {
	    cout << "Mismatch size or empty" << endl;
	}

	//return the result
	return rmse;
}

void Tools::NormalizeAngle(double& phi)
{
  phi = atan2(sin(phi), cos(phi));
}