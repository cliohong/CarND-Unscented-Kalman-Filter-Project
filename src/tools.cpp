#include <iostream>
#include "tools.h"

namespace tools{

VectorXd CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   TODO:
   * Calculate the RMSE here.
   */
  VectorXd rmse(4);
  rmse<<0,0,0,0;
  //set the conditions
  if(estimations.size()==0||estimations.size()!=ground_truth.size()){
    cout<<"Invalid estimation or ground_truth data"<<endl;
    return rmse;
  }
  
  for (unsigned int i=0;i<estimations.size();++i){
    VectorXd residual=estimations[i]-ground_truth[i];
    residual=residual.array()*residual.array();
    rmse+=residual;
  }
  rmse=rmse/estimations.size();
  rmse=rmse.array().sqrt();
  return rmse;
  
}
}
