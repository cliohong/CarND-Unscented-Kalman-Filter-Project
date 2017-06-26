#include "ukf.h"
#include <iostream>
using namespace::std;
using Eigen::VectorXd;
using Eigen::MatrixXd;


///**
// * Initializes Unscented Kalman filter
// */
UKF::UKF() {

  is_initialized_=false;
  // if this is false, laser measurements will be ignored(excezduring init)
  use_laser_ = false;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  //initial augmented state vector
  x_aug=VectorXd(7);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;


  //time when the state is true
  time_us_=0;
  
  //state dimension
  n_x_=5;
  
  //augmentedd state dimension
  n_aug_=7;
  
  //sigma point spreading parameter
  lambda_=3-n_aug_;
//  lambda_=0.0001*0.0001*(3)-n_aug_;
  
  //create sigma point matrix
  Xsig_pred_=MatrixXd(n_x_, 2*n_aug_+1);
  
  //weights of sigma points
  weights_=VectorXd(2*n_aug_+1);
  weights_(0)=lambda_/(lambda_+n_aug_);
  for(int i=1;i<2*n_aug_+1;++i){
    weights_(i)=0.5/(lambda_+n_aug_);
  }

  NIS_radar_=0;
  
  ///* the current NIS for laser
  NIS_laser_=0;

  
  // Process noise standard deviation longitudinal acceleration in m/s^2
//  std_a_ = 3;
//
//  // Process noise standard deviation yaw acceleration in rad/s^2
//  std_yawdd_ = 0.5;
//
//  // Laser measurement noise standard deviation position1 in m
////  std_laspx_ = 0.0015;
//  std_laspx_ = 0.125;
//  // Laser measurement noise standard deviation position2 in m
//  std_laspy_ = 0.125;
//
//  // Radar measurement noise standard deviation radius in m
//  std_radr_ = 0.003;
//
//  // Radar measurement noise standard deviation angle in rad
//  std_radphi_ = 0.03;
//
//  // Radar measurement noise standard deviation radius change in m/s
   std_a_ = 0.1;
  
  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.8;
  
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;
  
  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;
  
  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;
  
  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;
  
  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  
  
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << std_laspx_*std_laspx_, 0,
  0, std_laspy_*std_laspy_;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements
   
  */
  if(!is_initialized_){
//    cout<<"UKF: "<<endl;
    double p_x=0;
    double p_y=0;
    
    if(meas_package.sensor_type_==MeasurementPackage::RADAR){
      double rho=meas_package.raw_measurements_[0];
      double phi=meas_package.raw_measurements_[1];

      p_x=rho*cos(phi);
      p_y=rho*sin(phi);
      
    }else if (meas_package.sensor_type_==MeasurementPackage::LASER){
      p_x=meas_package.raw_measurements_[0];
      p_y=meas_package.raw_measurements_[1];
    }
    
    x_<<p_x,p_y,0,0,0;
    time_us_=meas_package.timestamp_;
    is_initialized_=true;
    return;
  }
  
  ///* Prediction
  double dt=(meas_package.timestamp_-time_us_)/1000000.0;
  time_us_=meas_package.timestamp_;
  

  Prediction(dt);
  
//  ///* Update
  
    if(meas_package.sensor_type_==MeasurementPackage::LASER){
      UpdateLidar(meas_package);
    }
  
    if(meas_package.sensor_type_==MeasurementPackage::RADAR){
    UpdateRadar(meas_package);
  }
 

  
//  cout << "x_ = " << x_ << endl;
//  cout << "P_ = " << P_ << endl;
//  cout<<"squred error L :  "<<NIS_laser_<<endl;
//  cout<<"squred error R :  "<<NIS_radar_<<endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  Xsig_aug=MatrixXd(n_aug_,2*n_aug_+1);
  x_aug.head(5)=x_;
  x_aug(5)=0;
  x_aug(6)=0;
  
  P_aug=MatrixXd(7,7);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5)=P_;
  P_aug(5,5)=std_a_*std_a_;
  P_aug(6,6)=std_yawdd_*std_yawdd_;
  //compute squre matrix
  MatrixXd L=P_aug.llt().matrixL();
  
  
  Xsig_aug.col(0)=x_aug;
  for(int i=0;i<n_aug_;++i){
    Xsig_aug.col(i+1)        =x_aug+sqrt(lambda_+n_aug_)*L.col(i);
    Xsig_aug.col(i+n_aug_+1) =x_aug-sqrt(lambda_+n_aug_)*L.col(i);
  }
  
  //predict sigma points
  for(int i=0;i<2*n_aug_+1;++i){
    double p_x=Xsig_aug(0,i);
    double p_y=Xsig_aug(1,i);
    double v=Xsig_aug(2,i);
    double yaw=Xsig_aug(3,i);
    
    double yawd=Xsig_aug(4,i);
    double nu_a=Xsig_aug(5,i);
    double nu_yawdd=Xsig_aug(6,i);
    //predicted state values
    double px_p;
    double py_p;
    if(fabs(yawd)>0.001){
      px_p=p_x+v/yawd*(sin(yaw+yawd*delta_t)-sin(yaw));
      py_p=p_y+v/yawd*(cos(yaw)-cos(yaw+yawd*delta_t));
    }else{
      px_p=p_x+v*cos(yaw)*delta_t;
      py_p=p_y+v*sin(yaw)*delta_t;
    }
    double v_p=v;
    double yaw_p=yaw+yawd*delta_t;
    
    double yawd_p=yawd;
    
    //add noise
    double tt=delta_t*delta_t;
    px_p=px_p+0.5*tt*cos(yaw)*nu_a;
    py_p=py_p+0.5*tt*sin(yaw)*nu_a;
    v_p=v_p+delta_t*nu_a;
    yaw_p=yaw_p+0.5*tt*nu_yawdd;
    
    yaw_p=angleNormalise(yaw_p);
    
    yawd_p=yawd_p+delta_t*nu_yawdd;
    
    Xsig_pred_(0,i)=px_p;
    Xsig_pred_(1,i)=py_p;
    Xsig_pred_(2,i)=v_p;
    Xsig_pred_(3,i)=yaw_p;
    Xsig_pred_(4,i)=yawd_p;
  }
//
//  //predict state mean and covariance

  x_.fill(0.0);
  for(int i=0;i<2*n_aug_+1;++i){
    x_=x_+weights_(i)*Xsig_pred_.col(i);
  }
  
  P_.fill(0.0);
  for(int i=0;i<2*n_aug_+1;++i){
    VectorXd x_diff=Xsig_pred_.col(i)-x_;
    x_diff(3)=angleNormalise(x_diff(3));
//    while(x_diff(3)>M_PI) x_diff(3)-=2.*M_PI;
//    while(x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
    P_=P_+weights_(i)*x_diff*x_diff.transpose();
      }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:
   
  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  VectorXd z = meas_package.raw_measurements_;
//  VectorXd z(2);
//  z<<meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];
  
  Zsig=MatrixXd(2,2*n_aug_+1);
  for(int i=0;i<2*n_aug_+1;++i){
    Zsig(0,i)=Xsig_pred_(0,i);
    Zsig(1,i)=Xsig_pred_(1,i);
  }
  
  VectorXd z_pred(2);
  z_pred.fill(0.0);
  for(int i=0;i<2*n_aug_+1;++i){
    z_pred=z_pred+weights_(i)*Zsig.col(i);
  }
  
  MatrixXd S=MatrixXd(2,2);
  S.fill(0.0);
  for(int i=0;i<2*n_aug_+1;++i){
    VectorXd z_diff=Zsig.col(i)-z_pred;
    S=S+weights_(i)*z_diff*z_diff.transpose();
  }
  
  //add R noise
  MatrixXd R(2,2);
  R<<std_laspx_*std_laspx_,0,
     0,   std_laspy_*std_laspy_;
  S=S+R;
  
  MatrixXd T=MatrixXd(5,2);
  T.fill(0);
  for(int i=0;i<2*n_aug_+1;++i){
    VectorXd z_diff=Zsig.col(i)-z_pred;
    VectorXd x_diff=Xsig_pred_.col(i)-x_;
    
    x_diff(3)=angleNormalise(x_diff(3));
//    while(x_diff(3)>M_PI)x_diff(3)-=2.*M_PI;
//    while(x_diff(3)<-M_PI)x_diff(3)+=2.*M_PI;
    T=T+weights_(i)*x_diff*z_diff.transpose();
  }
  
  
  MatrixXd K=T*S.inverse();
  VectorXd z_diff=z-z_pred;
  
  x_=x_+K*z_diff;
  P_=P_-K*S*K.transpose();
 
////
  //Lidar NIS
  NIS_laser_=z_diff.transpose()*S.inverse()*z_diff;
  
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  VectorXd z=meas_package.raw_measurements_;
//  cout<<"z: "<<z<<endl;
//  cout<<z<<endl;
//  VectorXd z(3);
//  z<<meas_package.raw_measurements_[0],meas_package.raw_measurements_[1],meas_package.raw_measurements_[2];
  
  Zsig=MatrixXd(3,2*n_aug_+1);
  for(int i=0;i<2*n_aug_+1;++i){
    double px=Xsig_pred_(0,i);
    double py=Xsig_pred_(1,i);
    double v=Xsig_pred_(2,i);
    double yaw=Xsig_pred_(3,i);
    double v1=cos(yaw)*v;
    double v2=sin(yaw)*v;
    
    double rho = sqrt(px*px+py*py);
    double phi = atan2(py, px);
    
    phi=angleNormalise(phi);
    
    double rho_dot = (px*v1+py*v2) / rho;
    //in the condition of zero division issue
    if(fabs(rho)<0.001){
      rho=0;
      phi=0;
      rho_dot=0;
    }

    Zsig(0,i)=rho;
    Zsig(1,i)=phi;
    Zsig(2,i)=rho_dot;
  }

  VectorXd z_pred=VectorXd(3);
  z_pred.fill(0.0);

  for(int i=0;i<2*n_aug_+1;++i){
    z_pred=z_pred+weights_(i)*Zsig.col(i);
  }
  
  //predicted measurement covarince
  MatrixXd S=MatrixXd(3,3);
  S.fill(0.0);
  for(int i=0;i<2*n_aug_+1;++i){
    VectorXd z_diff=Zsig.col(i)-z_pred;
    z_diff(1)=angleNormalise(z_diff(1));
    
//    while(z_diff(1)>M_PI) z_diff(1)-=2.*M_PI;
//    while(z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    S=S+weights_(i)*z_diff*z_diff.transpose();
  }
  //add R noise
  MatrixXd R=MatrixXd(3,3);
  R<<std_radr_*std_radr_,0,0,
  0,  std_radphi_*std_radphi_, 0,
  0,    0,          std_radrd_*std_radrd_;
  
  S=S+R;

  MatrixXd T=MatrixXd(5,3);
  T.fill(0.0);
  for(int i=0;i<2*n_aug_+1;++i){
    
    VectorXd x_dff=Xsig_pred_.col(i)-x_;
    x_dff(3)=angleNormalise(x_dff(3));
//    while(x_dff(3)>M_PI)x_dff(3)-=2.*M_PI;
//    while(x_dff(3)<-M_PI)x_dff(3)+=2.*M_PI;
//    
    VectorXd z_dff=Zsig.col(i)-z_pred;
    z_dff(1)=angleNormalise(z_dff(1));
//    while(z_dff(1)>M_PI)z_dff(1)-=2.*M_PI;
//    while(z_dff(1)<-M_PI)z_dff(1)+=2.*M_PI;
    T=T+weights_(i)*x_dff*z_dff.transpose();
  }

  MatrixXd K=T*S.inverse();
  VectorXd z_diff=z-z_pred;
  z_diff(1)=angleNormalise(z_diff(1));
//  while(z_diff(1)>M_PI)z_diff(1)-=2.*M_PI;
//  while(z_diff(1)<-M_PI)z_diff(1)+=2.*M_PI;

  x_=x_+K*z_diff;
  P_=P_-K*S*K.transpose();
//
//  //Radar NIS
  NIS_radar_=z_diff.transpose()*S.inverse()*z_diff;
}

double UKF::angleNormalise(double phi){
  const double max=M_PI;
  const double min=-M_PI;
  
  return phi<min ? max+fmod(phi-min,2.*M_PI):fmod((phi-min),2.*M_PI)+min;
}
