/**
 * @file msckf.cpp
 * @author Byeongpil Choi (bpc1224@snu.ac.kr)
 * @brief Multi-state constraint Kalman filter with point and line features class
 * @date 2025-01-07
 *
 * @copyright Copyright (c) 2025 Byeongpil Choi
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include "msckf.h"

ros::Publisher pub_state;
ros::Publisher pub_path;
ros::Publisher pub_points;
nav_msgs::Path path;

// static member variables
std::map<int, double> MSCKF::chi_squared_test_table_;

// from smsckf - Sun et al.
using namespace msckf_vio;

MSCKF::MSCKF()
{
  DXI_ = 15; // dim of IMU-state
  DXA_ = 1;  // dim of calibration state
  DXS_ = 6;  // dim of SLW-state
}

M3d MSCKF::getK() {
    M3d K;
    K << cam_intrinsics_[0],      0,        cam_intrinsics_[2],
            0,         cam_intrinsics_[1],     cam_intrinsics_[3],
            0,                0,               1;
    
    return K;
}

void MSCKF::LoadParameters(ros::NodeHandle &n){

  n.getParam("/mono_tracker/cam0/intrinsics", cam_intrinsics_);

  // MSCKF params
  n.param<int>("/estimator_node/MinSLW", MIN_SLW_, 3);
  n.param<int>("/estimator_node/MaxSLW", MAX_SLW_, 5);

  // P0
  double p_att, p_pos, p_vel, p_ba, p_bg;
  n.param<double>("/estimator_node/P0/attitude", p_att, 1e-8);
  n.param<double>("/estimator_node/P0/position", p_pos, 1e-10);
  n.param<double>("/estimator_node/P0/velocity", p_vel, 1e-5);
  n.param<double>("/estimator_node/P0/ba", p_ba, 2e-1);
  n.param<double>("/estimator_node/P0/bg", p_bg, 1e-2);

  // Q - process noise
  double q_att, q_vel, q_ba, q_bg;
  n.param<double>("/estimator_node/Q/attitude", q_att, 1e-3);
  n.param<double>("/estimator_node/Q/velocity", q_vel, 2e-3);
  n.param<double>("/estimator_node/Q/ba", q_ba, 2e-4);
  n.param<double>("/estimator_node/Q/bg", q_bg, 7e-6);
  n.param<int>("/estimator_node/num_init_samples", NUM_INIT_SAMPLES_, 500);

  // R - meas. noise
  double R_f;
  double R_e;
  double R_e_L;

  double R_zupt_v_;
  double R_zupt_p_;
  double R_zupt_q_;
  n.param<double>("/estimator_node/R_f", R_f, 1.0);
  n.param<double>("/estimator_node/R_e", R_e, 1.0);
  n.param<double>("/estimator_node/R_e_L", R_e_L, 1.0);

  n.param<double>("/estimator_node/R_zupt_v_", R_zupt_v_, 0.01);
  n.param<double>("/estimator_node/R_zupt_p_", R_zupt_p_, 0.01);
  n.param<double>("/estimator_node/R_zupt_q_", R_zupt_q_, 0.0175);

  // normalize R by focal length
  // load focal length from yaml file
  std::vector<double> cam0_intrinsics(4);
  n.getParam("/estimator_node/cam0/intrinsics", cam0_intrinsics);
  // Normalize std_dev of meas. by mean of focal length
  R_f /= 0.5*(cam0_intrinsics[0] + cam0_intrinsics[1]);
  R_e /= 0.5*(cam0_intrinsics[0] + cam0_intrinsics[1]);
  //R_e_L /= 0.5*(cam0_intrinsics[0] + cam0_intrinsics[1]);

  // Load cam-imu time-offset
  n.param<double>("/estimator_node/cam0/timeshift_cam_imu", td_, 0.0);
  
  // Chi-squared test confidence
  double chi;
  n.param<double>("/estimator_node/chi", chi, 0.95);

  // Feature triangulation parameters
  n.param<double>("estimator_node/MaxDist", MAX_DIST_, 40.0);
  n.param<double>("estimator_node/MinDist", MIN_DIST_, 1.0);
  n.param<double>("estimator_node/CostProp", COST_PROP_, 5.0);

  // Sensor option
  n.param<bool>("estimator_node/use_vision", use_vision_, true);
  n.param<bool>("estimator_node/use_gps", use_gps_, false);

  // R^{vehicle}_{IMU}
  double roll_iv, pitch_iv, yaw_iv;
  n.param<double>("estimator_node/roll_imu_vehicle", roll_iv, 0.0);
  n.param<double>("estimator_node/pitch_imu_vehicle", pitch_iv, 0.0);
  n.param<double>("estimator_node/yaw_imu_vehicle", yaw_iv, 0.0);
  V3d euler_iv(D2R*roll_iv, D2R*pitch_iv, D2R*yaw_iv);
  Rvi_ = nesl::euler2dcm(euler_iv);

  // Print
  ROS_INFO("MSCKF parameters.MinSLW: %d", MIN_SLW_);
  ROS_INFO("MSCKF parameters.MaxSLW: %d", MAX_SLW_);
  ROS_INFO("Meas. parameters.R_f (std_dev): %f [pixel/f]", R_f);
  ROS_INFO("Meas. parameters.R_e (std_dev): %f [pixel/f]", R_e);
  ROS_INFO("Meas. parameters.R_e_L (std_dev): %f [pixel/f]", R_e_L);
  ROS_INFO("Chi-squared test confidence: %f", chi);

  ROS_INFO("IMU parameters.P0.attitude (std_dev): %f", p_att);
  ROS_INFO("IMU parameters.P0.position (std_dev): %f", p_pos);
  ROS_INFO("IMU parameters.P0.velocity (std_dev): %f", p_vel);
  ROS_INFO("IMU parameters.P0.ba (std_dev): %f", p_ba);
  ROS_INFO("IMU parameters.P0.bg (std_dev): %f", p_bg);

  ROS_INFO("IMU parameters.Q.attitude (std_dev): %f", q_att);
  ROS_INFO("IMU parameters.Q.velocity (std_dev): %f", q_vel);
  ROS_INFO("IMU parameters.Q.ba (std_dev): %f", q_ba);
  ROS_INFO("IMU parameters.Q.bg(std_dev): %f", q_bg);

  ROS_INFO("Triangulation parameters.MaxDist: %f", MAX_DIST_);
  ROS_INFO("Triangulation parameters.MinDist: %f", MIN_DIST_);
  ROS_INFO("Triangulation parameters.Cost_Prop: %f", COST_PROP_);

  ROS_INFO("Predefined CAM-IMU time-offset: %f", td_);

  cov_X_ = Eigen::MatrixXd::Zero(DXI_ + DXA_, DXI_ + DXA_);
  cov_X_.block<3,3>(0,0) = pow(p_att,2)*M3d::Identity();
  if (use_gps_) {
    cov_X_(2,2) = 0.0873*0.0873;  // 5 [deg]
  }
  cov_X_.block<3,3>(3,3) = pow(p_pos,2)*M3d::Identity();
  cov_X_.block<3,3>(6,6) = pow(p_vel,2)*M3d::Identity();
  cov_X_.block<3,3>(9,9) = pow(p_ba,2)*M3d::Identity();
  cov_X_.block<3,3>(12,12) = pow(p_bg,2)*M3d::Identity();
  cov_X_(15, 15) = pow(0.01, 2);  // 10 ms

  Q_ = M12d::Zero();
  Q_.block<3,3>(0,0) = pow(q_vel,2)*M3d::Identity();
  Q_.block<3,3>(3,3) = pow(q_att,2)*M3d::Identity();
  Q_.block<3,3>(6,6) = pow(q_ba,2)*M3d::Identity();
  Q_.block<3,3>(9,9) = pow(q_bg,2)*M3d::Identity();

  R_pix_f = pow(R_f,2);
  R_pix_e = pow(R_e,2);
  R_pix_e_L = pow(R_e_L,2);

  R_zupt_v = pow(R_zupt_v_,2);
  R_zupt_p = pow(R_zupt_p_,2);
  R_zupt_q = pow(R_zupt_q_,2);

  // read from yaml file
  // Cam/IMU extrinsics
  Eigen::Isometry3d T_cam_imu =
      utils::getTransformEigen(n, "/estimator_node/cam0/T_cam_imu");
  M4d Tcb = M4d::Zero();
  Tcb.block<3,3>(0,0) = T_cam_imu.linear() * Rvi_.transpose();;
  Tcb.block<3,1>(0,3) = T_cam_imu.translation();
  Tcb(3,3) = 1.0;
  Tbc_ = Tcb.inverse();

  // Cam left-right extrinsic
  Eigen::Isometry3d T_cam_rl =
      utils::getTransformEigen(n, "/estimator_node/cam1/T_cn_cnm1");
  M4d T_rl = M4d::Zero();
  T_rl.block<3,3>(0,0) = T_cam_rl.linear();
  T_rl.block<3,1>(0,3) = T_cam_rl.translation();
  T_rl(3,3) = 1.0;
  Trl_ = T_rl;
  Tlr_ = T_rl.inverse();

  // Initialize chi squared test table with
  // given confidence <from SMSCKF, Sun et al.>
  for (int i=1; i<100; i++){
    boost::math::chi_squared chi_squared_dist(i);
    double temp_chi = 1.0 - chi;
    chi_squared_test_table_[i] =
        boost::math::quantile(chi_squared_dist, chi);
  }

}

void MSCKF::InitializeState(
  std::vector<sensor_msgs::ImuConstPtr> &imu_buf,
  std::vector<sensor_msgs::NavSatFixConstPtr>& gps_buf,
  std::vector<sensor_msgs::MagneticFieldConstPtr>& mag_buf) {

  // Coarse alignment with IMU window
  V3d f_mean = V3d::Zero();
  V3d w_mean = V3d::Zero();
  for (unsigned int i = 0; i < imu_buf.size(); i++) {
    V3d f_i(imu_buf[i]->linear_acceleration.x,
      imu_buf[i]->linear_acceleration.y,
      imu_buf[i]->linear_acceleration.z);
    V3d w_i(imu_buf[i]->angular_velocity.x,
      imu_buf[i]->angular_velocity.y,
      imu_buf[i]->angular_velocity.z);
    f_mean += f_i;
    w_mean += w_i;
  }
  f_mean /= imu_buf.size();
  w_mean /= imu_buf.size();

  double roll0 = std::atan2(-f_mean(1), -f_mean(2));
  double pitch0 = std::atan(f_mean(0) /
    std::sqrt(f_mean(1)*f_mean(1) + f_mean(2)*f_mean(2)));

  // Heading initialization
  // TODO: heading coordinate frame should be checked !
  V3d mag_mean = V3d::Zero();
  for (unsigned int i = 0; i < mag_buf.size(); i++) {
    V3d mag_i(mag_buf[i]->magnetic_field.x,
      mag_buf[i]->magnetic_field.y,
      mag_buf[i]->magnetic_field.z);
    mag_mean += mag_i;
  }
  mag_mean /= mag_buf.size();
  double heading0 = 0;
  if (!mag_buf.empty() && use_gps_) {
    // Magnetic North to Geodetic North
    // @ SNU, 2021. 10. 01
    heading0 = -atan2(mag_mean(1), mag_mean(0)) + D2R*13.7;
  }

  Pose init_pose;
  V3d euler0(roll0, pitch0, heading0);
  init_pose.p << 0.0, 0.0, 0.0;
  init_pose.q = nesl::euler2quat(euler0);
  init_pose.R = nesl::quat2dcm(init_pose.q);
  Xi_.T = init_pose;
  Xi_.v << 0.0, 0.0, 0.0;
  Xi_.ba << 0.0, 0.0, 0.0;
  Xi_.bg = w_mean;  // when stationary at begining
  Xi_.time_stamp = 0.0;
  grv_ << 0.0, 0.0, 9.81;  // referenced at ned frame

  // Initialize the first ecef position
  V3d dllh_mean = V3d::Zero();
  for (unsigned int i = 0; i < gps_buf.size(); i++) {
    V3d llh_i(gps_buf[i]->latitude,
      gps_buf[i]->longitude, gps_buf[i]->altitude);
    dllh_mean += llh_i;
  }
  dllh_mean /= gps_buf.size();
  V3d llh_mean(D2R*dllh_mean(0), D2R*dllh_mean(1), dllh_mean(2));
  if (!gps_buf.empty() && use_gps_) {
    ecef0_ = NAV::llh2xyz(llh_mean);
  }
  else {
    ecef0_ = V3d::Zero();
  }

  // Print the results
  std::cout << "Mean fm: " << f_mean(0) << ", " << f_mean(1) << ", " 
    << f_mean(2) << std::endl;

  std::cout << "Initial Llh [deg, deg, m]: " << dllh_mean(0) << ", "
    << dllh_mean(1) << ", " << dllh_mean(2) << " with "
    << gps_buf.size() << std::endl;

  std::cout << "Initial roll, pith, heading [deg]: " << R2D*roll0 << ", " 
    << R2D*pitch0 << ", " << R2D*heading0 << std::endl;

  std::cout << "Initial gyro bias [rad/s]: " << w_mean(0) << ", "
    << w_mean(1) << ", " << w_mean(2) << std::endl;

  ROS_INFO("State initialization completed !");
}

void MSCKF::propagate(const sensor_msgs::ImuConstPtr &imu_msg,
                      double dt){
  // Nominal propagation - Euler integration
  // retrieve IMU measurements
  double ax = imu_msg->linear_acceleration.x;
  double ay = imu_msg->linear_acceleration.y;
  double az = imu_msg->linear_acceleration.z;
  V3d fib_b;
  fib_b << ax, ay, az;
  fib_b -=  Xi_.ba;

  double wx = imu_msg->angular_velocity.x;
  double wy = imu_msg->angular_velocity.y;
  double wz = imu_msg->angular_velocity.z;
  V3d wib_b;
  wib_b << wx, wy, wz;
  wib_b -= Xi_.bg;

  // IMU-state propagation
  if (dt > 0){
      NavState Xi_new =
          propagateIMU_Euler(fib_b, wib_b, dt);
      Eigen::MatrixXd cov_new =
          propagateCovX(fib_b, dt);
      Xi_ = Xi_new;
      cov_X_ = cov_new;
  }

}

Eigen::MatrixXd MSCKF::propagateCovX(V3d fib_b, double dt){
  M15d F = M15d::Zero();
  F.block<3,3>(0,12) = -Xi_.T.R;
  V3d Rfx;
  Rfx = Xi_.T.R*fib_b;
  F.block<3,3>(6,0) = -nesl::Vec2SkewMat(Rfx);
  F.block<3,3>(6,9) = -Xi_.T.R;
  F.block<3,3>(3,6) = M3d::Identity();

  Eigen::Matrix<double, 15, 12> G =
      Eigen::Matrix<double, 15, 12>::Zero();
  G.block<3,3>(0,3) = -Xi_.T.R;
  G.block<3,3>(6,0) = -Xi_.T.R;
  G.block<3,3>(9,6) = M3d::Identity();
  G.block<3,3>(12,9) = M3d::Identity();

  M15d Phi;
  Phi = M15d::Identity() + F*dt + 0.5*F*F*pow(dt,2);

  M15d cov_Xi_pred;
  M15d cov_Xi = cov_X_.block<15,15>(0,0);

  cov_Xi_pred = Phi*cov_Xi*Phi.transpose()
      + (Phi*G*Q_*G.transpose()*Phi.transpose()
         + G*Q_*G.transpose())*dt*0.5;

  Eigen::MatrixXd cov_pred = cov_X_;
  cov_pred.block<15, 15>(0, 0) = cov_Xi_pred;
  cov_pred(15, 15) += 1e-8;  // t_d random walk
  if (cov_X_.rows() == DXI_ + DXA_){
    return cov_pred;
  }
  else{
    // cov mtx propagation with SLW
    Eigen::MatrixXd cov_IS;
    Eigen::MatrixXd cov_diag;
    Eigen::MatrixXd cov_SLW;

    int Nslw = Xs_.size();
    cov_SLW = cov_X_.bottomRightCorner(6*Nslw+DXA_, 6*Nslw+DXA_);
    cov_IS = cov_X_.topRightCorner(DXI_, 6*Nslw+DXA_);
    cov_diag = Phi*cov_IS;

    // Initialize size of cov_new
    Eigen::MatrixXd cov_new(DXI_ + DXA_ + 6*Nslw, 
      DXI_ + DXA_ + 6*Nslw);
    cov_new << cov_Xi_pred, cov_diag,
        cov_diag.transpose(), cov_SLW;
    return cov_new;
  }

}

NavState MSCKF::propagateIMU_Euler
  (V3d fib_b, V3d wib_b, double dt){
  // Attitude update
  V3d theta = wib_b*dt;
  V4d del_quat = nesl::rvec2quat(theta);
  V4d q_gb = nesl::quatMultiply(Xi_.T.q, del_quat);

  // Acceleration transform
  M3d Cgb1 = nesl::quat2dcm(Xi_.T.q);
  M3d Cgb2 = nesl::quat2dcm(q_gb);
  M3d C12 = Cgb1.transpose()*Cgb2;

  V3d s = 0.5*(fib_b + C12*fib_b)*dt;
  V3d y = 0.5*s*dt;
  V3d g = grv_;

  // State update
  V3d p = Xi_.T.p + Xi_.v*dt + 0.5*g*pow(dt,2) + Cgb1*y;
  V3d v = Xi_.v + Cgb1*s + g*dt;
  V4d q = q_gb;
  V3d ba = Xi_.ba;
  V3d bg = Xi_.bg;

  Pose pos;
  pos.p = p;
  pos.q = q;
  pos.R = nesl::quat2dcm(q);
  pos.v = v;
  NavState state_updated;
  state_updated.T = pos;
  state_updated.v = v;
  state_updated.ba = ba;
  state_updated.bg = bg;

  return state_updated;
}

void MSCKF::correct_position(
  const sensor_msgs::NavSatFixConstPtr& gps_msg,
  const int& frame_id) {

  // Obtain measurements
  V3d llh_i(D2R*gps_msg->latitude, D2R*gps_msg->longitude,
    gps_msg->altitude);

  // Transform geodetic to local coordinate
  V3d ecef_i = NAV::llh2xyz(llh_i);
  V3d ned_i = NAV::xyz2ned(ecef_i, ecef0_);
  double std_p = 2;  // [m]
  double std_h = 8;  // [m]

  // filter update implementation
  int Nstate = cov_X_.cols();
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(3,Nstate);
  H.block<3,3>(0,3) = M3d::Identity();
  M3d R = M3d::Identity();
  R(0,0) = std_p*std_p;
  R(1,1) = std_p*std_p;
  R(2,2) = std_h*std_h;
  Eigen::MatrixXd S = H*cov_X_*H.transpose()+R;
  Eigen::MatrixXd K_t = S.ldlt().solve(H*cov_X_);
  Eigen::MatrixXd K = K_t.transpose();
  Eigen::MatrixXd I_KH =
      Eigen::MatrixXd::Identity(K.rows(),H.cols()) - K*H;
  Eigen::MatrixXd cov_new = I_KH*cov_X_*I_KH.transpose() +
      K*R*K.transpose();
  V3d r = ned_i - Xi_.T.p;
  Eigen::VectorXd delx_hat = K*r;

  // filter update
  cov_X_ = cov_new;

  // IMU/SLW-state update
  UpdateIMUstate(delx_hat.segment<15>(0));
  td_ += delx_hat(15);
  UpdateSLWstate(
    delx_hat.segment(DXI_+DXA_, DXS_*Xs_.size()), frame_id);
  ROS_INFO("GPS position update with NED position meas [%.2f  %.2f  %.2f]",
    ned_i(0), ned_i(1), ned_i(2));
}

void MSCKF::AugmentSLW(int frame_id, const sensor_msgs::ImuConstPtr &imu_msg){
  int Nslw = Xs_.size();
  int Nmsckf = DXI_ + DXA_ + DXS_*Nslw;

  // Extract IMU meas.
  V3d w_hat;
  double wx = imu_msg->angular_velocity.x;
  double wy = imu_msg->angular_velocity.y;
  double wz = imu_msg->angular_velocity.z;
  w_hat << wx, wy, wz;
  w_hat -= Xi_.bg;

  // Augmented state Jacobian
  Eigen::MatrixXd J_aug;
  J_aug = Eigen::MatrixXd::Zero(DXS_, Nmsckf);
  // SLW-pose Jacobian block
  J_aug.block<3,3>(0,0) = M3d::Identity();
  J_aug.block<3,3>(3,3) = M3d::Identity();
  // Calibration Jacobian block
  J_aug.block<3,1>(0,15) = Xi_.T.R*w_hat;
  J_aug.block<3,1>(3,15) = Xi_.v; 

  Eigen::MatrixXd cov_Aug(Nmsckf+DXS_, Nmsckf+DXS_);
  cov_Aug << cov_X_, cov_X_*J_aug.transpose(),
          J_aug*cov_X_, J_aug*cov_X_*J_aug.transpose();

  // Calculate augmented cov mtx.
  // and state augmentation.
  // if MaxSLW is not reached
  // only occurs at beginning
  if (Nslw < MAX_SLW_){
    Xs_.emplace(std::make_pair(frame_id, Xi_.T));
    cov_X_ = cov_Aug;
  }
  // when SLW is full (FIFO)
  else {
    Xs_.erase(frame_id - MAX_SLW_);
    Xs_.emplace(std::make_pair(frame_id, Xi_.T));
    Eigen::MatrixXd reduced_cov;
    reduced_cov = Eigen::MatrixXd::Zero(Nmsckf, Nmsckf);
    reduced_cov << cov_Aug.topLeftCorner(DXI_ + DXA_, DXI_ + DXA_),
            cov_Aug.topRightCorner(DXI_ + DXA_, DXS_*Nslw),
            cov_Aug.bottomLeftCorner(DXS_*Nslw, DXI_ + DXA_),
            cov_Aug.bottomRightCorner(DXS_*Nslw, DXS_*Nslw);
    cov_X_ = reduced_cov;
  }

}


Eigen::MatrixXd MSCKF::validateCovMtx
(Eigen::MatrixXd P_mtx){

  Eigen::MatrixXd P_out;
  // force to make symmetric
  P_out = 0.5 * (P_mtx + P_mtx.transpose());
  Eigen::LLT<Eigen::MatrixXd> lltOfP(P_out);

  if (lltOfP.info() == Eigen::NumericalIssue){
    ROS_WARN("System Cov Mtx is not Positive Definite !");
    Eigen::SelfAdjointEigenSolver<MatrixXd> sv(P_out);
    Eigen::MatrixXd D = sv.eigenvalues().asDiagonal();
    Eigen::MatrixXd V = sv.eigenvectors();

    for (int i = 0; i < D.rows(); i++){
      if (D(i,i) < 2.2204e-16) D(i,i) = 1e-12;
    }
    P_out = V*D*V.transpose();
  }

  return P_out;
}


NavState MSCKF::getNavState(){
  return Xi_;
}

void MSCKF::registerPub(ros::NodeHandle &n){
  pub_state = n.advertise<nav_msgs::Odometry>
    ("/smsckf_nesl/odom", 1000);
  pub_path = n.advertise<nav_msgs::Path>
    ("/smsckf_nesl/path", 1000);
  //pub_points = n.advertise<sensor_msgs::PointCloud2>
    //("/smsckf_nesl/points", 1000);
}

void MSCKF::pubFilterState(const std_msgs::Header &header){
  // publish tf (imu pose)
  static tf::TransformBroadcaster tf_pub;
  geometry_msgs::TransformStamped odom_trans;
  odom_trans.header = header;
  odom_trans.header.frame_id = "world";
  odom_trans.child_frame_id = "IMU";

  odom_trans.transform.translation.x = Xi_.T.p(0,0);
  odom_trans.transform.translation.y = Xi_.T.p(1,0);
  odom_trans.transform.translation.z = Xi_.T.p(2,0);
  odom_trans.transform.rotation.w = Xi_.T.q(0,0);
  odom_trans.transform.rotation.x = Xi_.T.q(1,0);
  odom_trans.transform.rotation.y = Xi_.T.q(2,0);
  odom_trans.transform.rotation.z = Xi_.T.q(3,0);

  tf_pub.sendTransform(odom_trans);

  geometry_msgs::TransformStamped cam0_trans;
  cam0_trans.header = header;
  cam0_trans.header.frame_id = "IMU";
  cam0_trans.child_frame_id = "cam0";
  M3d R_bc = Tbc_.block<3,3>(0,0);
  V4d q_bc = nesl::dcm2quat(R_bc);
  // T^{header.frame_id}_{child_frame_id}
  cam0_trans.transform.translation.x = Tbc_(0,3);
  cam0_trans.transform.translation.y = Tbc_(1,3);
  cam0_trans.transform.translation.z = Tbc_(2,3);
  cam0_trans.transform.rotation.w = q_bc(0,0);
  cam0_trans.transform.rotation.x = q_bc(1,0);
  cam0_trans.transform.rotation.y = q_bc(2,0);
  cam0_trans.transform.rotation.z = q_bc(3,0);
  tf_pub.sendTransform(cam0_trans);

  // publish odometry
  nav_msgs::Odometry odometry;
  odometry.header = header;
  odometry.header.frame_id = "world";
  odometry.child_frame_id = "IMU";
  odometry.pose.pose.position.x = Xi_.T.p(0,0);
  odometry.pose.pose.position.y = Xi_.T.p(1,0);
  odometry.pose.pose.position.z = Xi_.T.p(2,0);
  odometry.pose.pose.orientation.w = Xi_.T.q(0,0);
  odometry.pose.pose.orientation.x = Xi_.T.q(1,0);
  odometry.pose.pose.orientation.y = Xi_.T.q(2,0);
  odometry.pose.pose.orientation.z = Xi_.T.q(3,0);

  V3d v_body = Xi_.T.R.transpose()*Xi_.v;
  odometry.twist.twist.linear.x = v_body(0,0);
  odometry.twist.twist.linear.y = v_body(1,0);
  odometry.twist.twist.linear.z = v_body(2,0);

  // odometry : covariance
//  M6d Cov_imu_pose = M6d::Zero();
//  // x = [theta p]';
//  Cov_imu_pose = cov_X_.block<6,6>(0,0);
//  for (int i = 0; i < 6; ++i)
//    for (int j = 0; j < 6; ++j)
//      odometry.pose.covariance[6*i+j] =
//          Cov_imu_pose(i, j);

  M16d Cov_imu = M16d::Zero();
  // x = [q p v b_a b_g td]';
  Cov_imu = cov_X_.block<16,16>(0,0);
  for (unsigned int i = 0; i < 16; i ++)
      odometry.pose.covariance[i] = Cov_imu(i, i);
  for (unsigned int i = 0; i < 3; i++)
      odometry.pose.covariance[16+i] = Xi_.v(i);
  for (unsigned int i = 0; i < 3; i++)
      odometry.pose.covariance[19+i] = Xi_.ba(i);
  for (unsigned int i = 0; i < 3; i++)
      odometry.pose.covariance[22+i] = Xi_.bg(i);
  odometry.pose.covariance[25] = td_;
  for (unsigned int i = 0; i < 10; i++)
      odometry.pose.covariance[26+i] = 0.0;

  pub_state.publish(odometry);

  // publish path
  geometry_msgs::PoseStamped pose_stamped;
  pose_stamped.header = header;
  pose_stamped.header.frame_id = "world";
  pose_stamped.pose = odometry.pose.pose;
  path.header = header;
  path.header.frame_id = "world";
  path.poses.push_back(pose_stamped);
  pub_path.publish(path);

  // publish msckf points
  if (msckf_points_.size() > 0) {
    sensor_msgs::PointCloud2 cloud;
    cloud.header = header;
    cloud.header.frame_id = "world";
    cloud.width = 3*msckf_points_.size();
    cloud.height = 1;
    cloud.is_bigendian = false;
    cloud.is_dense = false;

    // Setup pointcloud fields
    sensor_msgs::PointCloud2Modifier modifier(cloud);
    modifier.setPointCloud2FieldsByString(1,"xyz");
    modifier.resize(3*msckf_points_.size());

    // Iterators
    sensor_msgs::PointCloud2Iterator<float> out_x(cloud, "x");
    sensor_msgs::PointCloud2Iterator<float> out_y(cloud, "y");
    sensor_msgs::PointCloud2Iterator<float> out_z(cloud, "z");

    // Fill our iterators
    for(const auto &pt : msckf_points_) {
        *out_x = pt(0); ++out_x;
        *out_y = pt(1); ++out_y;
        *out_z = pt(2); ++out_z;
    }
    pub_points.publish(cloud);
  }

}

void MSCKF::correct
  (const std::map<int, std::vector<std::pair<int, V2d>>> &deadTracks_f, 
      const std::map<int, std::vector<std::pair<int, V2d>>> &deadTracks_e, 
        const std::map<int, std::vector<std::pair<int, V4d>>> &deadTracks_e_L, int frame_id, const sensor_msgs::ImuConstPtr &imu_msg){
  // Extract IMU meas
  V3d w_hat;
  double wx = imu_msg->angular_velocity.x;
  double wy = imu_msg->angular_velocity.y;
  double wz = imu_msg->angular_velocity.z;
  w_hat << wx, wy, wz;
  w_hat -= Xi_.bg;

  int chi_check = 0;
  int f_size, e_size, e_L_size = 0;

  // Initialize big matrix
  Eigen::MatrixXd Ho = Eigen::MatrixXd::Zero
      (2*(deadTracks_f.size()+deadTracks_e.size()+deadTracks_e_L.size())*MAX_SLW_,
       DXI_ + DXA_ + DXS_*Xs_.size());
  Eigen::VectorXd ro = Eigen::VectorXd::Zero
      (2*(deadTracks_f.size()+deadTracks_e.size()+deadTracks_e_L.size())*MAX_SLW_);

  int stack_cnt = 0;
  ros::Time start_time = ros::Time::now();
  msckf_points_.clear();

  if (!deadTracks_f.empty()) {
    for (auto &track : deadTracks_f){
//    feature_id : track.first
//    0th_view : track.second[0].first
//    0th_pt : track.second[0].second
    // two-view triangulation
      V3d f_c0;
      f_c0 = triangulateTwoViews(track);
    // multi-view triangulation
      if(triangulateTrackWithGN(track,f_c0, R_pix_f)){
      // Compute Jacobians for a deadTrack
        Eigen::MatrixXd Hx_j;
        Eigen::MatrixXd Hf_j;
        Eigen::VectorXd r_j = track_info_.residual;
      // std::cout << r_j << std::endl;
      // for stereo-update
    //  Eigen::MatrixXd R_j = R_pix_*
    //      Eigen::MatrixXd::Identity
    //      (2*track.second.size(),2*track.second.size());
    //  calculatesMeasJacobian(track, Hx_j, Hf_j, frame_id, w_hat);
      // for mono-update
        Eigen::MatrixXd R_j = R_pix_f*
            Eigen::MatrixXd::Identity
            (2*track.second.size(),2*track.second.size());
        calculateMeasJacobian(track, Hx_j, Hf_j, frame_id, w_hat);

      // Nullspace projection
        Eigen::JacobiSVD<Eigen::MatrixXd>
            svd(Hf_j, ComputeFullU | ComputeThinV);
        Eigen::MatrixXd A =
            svd.matrixU().rightCols(Hf_j.rows()-3);
     // std::cout << Hf_j << std::endl;
      //std::cout << A << std::endl;
        Eigen::MatrixXd Ho_j = A.transpose()*Hx_j;
        Eigen::VectorXd ro_j = A.transpose()*r_j;
        Eigen::MatrixXd Ro_j = A.transpose()*R_j*A;
        double gamma = ro_j.transpose()*
            (Ho_j*cov_X_*Ho_j.transpose()+Ro_j).ldlt().solve(ro_j);
        int dof = ro_j.size();

        // pass chi-squared test
        if(ChiSquaredTest(gamma, dof)){
          // stack Jacobian for update
          Ho.block(stack_cnt, 0, Ho_j.rows(), Ho_j.cols())
              = Ho_j;
          ro.segment(stack_cnt, ro_j.rows()) = ro_j;
          stack_cnt += Ho_j.rows();
          msckf_points_.push_back(track_info_.f_g);

          chi_check++;
        }
      }
    }
  }
  f_size = stack_cnt;

  if (!deadTracks_e.empty()) {
    for (auto &track : deadTracks_e){
//    feature_id : track.first
//    0th_view : track.second[0].first
//    0th_pt : track.second[0].second
    // two-view triangulation
      V3d f_c0;
      f_c0 = triangulateTwoViews(track);
    // multi-view triangulation
      if(triangulateTrackWithGN(track,f_c0, R_pix_e)){
      // Compute Jacobians for a deadTrack
        Eigen::MatrixXd Hx_j;
        Eigen::MatrixXd Hf_j;
        Eigen::VectorXd r_j = track_info_.residual;
      // std::cout << r_j << std::endl;
      // for stereo-update
    //  Eigen::MatrixXd R_j = R_pix_*
    //      Eigen::MatrixXd::Identity
    //      (2*track.second.size(),2*track.second.size());
    //  calculatesMeasJacobian(track, Hx_j, Hf_j, frame_id, w_hat);
      // for mono-update
        Eigen::MatrixXd R_j = R_pix_e*
            Eigen::MatrixXd::Identity
            (2*track.second.size(),2*track.second.size());
        calculateMeasJacobian(track, Hx_j, Hf_j, frame_id, w_hat);

      // Nullspace projection
        Eigen::JacobiSVD<Eigen::MatrixXd>
            svd(Hf_j, ComputeFullU | ComputeThinV);
        Eigen::MatrixXd A =
            svd.matrixU().rightCols(Hf_j.rows()-3);
     // std::cout << Hf_j << std::endl;
      //std::cout << A << std::endl;
        Eigen::MatrixXd Ho_j = A.transpose()*Hx_j;
        Eigen::VectorXd ro_j = A.transpose()*r_j;
        Eigen::MatrixXd Ro_j = A.transpose()*R_j*A;
        double gamma = ro_j.transpose()*
            (Ho_j*cov_X_*Ho_j.transpose()+Ro_j).ldlt().solve(ro_j);
        int dof = ro_j.size();

        // pass chi-squared test
        if(ChiSquaredTest(gamma, dof)){
          // stack Jacobian for update
          Ho.block(stack_cnt, 0, Ho_j.rows(), Ho_j.cols())
              = Ho_j;
          ro.segment(stack_cnt, ro_j.rows()) = ro_j;
          stack_cnt += Ho_j.rows();
          msckf_points_.push_back(track_info_.f_g);

          chi_check++;
        }
      }
    }
  }

  e_size = stack_cnt - f_size;

  if (!deadTracks_e_L.empty()) {
    for (auto &track : deadTracks_e_L){
//    feature_id : track.first
//    0th_view : track.second[0].first
//    0th_pt : track.second[0].second
    // two-view triangulation


      // for robust line triangulation /////////////// 2024. 11. 28
      // if (track.second.size() < 7) {
      //   continue;
      // }
      ////////////////////////////////////////////////////////////////

      V6d L_c1_0;
      L_c1_0 = triangulateTwoViews_L(track);
    // multi-view triangulation
      if(triangulateTrackWithGN_L(track, L_c1_0, R_pix_e_L)){
      // Compute Jacobians for a deadTrack
        Eigen::MatrixXd Hx_j;
        Eigen::MatrixXd Hf_j;
        Eigen::VectorXd r_j = track_info_L.residual;
      // std::cout << r_j << std::endl;
      // for stereo-update
    //  Eigen::MatrixXd R_j = R_pix_*
    //      Eigen::MatrixXd::Identity
    //      (2*track.second.size(),2*track.second.size());
    //  calculatesMeasJacobian(track, Hx_j, Hf_j, frame_id, w_hat);
      // for mono-update
        Eigen::MatrixXd R_j = R_pix_e_L*
            Eigen::MatrixXd::Identity
            (2*track.second.size(),2*track.second.size());
        calculateMeasJacobian_L(track, Hx_j, Hf_j, frame_id, w_hat);

      // Nullspace projection
        Eigen::JacobiSVD<Eigen::MatrixXd>
            svd(Hf_j, ComputeFullU | ComputeThinV);
        Eigen::MatrixXd A =
            svd.matrixU().rightCols(Hf_j.rows()-4);
     // std::cout << Hf_j << std::endl;
      //std::cout << A << std::endl;
        Eigen::MatrixXd Ho_j = A.transpose()*Hx_j;
        Eigen::VectorXd ro_j = A.transpose()*r_j;
        Eigen::MatrixXd Ro_j = A.transpose()*R_j*A;
        double gamma = ro_j.transpose()*
            (Ho_j*cov_X_*Ho_j.transpose()+Ro_j).ldlt().solve(ro_j);
        int dof = ro_j.size();
        // pass chi-squared test
        if(ChiSquaredTest(gamma, dof)){
          // stack Jacobian for update
          Ho.block(stack_cnt, 0, Ho_j.rows(), Ho_j.cols())
              = Ho_j;
          ro.segment(stack_cnt, ro_j.rows()) = ro_j;
          stack_cnt += Ho_j.rows();
          
          //ROS_WARN("chi pass!!!!");

          chi_check++;
        }
      }
    }
  }

  e_L_size = stack_cnt - f_size - e_size;

  //ROS_WARN("before chi: %zu", deadTracks_f.size() + deadTracks_e.size());
  //ROS_WARN("after chi: %d", chi_check);

  start_time = ros::Time::now();
  // Resize Jacobian and residual for less execution time
  if (stack_cnt > 400) {
     stack_cnt = 400;
     if (f_size > 400) {
      f_size = 400;
      e_size = 0;
      e_L_size = 0;
     }
     else {
      if (f_size + e_size > 400) {
        e_size = 400 - f_size;
        e_L_size = 0;
      }
      else {
        e_L_size = 400 - (f_size + e_size);
      }
     }
  }
  Ho.conservativeResize(stack_cnt, Ho.cols());
  ro.conservativeResize(stack_cnt);

  //ROS_WARN("size f, e: %d, %d", f_size, e_size);

  MeasurementUpdate(Ho, ro, frame_id, f_size, e_size, e_L_size);
  double update_processing_time =
      (ros::Time::now()-start_time).toSec();
}

void MSCKF::MeasurementUpdate(const Eigen::MatrixXd &Ho,
                              const Eigen::VectorXd &ro,
                              int frame_id, const int &f_size, const int &e_size, const int &e_L_size){
  if (Ho.rows() == 0 || ro.rows() == 0 ) return;

  Eigen::MatrixXd Th;
  Eigen::VectorXd rn;
  Eigen::MatrixXd Rn;

  Eigen::MatrixXd Ro_f = R_pix_f*Eigen::MatrixXd::Identity
      (f_size, f_size);
  Eigen::MatrixXd Ro_e = R_pix_e*Eigen::MatrixXd::Identity
      (e_size, e_size);
  Eigen::MatrixXd Ro_e_L = R_pix_e_L*Eigen::MatrixXd::Identity
      (e_L_size, e_L_size);

  // Create the combined diagonal matrix Ro
  Eigen::VectorXd combined_diagonal(f_size + e_size + e_L_size);

  // Copy the diagonal elements of Ro_f
  combined_diagonal.segment(0, f_size) = Ro_f.diagonal();

  // Copy the diagonal elements of Ro_e
  combined_diagonal.segment(f_size, e_size) = Ro_e.diagonal();

  // Copy the diagonal elements of Ro_e_L
  combined_diagonal.segment(f_size + e_size, e_L_size) = Ro_e_L.diagonal();

  // Create the final diagonal matrix
  Eigen::MatrixXd Ro = combined_diagonal.asDiagonal();

  if (Ho.rows() > Ho.cols()){
  // QR decomposition
  Eigen::HouseholderQR<Eigen::MatrixXd> qr_helper(Ho);
  Eigen::MatrixXd Q = qr_helper.householderQ();
  Eigen::MatrixXd Q1 = Q.leftCols(DXI_ + DXA_ + DXS_*Xs_.size());
  Th = Q1.transpose()*Ho;
  rn = Q1.transpose()*ro;
  Rn = Q1.transpose()*Ro*Q1;
  }
  else{
    Th = Ho;
    rn = ro;
    Rn = Ro;
  }
  //Rn = R_pix_f*Eigen::MatrixXd::Identity
  //    (Th.rows(),Th.rows());

  // Compute Kalman gain
  // const Eigen::MatrixXd &P = cov_X_;
  Eigen::MatrixXd Sn = Th*cov_X_*Th.transpose() + Rn;
  Eigen::MatrixXd Kn_transpose = Sn.ldlt().solve(Th*cov_X_);
  Eigen::MatrixXd Kn = Kn_transpose.transpose();
  Eigen::VectorXd delx_hat = Kn*rn;
  Eigen::MatrixXd I_KH = Eigen::MatrixXd::Identity
      (delx_hat.rows(), delx_hat.rows()) - Kn*Th;
  Eigen::MatrixXd P_hat = I_KH*cov_X_*I_KH.transpose() +
      Kn*Rn*Kn.transpose();
  P_hat = validateCovMtx(P_hat);
  cov_X_ = P_hat;

  // IMU/SLW-state update
  UpdateIMUstate(delx_hat.segment<15>(0));
  td_ += delx_hat(15);
  std::cout << "time offset: " << td_ << ", std_dev: " << sqrt(cov_X_(15,15)) << std::endl;
  UpdateSLWstate(delx_hat.segment(DXI_ + DXA_, DXS_*Xs_.size()), frame_id);

  ROS_INFO("Filter update using %lu measurements",
           ro.size());

}

void MSCKF::UpdateSLWstate(Eigen::VectorXd delx,
                           int frame_id){
  for (int i=0; i<Xs_.size(); i++){
    int slw_idx = frame_id - Xs_.size() + i;
    V6d del_slw = delx.segment<6>(DXS_*i);
    V3d del_r = del_slw.segment<3>(0);
    V4d del_q = nesl::rvec2quat(del_r);
    V3d del_p = del_slw.segment<3>(3);

    auto old_slw = Xs_.find(slw_idx);
    Pose old_pose = old_slw->second;
    Pose new_pose;
    new_pose.q =
        nesl::quatMultiply(del_q, old_pose.q);
    new_pose.R = nesl::quat2dcm(new_pose.q);
    new_pose.p = old_pose.p + del_p;
    new_pose.v = old_pose.v;

    Xs_.at(slw_idx) = new_pose;
  }
}

void MSCKF::UpdateIMUstate(V15d delx){
  V3d delr = delx.segment<3>(0);
  V4d delq = nesl::rvec2quat(delr);
  V3d delp = delx.segment<3>(3);
  V3d delv = delx.segment<3>(6);
  V3d delba = delx.segment<3>(9);
  V3d delbg = delx.segment<3>(12);

  Xi_.T.q = nesl::quatMultiply(delq, Xi_.T.q);
  Xi_.T.R = nesl::quat2dcm(Xi_.T.q);
  Xi_.T.p += delp;
  Xi_.v += delv;
  Xi_.ba += delba;
  Xi_.bg += delbg;
}

bool MSCKF::ChiSquaredTest(double gamma, int dof){
  if (chi_squared_test_table_[dof] > gamma){
    return true;
  }
  else{
    // ROS_WARN("Failed to pass Chi-squared test !");
    return false;
  }
}

void MSCKF::calculateMeasJacobian
  (const std::pair<int,
    std::vector<std::pair<int, V2d>>> &track,
   Eigen::MatrixXd &H_x,
   Eigen::MatrixXd &H_f,
   int frame_id, V3d w_hat){

  // Matrix size initialization
  int Ntrack = track.second.size();
  int Nslw = Xs_.size();
  H_x = Eigen::MatrixXd::Zero(2*Ntrack, DXI_ + DXA_ + DXS_*Nslw);
  H_f = Eigen::MatrixXd::Zero(2*Ntrack, 3);

  for (int i=0; i<Ntrack; i++){
    int slw_id = track.second[i].first;
    auto slw = Xs_.find(slw_id);
    Pose pose_slw = slw->second;
    M4d Tgb_i;
    Tgb_i << pose_slw.R, pose_slw.p, 0.0, 0.0, 0.0, 1.0;
    M4d Tgc_i = Tgb_i*Tbc_;
    V4d hpf_c;
    V4d hpg_c;
    hpg_c << track_info_.f_g, 1.0;
    hpf_c = Tgc_i.inverse()*hpg_c;
    V3d pf_c = hpf_c.block<3,1>(0,0);
    double X = pf_c(0);
    double Y = pf_c(1);
    double Z = pf_c(2);

    Eigen::Matrix<double,2,3> J_pi;
    J_pi << 1.0/Z, 0.0, -X/pow(Z,2),
            0.0, 1.0/Z, -Y/pow(Z,2);
    Eigen::Matrix<double,2,3> Hfi_j =
        J_pi*(Tgc_i.block<3,3>(0,0)).transpose();

    H_f.block<2,3>(2*i,0) =Hfi_j;

    // slw index btw 1~MaxSLW
    int slw_idx = Nslw-(frame_id-slw_id);
    V3d temp = track_info_.f_g - Tgb_i.block<3,1>(0,3);
    H_x.block<2,3>(2*i, DXI_ + DXA_ + DXS_*slw_idx) =
        Hfi_j*nesl::Vec2SkewMat(temp);
    H_x.block<2,3>(2*i, DXI_+DXA_+(DXS_*slw_idx)+3) = -Hfi_j;

    // Calibration-related Jacobian
    V2d H_c;
    V3d w_g = Tgb_i.block<3,3>(0,0)*w_hat;
    V3d M = -(nesl::Vec2SkewMat(w_g)*temp - pose_slw.v);
    H_c << J_pi*(Tgc_i.block<3,3>(0,0)).transpose()*M;
    H_x.block<2,1>(2*i, DXI_) = H_c;
  }
}

void MSCKF::calculatesMeasJacobian
  (const std::pair<int,
    std::vector<std::pair<int, V4d>>> &track,
   Eigen::MatrixXd &H_x,
   Eigen::MatrixXd &H_f,
   int frame_id, V3d w_hat){

  // Matrix size initialization
  int Ntrack = track.second.size();
  int Nslw = Xs_.size();
  H_x = Eigen::MatrixXd::Zero(4*Ntrack, DXI_ + DXA_ + DXS_*Nslw);
  H_f = Eigen::MatrixXd::Zero(4*Ntrack, 3);

  for (int i=0; i<Ntrack; i++){
    int slw_id = track.second[i].first;
    auto slw = Xs_.find(slw_id);
    Pose pose_slw = slw->second;
    M4d Tgb_i;
    Tgb_i << pose_slw.R, pose_slw.p, 0.0, 0.0, 0.0, 1.0;
    M4d Tgcl_i = Tgb_i*Tbc_;
    V4d hpf_cl;
    V4d hpg_c;
    M4d Tgcr_i = Tgcl_i*Tlr_;
    V4d hpf_cr;
    hpg_c << track_info_.f_g, 1.0;
    hpf_cl = Tgcl_i.inverse()*hpg_c;
    hpf_cr = Tgcr_i.inverse()*hpg_c;
    V3d pf_cl = hpf_cl.block<3,1>(0,0);
    V3d pf_cr = hpf_cr.block<3,1>(0,0);
    double Xl = pf_cl(0);
    double Yl = pf_cl(1);
    double Zl = pf_cl(2);
    double Xr = pf_cr(0);
    double Yr = pf_cr(1);
    double Zr = pf_cr(2);

    Eigen::Matrix<double,2,3> Jl_pi, Jr_pi;
    Jl_pi << 1.0/Zl, 0.0, -Xl/pow(Zl,2),
             0.0, 1.0/Zl, -Yl/pow(Zl,2);
    Jr_pi << 1.0/Zr, 0.0, -Xr/pow(Zr,2),
             0.0, 1.0/Zr, -Yr/pow(Zr,2);

    // feature-related Jacobian
    Eigen::Matrix<double,4,3> Hfi_j;
    Hfi_j.block<2,3>(0,0) =
        Jl_pi*(Tgcl_i.block<3,3>(0,0)).transpose();
    Hfi_j.block<2,3>(2,0) =
        Jr_pi*Trl_.block<3,3>(0,0)*(Tgcl_i.block<3,3>(0,0)).transpose();
    H_f.block<4,3>(4*i,0) = Hfi_j;

    // pose-related Jacobian
    // slw index btw 1~MaxSLW
    int slw_idx = Nslw-(frame_id-slw_id);
    V3d temp = track_info_.f_g - Tgb_i.block<3,1>(0,3);
    H_x.block<4,3>(4*i, DXI_ + DXA_ + DXS_*slw_idx) =
        Hfi_j*nesl::Vec2SkewMat(temp);
    H_x.block<4,3>(4*i, DXI_+DXA_+(DXS_*slw_idx)+3) = -Hfi_j;

    // Calibration-related Jacobian
    V4d H_c;
    V3d w_g = Tgb_i.block<3,3>(0,0)*w_hat;
    V3d M = -(nesl::Vec2SkewMat(w_g)*temp - pose_slw.v);
    H_c << Jl_pi*(Tgcl_i.block<3,3>(0,0)).transpose()*M,
           Jr_pi*(Tgcr_i.block<3,3>(0,0)).transpose()*M;
    H_x.block<4,1>(4*i, DXI_) = H_c;
  }

}

V3d MSCKF::triangulateTwoViews
  (const std::pair<int,std::vector<std::pair<int, V2d>>>
    &track){
  // Search for corresponding SLW poses
  // for first and last view
  int first_id = track.second.front().first;
  int last_id = track.second.back().first;
  auto first_SLW = Xs_.find(first_id);
  auto last_SLW = Xs_.find(last_id);
  Pose T1 = first_SLW->second; // first imu
  Pose T2 = last_SLW->second; // last imu

  // Relative pose
  M4d T_gb1;  M4d T_gb2;
  T_gb1 << T1.R, T1.p, 0.0, 0.0, 0.0, 1.0;
  T_gb2 << T2.R, T2.p, 0.0, 0.0, 0.0, 1.0;

  M4d T_gc1;
  M4d T_gc2;
  T_gc1 = T_gb1*Tbc_;
  T_gc2 = T_gb2*Tbc_;

  // T12 = T^{c1}_{c2}
  M4d T12 = T_gc1.inverse()*T_gc2;
  M3d R12 = T12.block<3,3>(0,0);
  V3d t12 = T12.block<3,1>(0,3);

  V3d p;
  p << track.second.front().second.head<2>(),
          1.0;
  V3d q;
  q << track.second.back().second.head<2>(),
          1.0;

  p = p/p.norm();  q = q/q.norm();

  Eigen::Matrix<double, 3,2> A;
  A << p, -R12*q;
  // solve least-squares problem
  V2d lambda = (A.transpose()*A).ldlt().
      solve(A.transpose()*t12);
  V3d f_c0 = lambda(0)*p;

  return f_c0;
}

bool MSCKF::triangulateTrackWithGN
(const std::pair<int,std::vector<std::pair<int, V2d>>>
 &track, V3d f_c0, const double &R_){
  // Multi-view triangulation
  int Ntracks = track.second.size();
  V3d xHat;
  xHat << f_c0(0)/f_c0(2),
          f_c0(1)/f_c0(2),
          1.0/f_c0(2);
  int maxIter = 10;
  double Cprev =
      std::numeric_limits<double>::infinity();

  int first_id = track.second.front().first;
  auto first_SLW = Xs_.find(first_id);
  Pose T1 = first_SLW->second;
  M4d T_gb1;
  T_gb1 << T1.R, T1.p, 0.0, 0.0, 0.0, 1.0;
  M4d T_gc1 = T_gb1*Tbc_;
  double Cnew;
  Eigen::MatrixXd b =
      Eigen::MatrixXd::Zero(2*Ntracks,1);

  // Optimization loop
  for (int i=0; i<maxIter; i++){
    // optimization cost
    double cost = 0.0;

    Eigen::MatrixXd A =
        Eigen::MatrixXd::Zero(2*Ntracks,3);

    // stack Jacobian and residuals
    for (int j=0; j<Ntracks; j++){
      // find j-th view pose
      int j_id = track.second[j].first;
      auto j_SLW = Xs_.find(j_id);
      Pose Tj = j_SLW->second;
      M4d T_gbj;
      T_gbj << Tj.R, Tj.p, 0.0, 0.0, 0.0, 1.0;
      M4d T_gcj = T_gbj*Tbc_;

      // cal. relative pose

      M4d T_j1 = T_gcj.inverse()*T_gc1;

      V4d hl;
      hl << xHat(0), xHat(1), 1.0, xHat(2);
      hl = T_j1*hl;

      // cal. residual
      V2d zHat;
      zHat << hl(0)/hl(2), hl(1)/hl(2);
      V2d z;
      z = track.second[j].second.head<2>();
      b.block<2,1>(2*j,0) = z-zHat;
      cost += (z-zHat).squaredNorm();

      // Jacobian
      // for left features
      Eigen::Matrix<double,2,3> dzl_dg;
      dzl_dg << 1/hl(2), 0.0, -hl(0)/pow(hl(2),2),
                0.0, 1/hl(2), -hl(1)/pow(hl(2),2);
      M3d dgl_dx;
      dgl_dx << T_j1.block<3,2>(0,0), T_j1.block<3,1>(0,3);

      Eigen::Matrix<double, 2,3> Aj;
      Aj << dzl_dg*dgl_dx;
      A.block<2,3>(2*j,0) = Aj;
    }
    Cnew = (0.5/R_)*cost;
    M3d AtA = (1/R_)*A.transpose()*A;
    M3d AtA_diag = M3d::Zero();
    AtA_diag(0,0) = AtA(0,0);
    AtA_diag(1,1) = AtA(1,1);
    AtA_diag(2,2) = AtA(2,2);

    V3d dx = (AtA + (1e-3)*AtA_diag).ldlt().
        solve((1/R_)*A.transpose()*b);
    xHat += dx;
    double Cderiv = std::fabs((Cnew-Cprev)/Cnew);
    Cprev = Cnew;
    if (Cderiv < 1e-6) break;
  }
  V3d f_c1;
  f_c1 << xHat(0)/xHat(2),
          xHat(1)/xHat(2),
          1.0/xHat(2);
  track_info_.f_c1 = f_c1;

  // global featue position
  track_info_.f_g = T_gc1.block<3,3>(0,0)*f_c1
        + T_gc1.block<3,1>(0,3);

  // residuals at each SLW
  track_info_.residual = b;
  // for mono-update
//  Eigen::MatrixXd b2 =
//      Eigen::MatrixXd::Zero(2*Ntracks,1);
//  for (int i=0; i<Ntracks; i++){
//      b2.block<2,1>(2*i,0) = b.block<2,1>(4*i,0);
//  }
//  track_info_.residual = b2;

  if (1/xHat(2) > MAX_DIST_ ||
      1/xHat(2) < MIN_DIST_ ||
      Cnew > 2.0*Ntracks*COST_PROP_){
    // ROS_WARN("Bad Triangulation : throw !");
    // std::cout << 1/xHat(2) << "  " << Cnew << "  " << 2.0*Ntracks*COST_PROP_ << std::endl;
    return false;
  }
  else return true;

}


void MSCKF::calculateMeasJacobian_L
  (const std::pair<int,
    std::vector<std::pair<int, V4d>>> &track,
   Eigen::MatrixXd &H_x,
   Eigen::MatrixXd &H_f,
   int frame_id, V3d w_hat){

  M3d K_;
    K_ << cam_intrinsics_[1],                         0,                                                          0,
                  0,                          cam_intrinsics_[0],                                                 0,
          -cam_intrinsics_[1]*cam_intrinsics_[2],       -cam_intrinsics_[0]*cam_intrinsics_[3],      cam_intrinsics_[0]*cam_intrinsics_[1];
  V3d n_g = track_info_L.L_g.block<3,1>(0,0);
  V3d d_g = track_info_L.L_g.block<3,1>(3,0);

  // Matrix size initialization
  int Ntrack = track.second.size();
  int Nslw = Xs_.size();
  H_x = Eigen::MatrixXd::Zero(2*Ntrack, DXI_ + DXA_ + DXS_*Nslw);
  H_f = Eigen::MatrixXd::Zero(2*Ntrack, 4);

  for (int i=0; i<Ntrack; i++){
    int slw_id = track.second[i].first;
    auto slw = Xs_.find(slw_id);
    Pose pose_slw = slw->second;
    M4d Tgb_i;
    Tgb_i << pose_slw.R, pose_slw.p, 0.0, 0.0, 0.0, 1.0;
    M4d Tgc_i = Tgb_i*Tbc_;
    M4d Tcg_i = Tgc_i.inverse();
    M6d Hcg_i;
    Hcg_i << Tcg_i.block<3, 3>(0, 0), nesl::Vec2SkewMat(Tcg_i.block<3, 1>(0, 3))*Tcg_i.block<3, 3>(0, 0), 
              M3d::Zero(), Tcg_i.block<3, 3>(0, 0);

    V6d L_ci = Hcg_i*track_info_L.L_g;
    V3d l_ci = K_*L_ci.block<3,1>(0,0);
    double ln = l_ci.block<2,1>(0,0).norm();

    V3d ps_, pe_;
    ps_ << track.second[i].second.head<2>(), 1;
    pe_ << track.second[i].second.tail<2>(), 1;
    V2d residuals;
    residuals << ps_.dot(l_ci)/ln, pe_.dot(l_ci)/ln;

    Eigen::Matrix<double,2,3> J_li;
    J_li << ps_[0]-(l_ci[0]*residuals[0])/ln, ps_[1]-(l_ci[1]*residuals[0])/ln, 1,
              pe_[0]-(l_ci[0]*residuals[1])/ln, pe_[1]-(l_ci[1]*residuals[1])/ln, 1;
    J_li *= -(1/ln);

    int slw_idx = Nslw-(frame_id-slw_id);
    H_x.block<2,3>(2*i, DXI_ + DXA_ + DXS_*slw_idx) = J_li*K_*(Tgc_i.block<3,3>(0,0).transpose()*nesl::Vec2SkewMat(n_g)+
                                                              Tgc_i.block<3,3>(0,0).transpose()*nesl::Vec2SkewMat(nesl::Vec2SkewMat(d_g)*Tgb_i.block<3,1>(0,3))+
                                                              nesl::Vec2SkewMat(Tbc_.inverse().block<3,1>(0,3))*Tgc_i.block<3,3>(0,0).transpose()*nesl::Vec2SkewMat(d_g));
    H_x.block<2,3>(2*i, DXI_+DXA_+(DXS_*slw_idx)+3) = J_li*K_*Tgc_i.block<3,3>(0,0).transpose()*nesl::Vec2SkewMat(d_g);


    V3d nd_cross = n_g.cross(d_g);
    M3d U;
    U.block<3, 1>(0, 0) = n_g/n_g.norm();
    U.block<3, 1>(0, 1) = d_g/d_g.norm();
    U.block<3, 1>(0, 2) = nd_cross/nd_cross.norm();

    M2d W;
    W << n_g.norm(), -d_g.norm(),
          d_g.norm(), n_g.norm();
    W *= (1/std::sqrt(n_g.norm()*n_g.norm() + d_g.norm()*d_g.norm()));

    V3d u1 = U.block<3,1>(0,0);
    V3d u2 = U.block<3,1>(0,1);
    V3d u3 = U.block<3,1>(0,2);

    double w1 = W(0, 0);
    double w2 = W(1, 0);

    Eigen::Matrix<double,6,4> dLg_df;
    dLg_df << V3d::Zero(), -w1*u3, w1*u2, -w2*u1,
                w2*u3, V3d::Zero(), -w2*u1, w1*u2;

    Eigen::Matrix<double,3,6> dl_dL;
    dl_dL << K_, M3d::Zero();


    Eigen::Matrix<double,2,4> Hfi_j = J_li*dl_dL*Hcg_i*dLg_df;

    H_f.block<2,4>(2*i,0) = Hfi_j;

  }
}

V6d MSCKF::triangulateTwoViews_L
  (const std::pair<int,std::vector<std::pair<int, V4d>>>
    &track){
  // Search for corresponding SLW poses

  M3d Kmtx = getK();

  V4d first_line = track.second.front().second;
  int first_id = track.second.front().first;
  auto first_SLW = Xs_.find(first_id);
  Pose T1 = first_SLW->second; // first imu
  M4d T_gb1;
  T_gb1 << T1.R, T1.p, 0.0, 0.0, 0.0, 1.0;
  M4d T_gc1 = T_gb1*Tbc_;

  V3d ps_c1 = Kmtx.inverse()*V3d(first_line[0], first_line[1], 1);
  V3d pe_c1 = Kmtx.inverse()*V3d(first_line[2], first_line[3], 1);
  Eigen::MatrixXd L_c1(6, track.second.size() - 1);

  for (size_t i = 1; i < track.second.size(); i++) {
    V4d ith_line = track.second[i].second;
    V3d ps_c2 = Kmtx.inverse()*V3d(ith_line[0], ith_line[1], 1);
    V3d pe_c2 = Kmtx.inverse()*V3d(ith_line[2], ith_line[3], 1);

    V3d n1_c1 = ps_c1.cross(pe_c1);
    V3d n2_c2 = ps_c2.cross(pe_c2);

    int ith_id = track.second[i].first;
    auto ith_SLW = Xs_.find(ith_id);
    Pose Ti = ith_SLW->second;
    M4d T_gbi;
    T_gbi << Ti.R, Ti.p, 0.0, 0.0, 0.0, 1.0;
    M4d T_gci = T_gbi*Tbc_;
    M4d T12 = T_gc1.inverse()*T_gci;

    V4d pi_1, pi_2;
    pi_1 << n1_c1, 0;
    V3d temp = T12.block<3, 3>(0, 0)*n2_c2;
    V3d t_ = T12.block<3, 1>(0, 3);
    pi_2 << temp, -(temp.dot(t_));

    M4d L_star_c1 = pi_1*pi_2.transpose() - pi_2*pi_1.transpose();
    V6d L_c1_temp;
    L_c1_temp << L_star_c1.block<3, 1>(0, 3), nesl::SkewMat2Vec(L_star_c1.block<3, 3>(0, 0));
    L_c1.col(i - 1) = L_c1_temp;
  }

  V6d L_c1_0 = -(L_c1.rowwise().mean());

  return L_c1_0;
}


bool MSCKF::triangulateTrackWithGN_L
(const std::pair<int,std::vector<std::pair<int, V4d>>>
 &track, V6d L_c1_0, const double &R_){
  // Multi-view triangulation
  M3d K_;
    K_ << cam_intrinsics_[1],                         0,                                                          0,
                  0,                          cam_intrinsics_[0],                                                 0,
          -cam_intrinsics_[1]*cam_intrinsics_[2],       -cam_intrinsics_[0]*cam_intrinsics_[3],      cam_intrinsics_[0]*cam_intrinsics_[1];
  V3d n_c1, d_c1, nd_cross, theta_L;
  n_c1 << L_c1_0[0], L_c1_0[1], L_c1_0[2];
  d_c1 << L_c1_0[3], L_c1_0[4], L_c1_0[5];
  nd_cross = n_c1.cross(d_c1);
  M3d temp_mat;
  temp_mat.block<3, 1>(0, 0) = n_c1/n_c1.norm();
  temp_mat.block<3, 1>(0, 1) = d_c1/d_c1.norm();
  temp_mat.block<3, 1>(0, 2) = nd_cross/nd_cross.norm();
  theta_L << nesl::logmap(temp_mat);

  M2d mat_L;
  mat_L << n_c1.norm(), -d_c1.norm(),
            d_c1.norm(), n_c1.norm();
  mat_L *= (1/std::sqrt(n_c1.norm()*n_c1.norm() + d_c1.norm()*d_c1.norm()));
  double phi_L = std::acos(mat_L(0, 0));

  int Ntracks = track.second.size();
  V4d xHat;
  xHat << theta_L, phi_L;
  int maxIter = 10;
  double Cprev =
      std::numeric_limits<double>::infinity();

  int first_id = track.second.front().first;
  auto first_SLW = Xs_.find(first_id);
  Pose T1 = first_SLW->second;
  M4d T_gb1;
  T_gb1 << T1.R, T1.p, 0.0, 0.0, 0.0, 1.0;
  M4d T_gc1 = T_gb1*Tbc_;
  double Cnew;
  Eigen::MatrixXd b =
      Eigen::MatrixXd::Zero(2*Ntracks,1);

  // Optimization loop
  for (int i=0; i<maxIter; i++){
    // optimization cost
    double cost = 0.0;

    Eigen::MatrixXd A =
        Eigen::MatrixXd::Zero(2*Ntracks,4);

    // stack Jacobian and residuals
    for (int j=0; j<Ntracks; j++){
      // find j-th view pose
      int j_id = track.second[j].first;
      auto j_SLW = Xs_.find(j_id);
      Pose Tj = j_SLW->second;
      M4d T_gbj;
      T_gbj << Tj.R, Tj.p, 0.0, 0.0, 0.0, 1.0;
      M4d T_gcj = T_gbj*Tbc_;

      // cal. relative pose
      M4d T_j1 = T_gcj.inverse()*T_gc1;
      M6d H_j1;
      H_j1 << T_j1.block<3, 3>(0, 0), nesl::Vec2SkewMat(T_j1.block<3, 1>(0, 3))*T_j1.block<3, 3>(0, 0), 
                M3d::Zero(), T_j1.block<3, 3>(0, 0);
      M3d U = nesl::expmap(xHat.block<3, 1>(0, 0));

      double w1 = std::cos(xHat[3]);
      double w2 = std::sin(xHat[3]);
      V3d n_hat = w1*U.block<3,1>(0,0);
      V3d d_hat = w2*U.block<3,1>(0,1);

      V6d L_c1;
      L_c1 << n_hat, d_hat;
      V6d L_ci = H_j1*L_c1;
      V3d l_ci = K_*L_ci.block<3,1>(0,0);

      // cal. residual
      V2d residuals;
      V3d ps_, pe_;
      ps_ << track.second[j].second.head<2>(), 1;
      pe_ << track.second[j].second.tail<2>(), 1;
      double ln = l_ci.block<2,1>(0,0).norm();
      residuals << ps_.dot(l_ci)/ln, pe_.dot(l_ci)/ln;
      b.block<2,1>(2*j,0) = residuals;
      cost += (residuals).squaredNorm();

      //ROS_WARN("cost: %f", cost);
      //ROS_WARN("res: %f, %f", residuals[0], residuals[1]);
      

      // Jacobian
      Eigen::Matrix<double,2,3> dr_dl;
      dr_dl << ps_[0]-(l_ci[0]*residuals[0])/ln, ps_[1]-(l_ci[1]*residuals[0])/ln, 1,
                pe_[0]-(l_ci[0]*residuals[1])/ln, pe_[1]-(l_ci[1]*residuals[1])/ln, 1;
      dr_dl *= -(1/ln);

      Eigen::Matrix<double,3,6> dl_dL;
      dl_dL << K_, M3d::Zero();

      Eigen::Matrix<double,6,6> dL_dLc1 = H_j1;

      V3d u1 = U.block<3,1>(0,0);
      V3d u2 = U.block<3,1>(0,1);
      V3d u3 = U.block<3,1>(0,2);

      Eigen::Matrix<double,6,4> dLc1_dL4;
      dLc1_dL4 << V3d::Zero(), -w1*u3, w1*u2, -w2*u1,
                    w2*u3, V3d::Zero(), -w2*u1, w1*u2;

      Eigen::Matrix<double, 2,4> Aj;
      Aj << dr_dl*dl_dL*dL_dLc1*dLc1_dL4;
      A.block<2,4>(2*j,0) = Aj;
    }
    Cnew = (0.5/R_)*cost;
    //ROS_WARN("Cnew, R, cost: %f, %f, %f", Cnew, R_, cost);

    M4d AtA = (1/R_)*A.transpose()*A;
    M4d AtA_diag = M4d::Zero();
    AtA_diag(0,0) = AtA(0,0);
    AtA_diag(1,1) = AtA(1,1);
    AtA_diag(2,2) = AtA(2,2);
    AtA_diag(3,3) = AtA(3,3);

    V4d dx = (AtA + (1e-3)*AtA_diag).ldlt().
        solve((1/R_)*A.transpose()*b);
    
    M3d U_star = nesl::expmap(xHat.block<3,1>(0,0))*nesl::expmap(dx.block<3,1>(0,0));
    M2d W_star, W_xHat, W_dx;
    W_xHat << std::cos(xHat[3]), -std::sin(xHat[3]),
                std::sin(xHat[3]), std::cos(xHat[3]);
    W_dx << std::cos(dx[3]), -std::sin(dx[3]),
                std::sin(dx[3]), std::cos(dx[3]);
    W_star = W_xHat*W_dx;

    V4d xHat_new;
    xHat_new << nesl::logmap(U_star), std::acos(W_star(0,0));

    xHat = xHat_new;

    double Cderiv = std::fabs((Cnew-Cprev)/Cnew);
    Cprev = Cnew;

    
    if (Cderiv < 1e-6) break;
  }


  M3d U = nesl::expmap(xHat.block<3,1>(0,0));
  double w1 = std::cos(xHat[3]);
  double w2 = std::cos(xHat[3]);
  V3d n_hat = w1*U.block<3,1>(0,0);
  V3d d_hat = w2*U.block<3,1>(0,1);

  V6d L_c1;
  L_c1 << n_hat, d_hat;
  track_info_L.L_c1 = L_c1;

  // global featue position
  M6d H_gc1;
      H_gc1 << T_gc1.block<3, 3>(0, 0), nesl::Vec2SkewMat(T_gc1.block<3, 1>(0, 3))*T_gc1.block<3, 3>(0, 0), 
                M3d::Zero(), T_gc1.block<3, 3>(0, 0);

  track_info_L.L_g = H_gc1*L_c1;

  // residuals at each SLW
  track_info_L.residual = b;


  double d = n_hat.norm()/d_hat.norm();


  if (Cnew > 2.0*Ntracks*COST_PROP_ || d > MAX_DIST_){
    // ROS_WARN("Bad Triangulation : throw !");
    // std::cout << 1/xHat(2) << "  " << Cnew << "  " << 2.0*Ntracks*COST_PROP_ << std::endl;
    return false;
  }
  else return true;

}


void MSCKF::zupt(int frame_id){
  // zero velocity update
  int Nstate = cov_X_.cols();
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(3,Nstate);
  H.block<3,3>(0,6) = M3d::Identity();
  M3d R = R_zupt_v*M3d::Identity();
  Eigen::MatrixXd S = H*cov_X_*H.transpose()+R;
  Eigen::MatrixXd K_t = S.ldlt().solve(H*cov_X_);
  Eigen::MatrixXd K = K_t.transpose();
  Eigen::MatrixXd I_KH =
      Eigen::MatrixXd::Identity(K.rows(),H.cols()) -
      K*H;
  Eigen::MatrixXd cov_new = I_KH*cov_X_*I_KH.transpose() +
      K*R*K.transpose();
  V3d r = -Xi_.v;
  Eigen::VectorXd delx_hat = K*r;
  cov_new = validateCovMtx(cov_new);

  // filter update
  cov_X_ = cov_new;

  // IMU/SLW-state update
  UpdateIMUstate(delx_hat.segment<15>(0));
  td_ += delx_hat(15);
  UpdateSLWstate
      (delx_hat.segment(DXI_+DXA_, DXS_*Xs_.size()), frame_id);
  ROS_INFO("Zero velocity updated... !");
}



void MSCKF::zupt_vpq(int frame_id){
  // zero velocity update
  int Nstate = cov_X_.cols();
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(9,Nstate);
  H.block<3,3>(0,6) = M3d::Identity();
  H.block<3,3>(3,3) = M3d::Identity();
  H.block<3,3>(3,Nstate - 3) = -M3d::Identity();
  H.block<3,3>(6,0) = -0.5*M3d::Identity();
  H.block<3,3>(6,Nstate - 6) = 0.5*M3d::Identity();

  //ROS_WARN("H %f, %f, %f, %f, %f", H(6,0), H(3,3), H(0,6), H(6,Nstate-6), H(3, Nstate-3));
  //ROS_WARN("v,p,q: %f, %f, %f", R_zupt_v, R_zupt_p, R_zupt_q);

  Eigen::MatrixXd R = Eigen::MatrixXd::Zero(9, 9); // Initialize a 9x9 matrix filled with zeros
  // Fill the diagonal blocks with the specified values
  R.block<3, 3>(0, 0) = R_zupt_v * Eigen::MatrixXd::Identity(3, 3);
  R.block<3, 3>(3, 3) = R_zupt_p * Eigen::MatrixXd::Identity(3, 3);
  R.block<3, 3>(6, 6) = R_zupt_q * Eigen::MatrixXd::Identity(3, 3);

  Eigen::MatrixXd S = H*cov_X_*H.transpose()+R;
  Eigen::MatrixXd K_t = S.ldlt().solve(H*cov_X_);
  Eigen::MatrixXd K = K_t.transpose();
  Eigen::MatrixXd I_KH =
      Eigen::MatrixXd::Identity(K.rows(),H.cols()) -
      K*H;
  Eigen::MatrixXd cov_new = I_KH*cov_X_*I_KH.transpose() +
      K*R*K.transpose();

  Eigen::VectorXd r(9);
  auto last_element = Xs_.rbegin();
  Pose last_pose = last_element->second;
  //ROS_WARN("last frame: %d", last_element->first);
  V4d dq = nesl::quatMultiply(Xi_.T.q, nesl::quatInverse(last_pose.q));
  r.segment(0, 3) << -Xi_.v;
  r.segment(3, 3) << -(Xi_.T.p - last_pose.p);
  r.segment(6, 3) << dq[1], dq[2], dq[3];

  Eigen::VectorXd delx_hat = K*r;
  cov_new = validateCovMtx(cov_new);

  // filter update
  cov_X_ = cov_new;

  // IMU/SLW-state update
  UpdateIMUstate(delx_hat.segment<15>(0));
  td_ += delx_hat(15);
  UpdateSLWstate
      (delx_hat.segment(DXI_+DXA_, DXS_*Xs_.size()), frame_id);
  ROS_INFO("Zero velocity updated... !");
}

V3d MSCKF::getGravity() {
  return grv_;
}
