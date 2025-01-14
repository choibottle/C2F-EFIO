/**
 * @file msckf.h
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


#ifndef MSCKF_H
#define MSCKF_H

#include <ros/ros.h>
#include <visualization_msgs/Marker.h>
#include <sensor_msgs/Imu.h>
#include <sensor_msgs/PointCloud.h>
#include <sensor_msgs/NavSatFix.h>
#include <sensor_msgs/MagneticField.h>
#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include <tf/transform_broadcaster.h>
#include <sensor_msgs/PointCloud.h>
#include <sensor_msgs/PointCloud2.h>
#include <sensor_msgs/point_cloud2_iterator.h>
#include <boost/math/distributions/chi_squared.hpp>

#include "../utils/define.h"
#include "../utils/Attitude.h"
#include "../utils/utils.h"
#include "../utils/NavEarth.h"

extern ros::Publisher pub_state;
extern ros::Publisher pub_path;

extern nav_msgs::Path path;

class MSCKF
{
private:
  NavState Xi_; // 15-dim IMU-state

  // 6-dim SLW-state<frame_id, Pose>
  std::map<int, Pose> Xs_;

  int DXI_; // dim. of IMU-state
  int DXS_; // dim. of SLW-state
  int DXA_; // dim. of Additional state
  V3d grv_; // gravity
  Eigen::MatrixXd cov_X_; // system's cov mtx
  M12d Q_; // process noise matrix
  double R_pix_f; // meas. noise [pixel^2]
  double R_pix_e; // meas. noise [pixel^2]
  double R_pix_e_L; // meas. noise [pixel^2]

  double R_zupt_v; // meas. noise [(m/s)^2]
  double R_zupt_p; // meas. noise [(m/s)^2]
  double R_zupt_q; // meas. noise [(m/s)^2]
  int MAX_SLW_; // max sliding-window
  int MIN_SLW_; // min sliding-window
  M4d Tbc_, Trl_, Tlr_;
  double MAX_DIST_, MIN_DIST_, COST_PROP_;
  double td_;  // cam-imu timeoffset: t_i = t_c + td_
  V3d ecef0_;  // Starting position in ECEF [m]

  struct trackINFO{
    V3d f_c1;
    V3d f_g; // global feature position
    Eigen::MatrixXd residual;
  };

  struct trackINFO_L{
    V6d L_c1;
    V6d L_g; // global feature position
    Eigen::MatrixXd residual;
  };

  // tracking information for a single feature
  trackINFO track_info_;

  trackINFO_L track_info_L;

  std::vector<double> cam_intrinsics_;

  // MSCKF points for drawing
  std::vector<V3d> msckf_points_;

  // Chi squared test table
  static std::map<int, double> chi_squared_test_table_;


  // IMU Euler integration
  NavState propagateIMU_Euler
    (V3d fib_b, V3d wib_b, double dt);

  // cov_X propagate
  Eigen::MatrixXd propagateCovX(V3d fib_b, double dt);

  // Validate cov matrix
  // if cov mtx is not positive definite
  // make it as positive definite
  // by modifying its eigenvalues
  Eigen::MatrixXd validateCovMtx(Eigen::MatrixXd P_mtx);

  V3d triangulateTwoViews
    (const std::pair<int,
      std::vector<std::pair<int, V2d>>> &track);

  // Feature triangulation within all SLWs
  // input : single feature track,
  //         initial guess
  // output : true/false meaning
  //          triangulated feature is usable or not
  bool triangulateTrackWithGN
    (const std::pair<int,
    std::vector<std::pair<int, V2d>>> &track,
   V3d f_c0, const double &R_);

  // calculate Measurement Jacobian
  void calculatesMeasJacobian
    (const std::pair<int,
      std::vector<std::pair<int, V4d>>> &track,
     Eigen::MatrixXd &H_x,
     Eigen::MatrixXd &H_f,
     int frame_id_, V3d w_hat);

  void calculateMeasJacobian
    (const std::pair<int,
      std::vector<std::pair<int, V2d>>> &track,
     Eigen::MatrixXd &H_x,
     Eigen::MatrixXd &H_f,
     int frame_id_, V3d w_hat);


  void calculateMeasJacobian_L
  (const std::pair<int,
    std::vector<std::pair<int, V4d>>> &track,
   Eigen::MatrixXd &H_x,
   Eigen::MatrixXd &H_f,
   int frame_id, V3d w_hat);


  V6d triangulateTwoViews_L(const std::pair<int,std::vector<std::pair<int, V4d>>> &track);

  bool triangulateTrackWithGN_L (const std::pair<int,std::vector<std::pair<int, V4d>>> &track, V6d L_c1_0, const double &R_);

  // Chi-squared test before correcting
  bool ChiSquaredTest(double gamma, int dof);

  // Filter update computation
  void MeasurementUpdate(const Eigen::MatrixXd &Ho,
                        const Eigen::VectorXd &ro,
                        int frame_id, const int &f_size, const int &e_size, const int &e_L_size);

  // Update IMU related state (15-dim)
  void UpdateIMUstate(V15d delx);

  // Update Sliding Window related state (6N-dim)
  void UpdateSLWstate(Eigen::VectorXd delx,
                      int frame_id);

public:
  // Constructor
  MSCKF();

  M3d getK();

  // Load params - init P0 and Q
  void LoadParameters(ros::NodeHandle &n);

  // Initialize filter states
  void InitializeState(std::vector<sensor_msgs::ImuConstPtr>& imu_buf,
    std::vector<sensor_msgs::NavSatFixConstPtr>& gps_buf,
    std::vector<sensor_msgs::MagneticFieldConstPtr>& mag_buf);

  // INS prediction
  void propagate(const sensor_msgs::ImuConstPtr &imu_msg,
                 double dt);

  void correct_position(const sensor_msgs::NavSatFixConstPtr& gps_msg,
    const int& frame_id);

  // Get filter states
  NavState getNavState();

  // Register ros publishers
  void registerPub(ros::NodeHandle &n);

  // Publish filter state
  void pubFilterState(const std_msgs::Header &header);

  // Augment SLW state
  void AugmentSLW(int frame_id, const sensor_msgs::ImuConstPtr &imu_msg);

  // MSCKF update
  void correct
    (const std::map<int, std::vector<std::pair<int, V2d>>> &deadTracks_f, 
      const std::map<int, std::vector<std::pair<int, V2d>>> &deadTracks_e, 
        const std::map<int, std::vector<std::pair<int, V4d>>> &deadTracks_e_L, int frame_id, const sensor_msgs::ImuConstPtr &imu_msg);

  // Zero velocity update
  void zupt(int frame_id);
  void zupt_vpq(int frame_id);

  V3d getGravity();

  inline double getTimeOffset() { return td_; }

  int NUM_INIT_SAMPLES_;

  M3d Rvi_;

  bool use_vision_;
  bool use_gps_;

  inline int getNumUpdatePoints() { return msckf_points_.size(); }

};

#endif // MSCKF_H
