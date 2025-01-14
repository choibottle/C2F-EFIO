/**
 * @file mtracker.h
 * @author Byeongpil Choi (bpc1224@snu.ac.kr)
 * @brief Feature tracking class for event and standrad frame image
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


#ifndef MTRACKER_H
#define MTRACKER_H

#include <ros/ros.h>
#include <dvs_msgs/Event.h>
#include <dvs_msgs/EventArray.h>
#include <sensor_msgs/Imu.h>
#include <sensor_msgs/Image.h>
#include <sensor_msgs/PointCloud.h>
#include <nav_msgs/Odometry.h>
#include <message_filters/subscriber.h>
#include <message_filters/time_synchronizer.h>
#include <cv_bridge/cv_bridge.h>

#include <opencv2/opencv.hpp>
#include <opencv2/features2d.hpp>
#include <opencv2/xfeatures2d.hpp>
#include <opencv2/imgproc.hpp>
#include <eigen3/Eigen/Dense>
#include <random>
#include <ctime>
#include <cmath>
#include <omp.h>
#include <unordered_map>
#include <tuple>

#include "../../msckf_estimator/src/utils/define.h"
#include "../../msckf_estimator/src/utils/Attitude.h"
#include "../../msckf_estimator/src/utils/utils.h"
#include "../../msckf_estimator/src/utils/NavEarth.h"

typedef struct {
    double first;
    std::vector<std::vector<double>> second;
} CM_pair;

typedef struct {
    V3d first;
    std::vector<cv::Point2f> second;
    int third;
} ransac_pair;

bool inBorder(const cv::Point2f &pt, int COL, int ROW);
void reduceVector(std::vector<cv::Point2f> &v, std::vector<uchar> status);
void reduceVector(std::vector<int> &v, std::vector<uchar> status);

bool isNonZero(const V4d& v);
std::unordered_map<int, V4d> convertToMap(const std::vector<std::pair<int, V4d>>& vec);     // convert vec to map
bool findMapValue(const std::unordered_map<int, V4d>& map, int j, V4d& value);        // find a value in a map

class mTracker
{
private:
    // left camera parameters: fu, fv, cu, cv
    std::vector<double> cam_intrinsics_;

    // left camera distortion coefficients
    // k1, k2, p1, p2
    std::vector<double> cam_dcf_;

    // Maximum number of features
    // set in .launch file
    int MAX_CNT_;

    // Minimum distant btw features
    // set in .launch file
    int MIN_DIST_;

    // Maximum number of events
    int MAX_EVENTS_;

    double LIFE_PIXEL_;


    // line detection parameters
    double ETA_;
    int PATCH_SIZE_;
    double PIX_THR_;
    int RANSAC_ITER_;

    // line matching parameters
    double DIST_THR_;
    double THETA_THR_;

    // radtan or equidistant: left and right should be the same
    std::string DISTORTION_MODEL_;

public:
  // constructor
  mTracker();
  int wsize_; // for drawing
  int w_img_; // width of image
  int h_img_; // height of image
  double cur_time_; // current ts of tracker
  double prev_time_; // previous ts of tracker

  double cur_time_e; // current ts of tracker
  double prev_time_fix; // for velocity calculation
  double prev_time_e; // previous ts of tracker

  bool several_windows;

  M4d Tbc_;

  // img0 = left, img1 = right
  cv::Mat cur_img_, prev_img_;
  std::vector<cv::Point2f> cur_pts_, prev_pts_;

  // tracked points at left
  std::vector<cv::Point2f> n_pts_;

  // mask for corner extraction
  cv::Mat mask_;

  // feature id
  std::vector<int> ids_;
  std::vector<int> track_cnt_;

  // undistorted feature points for publish
  std::vector<cv::Point2f> cur_un_pts_;

  // for calculating optical flow based on left
  std::map<int, cv::Point2f> cur_un_pts_map_, prev_un_pts_map_;
  std::vector<cv::Point2f> pts_velocity_;


  // for eventFrame ////////////////////////////////////////////////////////////////////////////

  cv::Mat makeEventFrame(const std::vector<dvs_msgs::Event>& evs, const std::vector<sensor_msgs::ImuConstPtr>& imus, const V3d& ba, const V3d& bg);

  // synchronize current previous lines according to ID
  std::vector<std::tuple<int, cv::Point2f, cv::Point2f, V4d, V4d>> synchronizeLines(const std::vector<std::tuple<int, cv::Point2f, V4d>> &currentLines,
                                                                                      const std::vector<std::tuple<int, cv::Point2f, V4d>> &previousLines);

  ransac_pair ransac_line(const std::vector<cv::Point2f>& data_points, const int& max_iter, const double& max_dist);

  CM_pair get_contrast(const cv::Point2f& x, const std::vector<dvs_msgs::Event>& ce, const std::vector<int>& x_length, const std::vector<int>& y_length);

  double line_search(const cv::Point2f& x, const cv::Point2f& dx, const cv::Point2f& df, const double& a, const double& b, 
                      const std::vector<dvs_msgs::Event>& ce, const std::vector<int>& x_length, const std::vector<int>& y_length);

  std::vector<std::vector<double>> CM(const std::vector<dvs_msgs::Event>& events, const cv::Point2f& feature_pos, const cv::Point2f& flow_init, const int& wsize);

  std::vector<cv::Point2f> sort_points(const std::vector<std::vector<double>>& sharp_frame, const int& wsize);

  std::vector<cv::Point2f> proj_endpoints(const std::vector<cv::Point2f>& inliers, const V3d& normal);

  std::tuple<cv::Point2f, cv::Point2f, int> extend_endpoints(const cv::Point2f& p_s, const cv::Point2f& p_e, const cv::Mat &img_raw, 
                                                            const std::vector<dvs_msgs::Event>& evs, const cv::Point2f& flow_init, const int& wsize, 
                                                            const int& max_iter, const V3d& normal_);

  std::vector<std::vector<V3d>> extend_and_intensity_single(const cv::Mat &img_raw, const cv::Point2f& p_s, const cv::Point2f& p_e, const double& intensity_threshold);

  std::vector<V3d> extend_line_single(const cv::Mat &img_raw, const cv::Point2f& point, const cv::Point2f& direction, const double& intensity_threshold);

  double get_intensity_mean(const cv::Mat &img_raw, const cv::Point2f& pointRounded, const cv::Point2f& direction);

  double get_intensity(const cv::Mat &img_raw, const cv::Point2f& pixel);

  bool is_outside_image(const cv::Point2f& point, const cv::Mat &img_raw);

  // calculate R_12
  V4d IMU_propa_q(const std::vector<sensor_msgs::ImuConstPtr>& IMUs, const double& before_t, const V3d& ba, const V3d& bg);


  cv::Mat cur_img_e, prev_img_e;
  std::vector<cv::Point2f> cur_pts_e, prev_pts_e;
  std::vector<cv::Point2f> before_pts_e, after_pts_e;   // for no count tracking
  std::vector<V4d> cur_lines_e;

  std::vector<std::tuple<int, cv::Point2f, V4d>> prev_id_point_line, cur_id_point_line;
  std::vector<std::pair<int, V4d>> final_lines;     // final matched lines

  // tracked points at left
  std::vector<cv::Point2f> n_pts_e;

  // mask for corner extraction
  cv::Mat mask_e;

  // feature id
  std::vector<int> ids_e;
  std::vector<int> track_cnt_e;

  // undistorted feature points for publish
  std::vector<cv::Point2f> cur_un_pts_e;
  std::vector<std::pair<int, V4d>> un_lines;      // undistored final_lines
  std::unordered_map<int, V4d> un_lines_map;      // undistored final_lines map

  // for calculating optical flow based on left
  std::map<int, cv::Point2f> cur_un_pts_map_e, prev_un_pts_map_e;
  std::vector<cv::Point2f> pts_velocity_e;

  double next_t_ev = -1;
  double t_window = -1;

  bool velocity_use = false;

  // for calculating time consumption in 'mtracker.cpp'  //////////// 2024. 11. 29
  int cnt_eventpoint = 0;
  double eventpoint_sum = 0;
  double mean_eventpointprocess = 0; 

  int cnt_line_ex = 0;
  double line_extract_sum = 0;
  double mean_line_extact = 0; 

  int cnt_line_mat = 0;
  double line_matching_sum = 0;
  double mean_line_matching = 0; 

  int cnt_eventframe = 0;
  double eventframe_sum = 0;
  double mean_eventframe = 0; 
  ////////////////////////////////////////////////////////////////////////////////

  M3d getK();

  int getMaxCnt();

  int getMaxEvents();

  // IMU-CAM ex/intrinsic params
  void LoadParameters(ros::NodeHandle &n);

  void readMonoImage(const cv::Mat &img_raw, double cur_time);

  void readEventImage(const cv::Mat &img_raw, double cur_time, const std::vector<dvs_msgs::Event>& evs, const V4d& q);

  void readEventImage_no_cnt(const cv::Mat &img_raw, double cur_time);

  void calTempWindow();

  // set mask for uniform distribution based on left image
  void setMask();

  void setMask_e();

  void setMask_e_point();

  void setMask_e_no_cnt();

  // stack tracked features
  void addPoints();

  void addPoints_e();

   void addPoints_e_no_cnt();

  // update feature ID
  bool updateID(unsigned int i);

  bool updateID_e(unsigned int i);

  // RANSAC: WHICH = 0:Left, 1:Right
  void rejectWithF();

  void rejectWithF_e();

  void rejectWithF_e_no_cnt();

  // undistortion wrapper
  void liftProjective(const Eigen::Vector2d& p,
                      Eigen::Vector3d& P) const;

  // distortion
  void distortion(const Eigen::Vector2d& p_u,
                  Eigen::Vector2d& d_u,
                  const Eigen::Vector4d& dist_coeff) const;
                  
  // feature undistortion
  cv::Point2f undistort_a_point(const cv::Point2f& p);

  void undistortedPoints();

  void undistortedPoints_e();

  void undistortedLines_e();

  static int n_id;

  static int n_id_e;

};

#endif // MTRACKER_H
