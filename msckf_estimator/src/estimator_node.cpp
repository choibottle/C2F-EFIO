/**
 * @file estimator_node.cpp
 * @author Byeongpil Choi (bpc1224@snu.ac.kr)
 * @brief Main node for Event-Frame-Inertial Odometry with point and line features
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

#include <stdio.h>
#include <thread>
#include <mutex>
#include <queue>
#include <condition_variable>
#include <fstream>
#include <tuple>

#include "core/msckf.h"
#include "core/mapmanager.h"

std::condition_variable con;
std::queue<sensor_msgs::ImuConstPtr> imu_buf;
std::queue<sensor_msgs::PointCloudConstPtr> feature_buf;
std::queue<sensor_msgs::PointCloudConstPtr> e_feature_buf;
std::queue<sensor_msgs::NavSatFixConstPtr> gps_buf;

std::vector<sensor_msgs::ImuConstPtr> imu_init_buf;
std::vector<sensor_msgs::NavSatFixConstPtr> gps_init_buf;
std::vector<sensor_msgs::MagneticFieldConstPtr> mag_init_buf;

std::mutex m_buf;  // imu and feature buffer mutex
std::mutex m_estimator;  // msckf mutex
std::mutex m_fast;  // mutex for fast propagation state
std::mutex m_mag;  // mutex for mag buffer

double last_imu_t = 0;
bool init_filter = false;
MSCKF Xmsckf; // set as a global variable
MapManager FeatureMap_f, FeatureMap_e;
int k = 0;
int cnt_img = 0;
bool flag_update = false;
double last_img_t = 0;


// for calcuating time consumption ///////////////////////// 2024. 11. 29
int cnt_it = 0;             // iteration count

double total_time = 0;
double mean_backend = 0;

///////////////////////////////////////////////////////////


bool is_f_dead = false;
bool is_e_dead = false;
bool is_e_dead_L = false;

// Global vars for fast propagator
ros::Publisher pub_fast;
V4d fast_q(1.0, 0.0, 0.0, 0.0);
V3d fast_p(0.0, 0.0, 0.0);
V3d fast_v(0.0, 0.0, 0.0);
V3d fast_ba(0.0, 0.0, 0.0);
V3d fast_bg(0.0, 0.0, 0.0);
V3d fast_grv(0.0, 0.0, 9.81);
double fast_t0 = 0.0;
bool is_fast = false;

// R^{Vehicle}_{IMU}
M3d R_vi;

// Tmp: save execution time
std::ofstream info_fout;

void predict_propagator(const sensor_msgs::ImuConstPtr &imu_msg) {
  V3d f_hat(imu_msg->linear_acceleration.x,
    imu_msg->linear_acceleration.y,
    imu_msg->linear_acceleration.z);
  f_hat -= fast_ba;
  V3d w_hat(imu_msg->angular_velocity.x,
    imu_msg->angular_velocity.y,
    imu_msg->angular_velocity.z);
  w_hat -= fast_bg;
  double t_cur = imu_msg->header.stamp.toSec();
  double dt = t_cur - fast_t0;
  fast_t0 = t_cur;

  M3d Rgb0 = nesl::quat2dcm(fast_q);
  fast_p = fast_p + fast_v*dt + 0.5*dt*dt*(Rgb0*f_hat + fast_grv);
  fast_v = fast_v + dt*(Rgb0*f_hat + fast_grv);
  V3d theta_hat = w_hat * dt;
  V4d del_quat = nesl::rvec2quat(theta_hat);
  fast_q = nesl::quatMultiply(fast_q, del_quat);
}

void update_propagator() {
  NavState Xi = Xmsckf.getNavState();
  fast_q = Xi.T.q;
  fast_p = Xi.T.p;
  fast_v = Xi.v;
  fast_ba = Xi.ba;
  fast_bg = Xi.bg;
  fast_grv = Xmsckf.getGravity();
  // std::cout << "[Update propagator] imu buffer size: " << imu_buf.size() << std::endl;
  std::queue<sensor_msgs::ImuConstPtr> tmp_imu_buf = imu_buf;
  for (sensor_msgs::ImuConstPtr tmp_imu_msg;
    !tmp_imu_buf.empty(); tmp_imu_buf.pop()) {
      predict_propagator(tmp_imu_buf.front());
  }
  is_fast = true;
}

void imu_callback(const sensor_msgs::ImuConstPtr &imu_msg) {
  //ROS_INFO("imu data received");
  if (imu_msg->header.stamp.toSec() <= last_imu_t){
    ROS_WARN("imu messages in disorder!");
    return;
  }

  // Rotate from imu to vehicle frame
  V3d f_raw(imu_msg->linear_acceleration.x,
    imu_msg->linear_acceleration.y,
    imu_msg->linear_acceleration.z);
  
  V3d w_raw(imu_msg->angular_velocity.x,
    imu_msg->angular_velocity.y,
    imu_msg->angular_velocity.z);

  V3d f_rot = R_vi * f_raw;
  V3d w_rot = R_vi * w_raw;

  sensor_msgs::Imu *imu_rot = new sensor_msgs::Imu();
  imu_rot->header = imu_msg->header;
  imu_rot->angular_velocity.x = w_rot(0);
  imu_rot->angular_velocity.y = w_rot(1);
  imu_rot->angular_velocity.z = w_rot(2);

  imu_rot->linear_acceleration.x = f_rot(0);
  imu_rot->linear_acceleration.y = f_rot(1);
  imu_rot->linear_acceleration.z = f_rot(2);

  sensor_msgs::ImuConstPtr imu_ptr(imu_rot);

  // protect imu_buf from other thread
  m_buf.lock();
  imu_buf.push(imu_ptr);
  m_buf.unlock();
  con.notify_one();

  // Fast propagation
  ros::Time start_time = ros::Time::now();
  m_fast.lock();
  if (is_fast) {
    // std::cout << "predict called !!!" << std::endl;
    predict_propagator(imu_msg);
    nav_msgs::Odometry fast_odom;
    fast_odom.header = imu_msg->header;
    fast_odom.header.frame_id = "world";
    fast_odom.child_frame_id = "IMU";
    fast_odom.pose.pose.position.x = fast_p(0);
    fast_odom.pose.pose.position.y = fast_p(1);
    fast_odom.pose.pose.position.z = fast_p(2);
    fast_odom.pose.pose.orientation.w = fast_q(0);
    fast_odom.pose.pose.orientation.x = fast_q(1);
    fast_odom.pose.pose.orientation.y = fast_q(2);
    fast_odom.pose.pose.orientation.z = fast_q(3);
    pub_fast.publish(fast_odom);
  }
  m_fast.unlock();
  double fast_processing_time =
            (ros::Time::now()-start_time).toSec();
  // std::cout << "###########################" << fast_processing_time << std::endl;
}


void feature_callback(const sensor_msgs::PointCloudConstPtr
                      &feature_msg) {
  m_buf.lock();
  feature_buf.push(feature_msg);
  m_buf.unlock();
  con.notify_one();
}

void e_feature_callback(const sensor_msgs::PointCloudConstPtr
                      &e_feature_msg) {
  m_buf.lock();
  e_feature_buf.push(e_feature_msg);
  m_buf.unlock();
  con.notify_one();
}

void gps_callback(const sensor_msgs::NavSatFixConstPtr& gps_msg) {
  // ROS_INFO("GPS called");
  m_buf.lock();
  gps_buf.push(gps_msg);
  m_buf.unlock();
  con.notify_one();
}

void mag_callback(const sensor_msgs::MagneticFieldConstPtr& mag_msg) {
  if (!init_filter) {
    // Use magnetometer just for heading initialization
    m_mag.lock();
    mag_init_buf.push_back(mag_msg);
    m_mag.unlock();
  }
}

std::vector<
  std::tuple<std::vector<sensor_msgs::ImuConstPtr>,
             sensor_msgs::PointCloudConstPtr,
             sensor_msgs::PointCloudConstPtr,
             sensor_msgs::NavSatFixConstPtr>> getMeasurements(){
  std::vector<
    std::tuple<std::vector<sensor_msgs::ImuConstPtr>,
               sensor_msgs::PointCloudConstPtr,
               sensor_msgs::PointCloudConstPtr,
               sensor_msgs::NavSatFixConstPtr>> measurements;

  while(true){

    if (imu_buf.empty() || feature_buf.empty() || e_feature_buf.empty()) break;

    if (!(imu_buf.back()->header.stamp.toSec() >
          feature_buf.front()->header.stamp.toSec() + Xmsckf.getTimeOffset())) break;

    if (!(imu_buf.front()->header.stamp.toSec() <
          feature_buf.front()->header.stamp.toSec() + Xmsckf.getTimeOffset())){
      ROS_WARN("throw img, only happen at the beginning");
      ROS_WARN("22222");
      feature_buf.pop();
      continue;
    }
    sensor_msgs::PointCloudConstPtr img_msg = feature_buf.front();
    double ts_img = img_msg->header.stamp.toSec() + Xmsckf.getTimeOffset();
    feature_buf.pop();

    sensor_msgs::PointCloudConstPtr event_msg = e_feature_buf.front();
    e_feature_buf.pop();

    // Get IMUs
    std::vector<sensor_msgs::ImuConstPtr> IMUs;
    while (imu_buf.front()->header.stamp.toSec() < ts_img){
      IMUs.emplace_back(imu_buf.front());
      imu_buf.pop();
    }
    if (IMUs.empty()) ROS_WARN("no imu btw two images");

    // Get the latest gps message
    sensor_msgs::NavSatFixConstPtr gps_msg;
    if (!gps_buf.empty()) {
      while (gps_buf.front()->header.stamp.toSec() < ts_img) {
        gps_msg = gps_buf.front();
        gps_buf.pop();
        if (gps_buf.empty()) break; 
      }
    }

    std::tuple<std::vector<sensor_msgs::ImuConstPtr>,
      sensor_msgs::PointCloudConstPtr,
      sensor_msgs::PointCloudConstPtr,
      sensor_msgs::NavSatFixConstPtr> temp_tuple
        = std::make_tuple(IMUs, img_msg, event_msg, gps_msg);
    measurements.emplace_back(temp_tuple);

  }
  return measurements;
}

void process(){
  while(true){
    std::vector<
      std::tuple<std::vector<sensor_msgs::ImuConstPtr>,
                 sensor_msgs::PointCloudConstPtr,
                 sensor_msgs::PointCloudConstPtr,
                 sensor_msgs::NavSatFixConstPtr>> measurements;

    // Block process thread until measurements are ready.
    // Escape when the below condition is true and con is notified
    // Lambda: [&]{return}: [&] captures all variables by reference.
    // Implemented in VINS-mono HKUST 2018.
    std::unique_lock<std::mutex> lk(m_buf);
    con.wait(lk, [&]
      { return (measurements = getMeasurements()).size() != 0; });
    lk.unlock();

    // Filter main loop
    m_estimator.lock();
    for (auto& measurement : measurements){
      // Retrieve buffer
      std::vector<sensor_msgs::ImuConstPtr> meas_imu 
        = std::get<0>(measurement);

      sensor_msgs::PointCloudConstPtr meas_img
        = std::get<1>(measurement);

      sensor_msgs::PointCloudConstPtr meas_event
        = std::get<2>(measurement);

      sensor_msgs::NavSatFixConstPtr meas_gps 
            = std::get<3>(measurement);

      // Do filter initialization
      if (!init_filter) {
        imu_init_buf.insert(std::end(imu_init_buf),
          std::begin(meas_imu), std::end(meas_imu));

        if (meas_gps != NULL) {
          gps_init_buf.push_back(meas_gps);
        }

        // Given a sufficient IMU meas
        if (imu_init_buf.size() > Xmsckf.NUM_INIT_SAMPLES_) {
          m_mag.lock();
          std::vector<sensor_msgs::MagneticFieldConstPtr> tmp_mag_buf
            = mag_init_buf;
          m_mag.unlock();
          Xmsckf.InitializeState(imu_init_buf, gps_init_buf, tmp_mag_buf);

          imu_init_buf.clear();
          gps_init_buf.clear();
          mag_init_buf.clear();

          last_imu_t = meas_imu.back()->header.stamp.toSec();
          init_filter = true;
        }
        last_img_t = meas_img->header.stamp.toSec()
          + Xmsckf.getTimeOffset();;
      } else {
        auto img_msg = meas_img;
        auto event_msg = meas_event;

        // INS prediction
        ros::Time start_time = ros::Time::now();
        for (auto &imu_msg : meas_imu){
          double now_imu_t = imu_msg->header.stamp.toSec();
          double imu_dt;
          // imu_dt = now_imu_t - last_imu_t;
          if (flag_update){
              flag_update = false;
              imu_dt = now_imu_t - last_img_t;
          }
          else {
              imu_dt = now_imu_t - last_imu_t;
          }
          last_imu_t = now_imu_t;
          Xmsckf.propagate(imu_msg, imu_dt);
          k++;

          // GPS update if available
          if (meas_gps != NULL && Xmsckf.use_gps_) {
            if (fabs(meas_gps->header.stamp.toSec()-now_imu_t) < 0.5*imu_dt) {
              Xmsckf.correct_position(meas_gps, FeatureMap_f.frame_id_);
            }
          }
        }
        double imu_processing_time =
            (ros::Time::now()-start_time).toSec();

        // Feature Manager
        start_time = ros::Time::now();
        FeatureMap_f.AddFeatures(img_msg);
        FeatureMap_e.AddFeatures_PL(event_msg);
        //ROS_WARN("f mapserver: %zu, e mapserver: %zu", FeatureMap_f.MapServer_.size(), FeatureMap_e.MapServer_.size());
        double feature_processing_time =
            (ros::Time::now()-start_time).toSec();

        // Find deadTracks & update filter
        start_time = ros::Time::now();
        // propagate IMU state to img stamp
        double img_t = img_msg->header.stamp.toSec() +
            Xmsckf.getTimeOffset();
        double imu_t = meas_imu.back()
            ->header.stamp.toSec();
        double gap = img_t - imu_t;
        if (FeatureMap_f.isStatic_){
           Xmsckf.zupt(FeatureMap_f.frame_id_);
           //ROS_WARN("#######  ZUPT  ######");
        }
        is_f_dead = FeatureMap_f.isdeadTracks();
        is_e_dead = FeatureMap_e.isdeadTracks();
        is_e_dead_L = FeatureMap_e.isdeadTracks_L();
        
        if ((is_f_dead || is_e_dead || is_e_dead_L) && Xmsckf.use_vision_){
          flag_update = true;
          Xmsckf.propagate(meas_imu.back(),gap);
          
          Xmsckf.correct(FeatureMap_f.deadTracks_, FeatureMap_e.deadTracks_, FeatureMap_e.deadTracks_L,
                         FeatureMap_f.frame_id_,
                         meas_imu.back());
          
          FeatureMap_f.deadTracks_.clear();
          FeatureMap_e.deadTracks_.clear();
          FeatureMap_e.deadTracks_L.clear();
        }
        else {
          //ROS_WARN("no deadtracks!!!!!!");
        }
        double update_processing_time =
            (ros::Time::now()-start_time).toSec();

        // Augment sliding-window
        start_time = ros::Time::now();
        Xmsckf.AugmentSLW(FeatureMap_f.frame_id_, meas_imu.back());
        ROS_INFO("data seq IMU: %d, IMG: %d", k, FeatureMap_f.frame_id_++);
        FeatureMap_e.frame_id_++;
        double augment_processing_time =
            (ros::Time::now()-start_time).toSec();

        //Publish updated state
        start_time = ros::Time::now();
        std_msgs::Header header = img_msg->header;
        Xmsckf.pubFilterState(header);
        double pub_processing_time =
            (ros::Time::now()-start_time).toSec();

        total_time += imu_processing_time + feature_processing_time
          + update_processing_time + augment_processing_time;
        //std::cout << imu_processing_time << ", " <<
        //            feature_processing_time << ", " <<
        //            update_processing_time << ", " <<
        //            augment_processing_time << ", " <<
        //            pub_processing_time << std::endl;

        // Save execution time
        double time_i = imu_processing_time + feature_processing_time
          + update_processing_time + augment_processing_time;

        last_img_t = img_t;

        cnt_it++;
        mean_backend = (double)total_time/cnt_it;

        ROS_WARN("Total process of backend: %f [ms]", mean_backend*1000);
      }
    }
    // Update fast propagator
    m_buf.lock();
    m_fast.lock();
    if (init_filter) {
      fast_t0 = last_img_t;
      update_propagator();
      // std::cout << "################ Num of IMU buffer AFTER update: " 
      //   << imu_buf.size() << std::endl;
    }
    m_fast.unlock();
    m_buf.unlock();
    m_estimator.unlock();
  }

}

int main(int argc, char **argv)
{
  info_fout.open("/root/catkin_ws/info/estimator_node_time.txt");
  ros::init(argc, argv, "estimator_node");
  ros::NodeHandle n;
  pub_fast = n.advertise<nav_msgs::Odometry>
    ("/smsckf_nesl/fast", 1000);

  Xmsckf.LoadParameters(n);
  R_vi = Xmsckf.Rvi_;
  Xmsckf.registerPub(n);

  FeatureMap_f.LoadParameters(n);
  FeatureMap_e.LoadParameters(n);

  ros::Subscriber sub_imu = n.subscribe("imu", 2000, imu_callback,
                                        ros::TransportHints().tcpNoDelay());

  ros::Subscriber sub_feature = n.subscribe("/mtracker/feature", 100, feature_callback);

  ros::Subscriber sub_feature_e = n.subscribe("/mtracker/feature_e", 100, e_feature_callback);


  std::thread VIO_thread{process};
  ros::spin();

  info_fout.close();
  std::cout << "estimator_node, Total excute time: " << total_time << "[sec]" << std::endl;
  return 0;
}
