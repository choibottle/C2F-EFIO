/**
 * @file mtracker_node.cpp
 * @author Byeongpil Choi (bpc1224@snu.ac.kr)
 * @brief Feature tracking node for event and standrad frame image
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

#include <ros/ros.h>
#include <fstream>
#include <thread>
#include <mutex>
#include <queue>
#include <condition_variable>
#include <algorithm>

#include "mtracker.h"

std::condition_variable con;
std::condition_variable con_pose;
std::mutex m_buf;  // imu and event buffer and img mutex
std::mutex m_frontend;
std::mutex m_pose;

std::queue <dvs_msgs::Event> eventBuffer;
std::queue <sensor_msgs::ImuConstPtr> imuBuffer;
std::queue <sensor_msgs::ImageConstPtr> imgBuffer;
//std::queue <nav_msgs::OdometryConstPtr> poseBuffer;

bool backend_start = false;
bool bias_update_flag = false;

mTracker mtracker;

ros::Publisher pub_match;
ros::Publisher pub_img;

ros::Publisher pub_match_e;
ros::Publisher pub_img_e;


double first_image_time;
bool first_image_flag = true;
double last_image_time = 0;
double last_image_time0 = 0;
bool init_pub = 0;

// for calcuating time consumption ///////////////////////// 2024. 11. 29
int cnt_it = 0;             // iteration count

double total_time = 0;
double mean_frontend = 0;

double framepoint_sum = 0;
double mean_framepointprocess = 0;

//////////////////////////////////////////////////////////


double first_event_time;
double last_event_time = 0;

int next_N = -1;
int real_N;
int win_num;

double t_before_propa;

V3d ba(0.0, 0.0, 0.0);
V3d bg(0.0, 0.0, 0.0);

std::ofstream time_fileOut;


void mono_callback(const sensor_msgs::ImageConstPtr& cam_img) {

    // Throw the first image, only save the timestamp
    if (first_image_flag){
        first_image_flag = false;
        first_image_time = cam_img->header.stamp.toSec();
        last_image_time = cam_img->header.stamp.toSec();
        return;
    }

    // detect unstable camera stream
    if (cam_img->header.stamp.toSec() - last_image_time > 1.0
            || cam_img->header.stamp.toSec() < last_image_time){
        ROS_WARN("Image discontinue! reset the feature tracker !");
        first_image_flag = true;
        last_image_time = 0;

        return;
    }
    last_image_time0 = last_image_time;
    last_image_time = cam_img->header.stamp.toSec();

    m_buf.lock();
    imgBuffer.push(cam_img);
    m_buf.unlock();
    con.notify_one();
    
}

void event_callback(const dvs_msgs::EventArrayConstPtr& cam_events) {
    
    //std::cout << "##### received event ####" << std::endl;

    last_event_time = cam_events->events[0].ts.toSec();

    //ROS_INFO("coming events size: %zu ", cam_events->events.size());

    m_buf.lock();
    for (const auto& event : cam_events->events) {
        eventBuffer.push(event);
    }

    //ROS_INFO("pushed events size: %zu ", eventBuffer.size());
    m_buf.unlock();
    con.notify_one();

}

void imu_callback(const sensor_msgs::ImuConstPtr& imu_msg) {

    //std::cout << "##### received imu ####" << std::endl;
    m_buf.lock();
    imuBuffer.push(imu_msg);
    m_buf.unlock();
    con.notify_one();

}

void pose_callback(const nav_msgs::OdometryConstPtr& pose_msg) {

    //std::cout << "##### received pose ####" << std::endl;

    std::lock_guard<std::mutex> lock(m_pose);
    t_before_propa = pose_msg->header.stamp.toSec();

    ba(0, 0) = pose_msg->pose.covariance[19];
    ba(1, 0) = pose_msg->pose.covariance[20];
    ba(2, 0) = pose_msg->pose.covariance[21];

    bg(0, 0) = pose_msg->pose.covariance[22];
    bg(1, 0) = pose_msg->pose.covariance[23];
    bg(2, 0) = pose_msg->pose.covariance[24];

    if (!backend_start) {
        backend_start = true;
        ROS_WARN("###########  First pose arrives ########### ");
    } 
    bias_update_flag = true;
    //ROS_INFO("arrived pose timestamp: %f", pose_msg->header.stamp.toSec());
    con_pose.notify_all();
}


std::vector<std::tuple<
    std::vector <dvs_msgs::Event>,
    std::vector <sensor_msgs::ImuConstPtr>,
    sensor_msgs::ImageConstPtr>>
getMeasurements() {

    std::vector<std::tuple<
        std::vector <dvs_msgs::Event>,
        std::vector <sensor_msgs::ImuConstPtr>,
        sensor_msgs::ImageConstPtr>>
        measurements;

    while (true) {
        if (imuBuffer.empty() || eventBuffer.empty() || imgBuffer.empty()) break;

        if (!(eventBuffer.back().ts.toSec() > imgBuffer.front()->header.stamp.toSec()) || 
            !(imuBuffer.back()->header.stamp.toSec() > imgBuffer.front()->header.stamp.toSec())) break;
        
        if (!(eventBuffer.front().ts.toSec() < imgBuffer.front()->header.stamp.toSec()) ||
            !(imuBuffer.front()->header.stamp.toSec() < imgBuffer.front()->header.stamp.toSec())) {
            ROS_WARN("throw img, only happen at the beginning");
            ROS_WARN("111111");
            imgBuffer.pop();
            continue;
        }


        // Get images
        sensor_msgs::ImageConstPtr img_msg = imgBuffer.front();
        double ts_img = img_msg->header.stamp.toSec();
        imgBuffer.pop();

        // Get IMUs
        std::vector<sensor_msgs::ImuConstPtr> IMUs;
        while (imuBuffer.front()->header.stamp.toSec() < ts_img) {
            IMUs.emplace_back(imuBuffer.front());
            imuBuffer.pop();
        }

        // Get events
        std::vector<dvs_msgs::Event> events;
        while (eventBuffer.front().ts.toSec() < ts_img) {
            events.emplace_back(eventBuffer.front());
            eventBuffer.pop();
        }
        
        std::tuple<std::vector<dvs_msgs::Event>,
                    std::vector<sensor_msgs::ImuConstPtr>,
                    sensor_msgs::ImageConstPtr> temp_tuple
            = std::make_tuple(events, IMUs, img_msg);
        measurements.emplace_back(temp_tuple);
    }
    return measurements;
}


void process() {
   
    while (true) {

        std::vector<std::tuple<
            std::vector <dvs_msgs::Event>,
            std::vector <sensor_msgs::ImuConstPtr>,
            sensor_msgs::ImageConstPtr>>
            measurements;

        std::unique_lock<std::mutex> lk(m_buf);
        con.wait(lk, [&] 
            { return (measurements = getMeasurements()).size() !=0; });       
        lk.unlock();
        m_frontend.lock();

        for (auto &measurement : measurements) {
            std::vector<dvs_msgs::Event> meas_ev = 
                std::get<0>(measurement);

            std::vector<sensor_msgs::ImuConstPtr> meas_imu =
                std::get<1>(measurement);

            sensor_msgs::ImageConstPtr meas_img =
                std::get<2>(measurement);

            ros::Time start_time = ros::Time::now();


            // Adaptive event accumulation
            // Determine the number of events to accumulate
            if (mtracker.next_t_ev == -1) {
                next_N = mtracker.getMaxEvents();
            }
            else {
                if (mtracker.next_t_ev >= meas_img->header.stamp.toSec()) {
                    next_N = meas_ev.size();
                    //ROS_INFO("next_N: %d ", next_N);
                }
                else {
                    std::vector<dvs_msgs::Event> temp_ev;
                    for (const auto& ev : meas_ev) {
                        double ev_time = ev.ts.toSec();                   
                        if (ev_time <= mtracker.next_t_ev) {
                            temp_ev.push_back(ev);
                        }
                    }
                    //ROS_INFO("temp_ev size: %zu ", temp_ev.size());
                    if (temp_ev.size() < 0.4*mtracker.getMaxEvents()) {
                        next_N = 0.4*mtracker.getMaxEvents();
                    }
                    else if (temp_ev.size() > mtracker.getMaxEvents()) {
                        next_N = mtracker.getMaxEvents();
                    }
                    else {
                        next_N = temp_ev.size();
                    }
                }
            }

            int M_events = meas_ev.size();
            if (next_N > M_events) {
                real_N = M_events;
                win_num = 1;
            }
            else {
                win_num = M_events/next_N;
                real_N = M_events/win_num;
            }

            //ROS_INFO("M, next_N, win_num, real_N: %d, %d, %d, %d ", M_events, next_N, win_num, real_N);
            //ROS_WARN("win_num: %d", win_num);


            // Synthesize eventframes and event feature tracking
            cv::Mat eventFrame;

            std::unique_lock<std::mutex> lock(m_pose);

            if (!backend_start) {
                if (win_num == 1) {
                    mtracker.several_windows = false;
                    eventFrame = mtracker.makeEventFrame(meas_ev, meas_imu, ba, bg);
                    V4d q_12(1, 0, 0, 0);
                    mtracker.readEventImage(eventFrame, meas_img->header.stamp.toSec(), meas_ev, q_12);
                    for (unsigned int i = 0;; i++){
                        bool completed = false;
                        completed |= mtracker.updateID_e(i);
                        if (!completed)
                            break;
                    }
                }
                else {
                    mtracker.several_windows = true;
                    for (int i = 0; i < win_num; i++) {
                        auto start_iter = meas_ev.begin() + i * real_N;
                        auto end_iter = (i == win_num - 1) ? meas_ev.end() : (meas_ev.begin() + (i + 1) * real_N);
                        std::vector<dvs_msgs::Event> meas_ev_i(start_iter, end_iter);

                        double start_time = meas_ev_i.front().ts.toSec();
                        double end_time = meas_ev_i.back().ts.toSec();

                        // Create a vector for the current slice of meas_imu
                        std::vector<sensor_msgs::ImuConstPtr> meas_imu_i;
                        for (const auto& imu : meas_imu) {
                            double imu_time = imu->header.stamp.toSec();
                            if (imu_time >= start_time && imu_time <= end_time) {
                                meas_imu_i.push_back(imu);
                            }
                        }
                        eventFrame = mtracker.makeEventFrame(meas_ev_i, meas_imu_i, ba, bg);
                        if (i == win_num - 1) {
                            V4d q_12(1, 0, 0, 0);
                            mtracker.readEventImage(eventFrame, meas_img->header.stamp.toSec(), meas_ev_i, q_12);
                        }
                        else {
                            mtracker.readEventImage_no_cnt(eventFrame, end_time);
                        }
                        for (unsigned int j = 0;; j++){
                            bool completed = false;
                            completed |= mtracker.updateID_e(j);
                            if (!completed)
                                break;
                        }
                    }
                }
            }
            else {
                con_pose.wait(lock, [] {return bias_update_flag;});

                if (win_num == 1) {
                    mtracker.several_windows = false;
                    eventFrame = mtracker.makeEventFrame(meas_ev, meas_imu, ba, bg);
                    V4d q_12 = mtracker.IMU_propa_q(meas_imu, t_before_propa, ba, bg);      // calculate R_12
                    mtracker.readEventImage(eventFrame, meas_img->header.stamp.toSec(), meas_ev, q_12);
                    for (unsigned int i = 0;; i++){
                        bool completed = false;
                        completed |= mtracker.updateID_e(i);
                        if (!completed)
                            break;
                    }
                }
                else {
                    mtracker.several_windows = true;
                    for (int i = 0; i < win_num; i++) {
                        auto start_iter = meas_ev.begin() + i * real_N;
                        auto end_iter = (i == win_num - 1) ? meas_ev.end() : (meas_ev.begin() + (i + 1) * real_N);
                        std::vector<dvs_msgs::Event> meas_ev_i(start_iter, end_iter);

                        double start_time = meas_ev_i.front().ts.toSec();
                        double end_time = meas_ev_i.back().ts.toSec();

                        // Create a vector for the current slice of meas_imu
                        std::vector<sensor_msgs::ImuConstPtr> meas_imu_i;
                        for (const auto& imu : meas_imu) {
                            double imu_time = imu->header.stamp.toSec();
                            if (imu_time >= start_time && imu_time <= end_time) {
                                meas_imu_i.push_back(imu);
                            }
                        }
                        eventFrame = mtracker.makeEventFrame(meas_ev_i, meas_imu_i, ba, bg);


                        if (i == win_num - 1) {
                            V4d q_12 = mtracker.IMU_propa_q(meas_imu, t_before_propa, ba, bg);          // calculate R_12
                            mtracker.readEventImage(eventFrame, meas_img->header.stamp.toSec(), meas_ev_i, q_12);
                        }
                        else {
                            mtracker.readEventImage_no_cnt(eventFrame, end_time);
                        }
                        for (unsigned int j = 0;; j++){
                            bool completed = false;
                            completed |= mtracker.updateID_e(j);
                            if (!completed)
                                break;
                        }
                    }
                }
                bias_update_flag = false;
            }


            // Compute next temporal window
            if (mtracker.velocity_use) {
                mtracker.calTempWindow();
                mtracker.next_t_ev = meas_img->header.stamp.toSec() + mtracker.t_window;
                //ROS_WARN("t_window, t_ev: %f, %f", mtracker.t_window, mtracker.next_t_ev);
            }

            

            // to opencv data structure            
            cv_bridge::CvImageConstPtr ptr;
            ptr = cv_bridge::toCvCopy(meas_img, sensor_msgs::image_encodings::MONO8);
            cv::Mat show_img = ptr->image;


            ros::Time start_point = ros::Time::now();
            

            // Standard frame image point feature tracking
            mtracker.readMonoImage(ptr->image.rowRange(0, mtracker.h_img_), meas_img->header.stamp.toSec());
            // update feature id
            for (unsigned int i = 0;; i++){
                bool completed = false;
                completed |= mtracker.updateID(i);
                if (!completed)
                    break;
            }

            double framepoint_processing_time =
                (ros::Time::now()-start_point).toSec();

            double img_processing_time =
                (ros::Time::now()-start_time).toSec();

            // Input to the filter
            start_time = ros::Time::now();

            sensor_msgs::PointCloudPtr feature_points(new sensor_msgs::PointCloud);
            sensor_msgs::ChannelFloat32 id_of_point;
            sensor_msgs::ChannelFloat32 u_of_point; // left undistorted u
            sensor_msgs::ChannelFloat32 v_of_point; // left undistorted v
            sensor_msgs::ChannelFloat32 velocity_x_of_point; // left of
            sensor_msgs::ChannelFloat32 velocity_y_of_point; // left
            feature_points->header = meas_img->header;
            feature_points->header.frame_id = "world";

            auto &un_pts = mtracker.cur_un_pts_;
            auto &cur_pts = mtracker.cur_pts_;
            auto &ids = mtracker.ids_;
            auto &pts_velocity = mtracker.pts_velocity_;
            
            for (unsigned int j = 0; j < ids.size(); j++)
            {
                if (mtracker.track_cnt_[j] > 1)
                {
                    int p_id = ids[j];
                    geometry_msgs::Point32 p;
                    p.x = un_pts[j].x;
                    p.y = un_pts[j].y;
                    p.z = 1;

                    feature_points->points.push_back(p);
                    id_of_point.values.push_back(p_id);
                    u_of_point.values.push_back(cur_pts[j].x);
                    v_of_point.values.push_back(cur_pts[j].y);
                    velocity_x_of_point.values.push_back(pts_velocity[j].x);
                    velocity_y_of_point.values.push_back(pts_velocity[j].y);
                }
            }
            feature_points->channels.push_back(id_of_point);
            feature_points->channels.push_back(u_of_point);
            feature_points->channels.push_back(v_of_point);
            feature_points->channels.push_back(velocity_x_of_point);
            feature_points->channels.push_back(velocity_y_of_point);


            // for eventframe ////////////////////////////////////////////////////////////////////////////////////////////////
            sensor_msgs::PointCloudPtr feature_points_e(new sensor_msgs::PointCloud);
            sensor_msgs::ChannelFloat32 id_of_point_e;
            sensor_msgs::ChannelFloat32 u_of_point_e; // left undistorted u
            sensor_msgs::ChannelFloat32 v_of_point_e; // left undistorted v
            sensor_msgs::ChannelFloat32 velocity_x_of_point_e; // left of
            sensor_msgs::ChannelFloat32 velocity_y_of_point_e; // left
            sensor_msgs::ChannelFloat32 un_ps_x;
            sensor_msgs::ChannelFloat32 un_ps_y;
            sensor_msgs::ChannelFloat32 un_pe_x;
            sensor_msgs::ChannelFloat32 un_pe_y;
            feature_points_e->header = meas_img->header;
            feature_points_e->header.frame_id = "world";

            auto &un_pts_e = mtracker.cur_un_pts_e;
            auto &cur_pts_e = mtracker.cur_pts_e;
            auto &ids_e = mtracker.ids_e;
            auto &pts_velocity_e = mtracker.pts_velocity_e; 
            mtracker.un_lines_map = convertToMap(mtracker.un_lines);
            auto &un_lines = mtracker.un_lines_map;

            for (unsigned int j = 0; j < ids_e.size(); j++)
            {
                if (mtracker.track_cnt_e[j] > 1)
                {
                    int p_id = ids_e[j];
                    geometry_msgs::Point32 p;
                    p.x = un_pts_e[j].x;
                    p.y = un_pts_e[j].y;
                    p.z = 1;

                    feature_points_e->points.push_back(p);
                    id_of_point_e.values.push_back(p_id);
                    u_of_point_e.values.push_back(cur_pts_e[j].x);
                    v_of_point_e.values.push_back(cur_pts_e[j].y);
                    velocity_x_of_point_e.values.push_back(pts_velocity_e[j].x);
                    velocity_y_of_point_e.values.push_back(pts_velocity_e[j].y);

                    V4d un_line;
                    if (findMapValue(un_lines, p_id, un_line)) {
                        un_ps_x.values.push_back(un_line[0]);
                        un_ps_y.values.push_back(un_line[1]);
                        un_pe_x.values.push_back(un_line[2]);
                        un_pe_y.values.push_back(un_line[3]);
                    }
                    else {
                        un_ps_x.values.push_back(-100);
                        un_ps_y.values.push_back(-100);
                        un_pe_x.values.push_back(-100);
                        un_pe_y.values.push_back(-100);
                    }
                }
            }
            feature_points_e->channels.push_back(id_of_point_e);
            feature_points_e->channels.push_back(u_of_point_e);
            feature_points_e->channels.push_back(v_of_point_e);
            feature_points_e->channels.push_back(velocity_x_of_point_e);
            feature_points_e->channels.push_back(velocity_y_of_point_e);
            feature_points_e->channels.push_back(un_ps_x);
            feature_points_e->channels.push_back(un_ps_y);
            feature_points_e->channels.push_back(un_pe_x);
            feature_points_e->channels.push_back(un_pe_y);


            // skip the first image; since no optical speed on frist image
            if (!init_pub)
            {
                init_pub = 1;
            }
            else
            {
                pub_img.publish(feature_points);
                pub_img_e.publish(feature_points_e);
            }
            double feature_processing_time =
                (ros::Time::now()-start_time).toSec();
            
            // draw left and right image with feature
            start_time = ros::Time::now();

            ptr = cv_bridge::cvtColor(ptr, sensor_msgs::image_encodings::BGR8);
            cv::Mat mono_img = ptr->image;
            cv::Mat tmp_img = mono_img.rowRange(0, mtracker.h_img_);
            cv::cvtColor(show_img, tmp_img, CV_GRAY2RGB);

            cv::Mat mono_img_e;
            cv::cvtColor(eventFrame, mono_img_e, cv::COLOR_GRAY2BGR);
            cv::Mat tmp_img_e = mono_img_e.rowRange(0, mtracker.h_img_);

            // draw points
            for (unsigned int j = 0; j < mtracker.cur_pts_.size(); j ++){
                double len = std::min(1.0, 1.0 * mtracker.track_cnt_[j] / mtracker.wsize_);
                cv::circle(tmp_img, mtracker.cur_pts_[j], 1, cv::Scalar(255*(1 - len), 255 * len, 0), -1);
            }
            for (unsigned int j = 0; j < mtracker.cur_pts_e.size(); j ++){
                double len = std::min(1.0, 1.0 * mtracker.track_cnt_e[j] / mtracker.wsize_);
                cv::circle(tmp_img_e, mtracker.cur_pts_e[j], 1, cv::Scalar(255*(1 - len), 255 * len, 0), -1);
            }

            // draw lines
            for (const auto& line : mtracker.cur_id_point_line) {
                V4d line_ = std::get<2>(line);
                if (isNonZero(line_)) {
                    cv::line(tmp_img_e, cv::Point2f(line_[0], line_[1]), cv::Point2f(line_[2], line_[3]), cv::Scalar(0, 0, 255), 1);
                }
            }
            // for (const auto& line : mtracker.final_lines) {
            //     V4d line_ = line.second;
            //     cv::line(tmp_img_e, cv::Point2f(line_[0], line_[1]), cv::Point2f(line_[2], line_[3]), cv::Scalar(255, 255, 0), 1);

            // }
            double draw_processing_time =
                (ros::Time::now()-start_time).toSec();
            
            // only for visualization
            pub_match.publish(ptr->toImageMsg());
            cv_bridge::CvImage cv_image(std_msgs::Header(), sensor_msgs::image_encodings::BGR8, mono_img_e);
            sensor_msgs::ImagePtr ros_image = cv_image.toImageMsg();
            pub_match_e.publish(ros_image);

            total_time += img_processing_time + feature_processing_time;
            //std::cout << img_processing_time << ", " <<
            //              feature_processing_time << ", " <<
            //              draw_processing_time << std::endl;
            framepoint_sum += framepoint_processing_time;
                        
            // Save execution time
            double time_i = img_processing_time + feature_processing_time;
            //time_fileOut.precision(6);
            //time_fileOut << time_i << std::endl;
            //time_fileOut << img_processing_time << " " << feature_processing_time << " " << win_num <<std::endl;

            cnt_it++;
            mean_frontend = (double)total_time/cnt_it;
            mean_framepointprocess = (double)framepoint_sum/cnt_it;

            ROS_WARN("Frame-point detection & tracking: %f [ms]", mean_framepointprocess*1000);
            ROS_WARN("Total process of frontend: %f [ms]", mean_frontend*1000);

        }

        m_frontend.unlock();

    }

}

int main(int argc, char **argv)
{
  time_fileOut.open("/root/catkin_ws/info/mtracker_node_time.txt");
  ros::init(argc, argv, "mtracker_node");
  ros::NodeHandle n;

  mtracker.LoadParameters(n);

  
  // Subscribe to pose from backend
  ros::Subscriber pose_sub = n.subscribe("/smsckf_nesl/odom", 1000, pose_callback);

  // Subscribe to events
  ros::Subscriber event_sub = n.subscribe("events", 1000, event_callback);

  // Subscribe to imu
  ros::Subscriber imu_sub = n.subscribe("imu", 2000, imu_callback,
                                        ros::TransportHints().tcpNoDelay());

  // Subscribe to mono image
  ros::Subscriber cam_sub = n.subscribe("left_img", 10, mono_callback);

  // Publish feature points to backend
  pub_img = n.advertise<sensor_msgs::PointCloud>("mtracker/feature", 1000);
  pub_img_e = n.advertise<sensor_msgs::PointCloud>("mtracker/feature_e", 1000);

  // Publish topics for debugging
  pub_match = n.advertise<sensor_msgs::Image>("mtracker/feature_img", 1000);
  pub_match_e = n.advertise<sensor_msgs::Image>("mtracker/feature_img_e", 1000);


  std::thread frontend_thread{process};
  ros::spin();

  time_fileOut.close();
  std::cout << "mtracker_node, Execution time: " << total_time << "[sec]" << std::endl;
  return 0;
}
