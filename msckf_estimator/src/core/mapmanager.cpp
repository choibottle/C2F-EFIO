/**
 * @file mapmanager.cpp
 * @author Byeongpil Choi (bpc1224@snu.ac.kr)
 * @brief Point and line feature track manager
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

#include "mapmanager.h"

MapManager::MapManager()
{
  frame_id_ = 1;
  isStatic_ = false;
}

void MapManager::LoadParameters(ros::NodeHandle &n){
  // MSCKF params
  n.param<int>("/estimator_node/MinSLW", MIN_SLW_, 3);
  n.param<int>("/estimator_node/MaxSLW", MAX_SLW_, 5);

  n.param<double>("estimator_node/thr_zupt",
                   thr_zupt_, 0.7);

  n.param<double>("estimator_node/fps", fps_, 20.0);

  std::vector<double> cam0_intrinsics(4);
  n.getParam("/estimator_node/cam0/intrinsics",
              cam0_intrinsics);

  fu_ = cam0_intrinsics[0];
  fv_ = cam0_intrinsics[1];
}

void MapManager::AddFeatures
(sensor_msgs::PointCloudConstPtr img_msg){
  double mean_of = 0.0;
  double size_valid = 0;
  unsigned int temp_left_num = img_msg->points.size();
  for (unsigned int i = 0; i < temp_left_num; i++){
    int feature_id = img_msg->channels[0].values[i];
    double xl = img_msg->points[i].x;
    double yl = img_msg->points[i].y;
//    double zr = img_msg_->points[i].z;
    double v_x = img_msg->channels[3].values[i];
    double v_y = img_msg->channels[4].values[i];
//    ROS_ASSERT(zr==1);
    V2d temp_pts;
    temp_pts << xl, yl;
    std::pair<int, V2d> temp_pair;
    temp_pair = std::make_pair(frame_id_, temp_pts);
//    std::cout << temp_pair.first << ", " << temp_pair.second << std::endl;
    MapServer_[feature_id].emplace_back(temp_pair);

    if (!(v_x == 0 && v_y == 0)) {
      V2d of;
      of << v_x, v_y;
      mean_of += of.norm();
      size_valid++;
    }
  }
  mean_of /= size_valid;
  // #TODO read from yaml file
  // if mean of Optical Flow is less than threshold
  // thr_zupt * (fps/focal length)
  if (mean_of < thr_zupt_*(fps_/(fu_+fv_)) && frame_id_ > 1 && size_valid > 10){
    isStatic_ = true;
  }
  else{
    isStatic_ = false;
  }
  //ROS_WARN("mean_of: %f", mean_of);
  //ROS_WARN("num of features: %f", size_valid);

//  auto sample = MapServer_.find(4); // find feature_id : 4
//  auto track = sample->second;
//  std::cout << "feature id : " << sample->first << ", "
//            << "frame id : " << track[0].first <<  ", "
//            << "track sample : " << track[0].second << std::endl;
//  std::cout << "feature id : " << sample->first << ", "
//            << "track length : " << track.size() << std::endl;

//  ROS_INFO("%d-frame stacked at MapServer_",frame_id);
}


void MapManager::AddFeatures_PL
(sensor_msgs::PointCloudConstPtr img_msg){

  unsigned int temp_left_num = img_msg->points.size();
  for (unsigned int i = 0; i < temp_left_num; i++){
    int feature_id = img_msg->channels[0].values[i];
    double xl = img_msg->points[i].x;
    double yl = img_msg->points[i].y;

    V2d temp_pts;
    temp_pts << xl, yl;
    std::pair<int, V2d> temp_pair;
    temp_pair = std::make_pair(frame_id_, temp_pts);
//    std::cout << temp_pair.first << ", " << temp_pair.second << std::endl;
    MapServer_[feature_id].emplace_back(temp_pair);

    // check line is valid
    if (img_msg->channels[5].values[i] != -100) {
      double ps_x = img_msg->channels[5].values[i];
      double ps_y = img_msg->channels[6].values[i];
      double pe_x = img_msg->channels[7].values[i];
      double pe_y = img_msg->channels[8].values[i];

      V4d temp_line;
      temp_line << ps_x, ps_y, pe_x, pe_y;
      std::pair<int, V4d> temp_pair_line;
      temp_pair_line = std::make_pair(frame_id_, temp_line);

      MapServer_L[feature_id].emplace_back(temp_pair_line);

    }
  }

  //ROS_WARN("mean_of: %f", mean_of);
  //ROS_WARN("num of features: %f", size_valid);

//  auto sample = MapServer_.find(4); // find feature_id : 4
//  auto track = sample->second;
//  std::cout << "feature id : " << sample->first << ", "
//            << "frame id : " << track[0].first <<  ", "
//            << "track sample : " << track[0].second << std::endl;
//  std::cout << "feature id : " << sample->first << ", "
//            << "track length : " << track.size() << std::endl;

//  ROS_INFO("%d-frame stacked at MapServer_",frame_id);
}

bool MapManager::isdeadTracks(){
  bool update_flag = false;
  std::vector<int> removed_ids;
  for (auto &track : MapServer_){
    int feature_id = track.first;
    // Erase garbage tracks
    // tracking failed && less than MinSLW
    if (track.second.size() < MIN_SLW_ &&
        track.second.back().first == frame_id_-1){
      removed_ids.push_back(feature_id);
      //ROS_WARN("track early fail");
      continue;
    }

    // 1st criteria for deadTracks
    // tracking failed after MinSLW
    if (track.second.back().first == frame_id_-1 &&
        track.second.size() >= MIN_SLW_ &&
        track.second.size() <= MAX_SLW_){
      // add to deadTracks
      auto temp_track = MapServer_.find(feature_id);
      deadTracks_.insert
          (std::pair<int, std::vector<std::pair<int, V2d>>>
           (temp_track->first, temp_track->second));

      // remove from MapServer_
      removed_ids.push_back(feature_id);

      update_flag = true;
      continue;
    }

    // 2nd critera for deadTracks
    // keep tracking but
    // track length exceeds MaxSLW
    if (track.second.back().first == frame_id_ &&
        track.second.size() == MAX_SLW_+1){
//      std::cout << track.first  << ", " << track.second.size() << std::endl;
      // add to deadTracks
      auto temp_track_ = MapServer_.find(feature_id);

      // #TODO : handling continued track
      // but now just remove the latest one
      temp_track_->second.pop_back();
      deadTracks_.insert
          (std::pair<int, std::vector<std::pair<int, V2d>>>
           (temp_track_->first, temp_track_->second));

      // remove from MapServer_
      removed_ids.push_back(feature_id);

      update_flag = true;
    }
  }

  // Removing deadTracks from MapServer_
  for (int i = 0; i < removed_ids.size(); i++){
    MapServer_.erase(removed_ids[i]);
  }

  return update_flag;
}


bool MapManager::isdeadTracks_L(){
  bool update_flag = false;
  std::vector<int> removed_ids;
  for (auto &track : MapServer_L){
    int feature_id = track.first;
    // Erase garbage tracks
    // tracking failed && less than MinSLW
    if (track.second.size() < MIN_SLW_ &&
        track.second.back().first == frame_id_-1){
      removed_ids.push_back(feature_id);
      //ROS_WARN("track early fail");
      continue;
    }

    // 1st criteria for deadTracks
    // tracking failed after MinSLW
    if (track.second.back().first == frame_id_-1 &&
        track.second.size() >= MIN_SLW_ &&
        track.second.size() <= MAX_SLW_){
      // add to deadTracks
      auto temp_track = MapServer_L.find(feature_id);
      deadTracks_L.insert
          (std::pair<int, std::vector<std::pair<int, V4d>>>
           (temp_track->first, temp_track->second));

      // remove from MapServer_
      removed_ids.push_back(feature_id);

      update_flag = true;
      continue;
    }

    // 2nd critera for deadTracks
    // keep tracking but
    // track length exceeds MaxSLW
    if (track.second.back().first == frame_id_ &&
        track.second.size() == MAX_SLW_+1){
//      std::cout << track.first  << ", " << track.second.size() << std::endl;
      // add to deadTracks
      auto temp_track_ = MapServer_L.find(feature_id);

      // #TODO : handling continued track
      // but now just remove the latest one
      temp_track_->second.pop_back();
      deadTracks_L.insert
          (std::pair<int, std::vector<std::pair<int, V4d>>>
           (temp_track_->first, temp_track_->second));

      // remove from MapServer_
      removed_ids.push_back(feature_id);

      update_flag = true;
    }
  }

  // Removing deadTracks from MapServer_
  for (int i = 0; i < removed_ids.size(); i++){
    MapServer_L.erase(removed_ids[i]);
  }

  return update_flag;
}









