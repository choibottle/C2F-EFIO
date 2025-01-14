/**
 * @file mapmanager.h
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

#ifndef MAPMANAGER_H
#define MAPMANAGER_H

#include "msckf.h"

class MapManager
{
private:
  // Minimum sliding-window
  int MIN_SLW_;

  // Maximum sliding-window
  int MAX_SLW_;

  // zupt threshold [pixel]
  double thr_zupt_;

  // camera fps
  double fps_;

  // focal lengh
  double fu_, fv_;

public:
  MapManager();

  // is the scene is static ?
  bool isStatic_;

  // frame sequence id from 1
  int frame_id_;

  // Set parameters
  void LoadParameters(ros::NodeHandle &n);

  // map(feature_id, vector(frame_id, normalized_pts))
  std::map<int, std::vector<std::pair<int, V2d>>> MapServer_;
  std::map<int, std::vector<std::pair<int, V2d>>> deadTracks_;

  std::map<int, std::vector<std::pair<int, V4d>>> MapServer_L;
  std::map<int, std::vector<std::pair<int, V4d>>> deadTracks_L;


  // Add features to the MapServer
  void AddFeatures
    (sensor_msgs::PointCloudConstPtr img_msg);

  void AddFeatures_PL
    (sensor_msgs::PointCloudConstPtr img_msg);


  // detect deadTracks
  // Output : true if available
  //          false if not
  bool isdeadTracks();

  bool isdeadTracks_L();
};

#endif // MAPMANAGER_H
