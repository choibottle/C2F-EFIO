/**
 * @file mtracker.cpp
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


#include "mtracker.h"

int mTracker::n_id = 0;
int mTracker::n_id_e = 0;

bool inBorder(const cv::Point2f &pt, int COL, int ROW)
{
    const int BORDER_SIZE = 1;
    int img_x = cvRound(pt.x);
    int img_y = cvRound(pt.y);
    return BORDER_SIZE <= img_x && img_x < COL - BORDER_SIZE && BORDER_SIZE <= img_y && img_y < ROW - BORDER_SIZE;
}

void reduceVector(std::vector<cv::Point2f> &v, std::vector<uchar> status)
{
    int j = 0;
    for (int i = 0; i < int(v.size()); i++)
        if (status[i])
            v[j++] = v[i];
    v.resize(j);
}

void reduceVector(std::vector<int> &v, std::vector<uchar> status)
{
    int j = 0;
    for (int i = 0; i < int(v.size()); i++)
        if (status[i])
            v[j++] = v[i];
    v.resize(j);
}

bool isNonZero(const V4d& v) {
    return v[0] != 0 || v[1] != 0 || v[2] != 0 || v[3] != 0;
}

std::unordered_map<int, V4d> convertToMap(const std::vector<std::pair<int, V4d>>& vec) {
    std::unordered_map<int, V4d> map;
    for (const auto& pair : vec) {
    map[pair.first] = pair.second;
    }
    return map;
}

bool findMapValue(const std::unordered_map<int, V4d>& map, int j, V4d& value) {
    auto it = map.find(j);
    if (it != map.end()) {
        value = it->second;
        return true;
    }
    return false;
}

mTracker::mTracker()
{
    ;
}

M3d mTracker::getK() {
    M3d K;
    K << cam_intrinsics_[0],      0,        cam_intrinsics_[2],
            0,         cam_intrinsics_[1],     cam_intrinsics_[3],
            0,                0,               1;
    
    return K;
}

int mTracker::getMaxCnt() {
    return MAX_CNT_;
}

int mTracker::getMaxEvents() {
    return MAX_EVENTS_;
}

std::vector<std::tuple<int, cv::Point2f, cv::Point2f, V4d, V4d>> mTracker::synchronizeLines(const std::vector<std::tuple<int, cv::Point2f, V4d>>& currentLines,
                                                                                        const std::vector<std::tuple<int, cv::Point2f, V4d>>& previousLines) {
    // Create a map for the previous lines using IDs as keys
    std::map<int, std::tuple<cv::Point2f, V4d>> prevLinesMap;
    for (const auto& prevLine : previousLines) {
        int id = std::get<0>(prevLine);
        prevLinesMap[id] = std::make_tuple(std::get<1>(prevLine), std::get<2>(prevLine));
    }

    // Vector to store matched and synchronized elements
    std::vector<std::tuple<int, cv::Point2f, cv::Point2f, V4d, V4d>> synchronizedLines;

    // Iterate through the current lines and synchronize with previous lines by ID
    for (const auto& curLine : currentLines) {
        int id = std::get<0>(curLine);
        if (prevLinesMap.find(id) != prevLinesMap.end()) {
            // Retrieve previous line data
            cv::Point2f prevPt = std::get<0>(prevLinesMap[id]);
            V4d prevLine = std::get<1>(prevLinesMap[id]);

            // Retrieve current line data
            cv::Point2f curPt = std::get<1>(curLine);
            V4d curLineData = std::get<2>(curLine);

            // Store synchronized lines (or adjust as needed)
            synchronizedLines.push_back(std::make_tuple(id, prevPt, curPt, prevLine, curLineData));
        }
    }

    return synchronizedLines;
}


ransac_pair mTracker::ransac_line(const std::vector<cv::Point2f>& data_points, const int& max_iter, const double& max_dist) {
    int N = data_points.size();
    int n1, n2;

    const unsigned int fixed_seed = 100;
    std::mt19937 gen(fixed_seed);
    std::uniform_int_distribution<int> dist(0, N-1);
    
    ransac_pair result;
    V3d best_normal(0.0, 0.0, 0.0);
    std::vector<cv::Point2f> best_inliers;
    int most_inliers = 0;

    cv::Point2f A1, A2, normal;
    double c, norm;

    std::vector<cv::Point2f> inLier;
    inLier.reserve(N);

    for (int i = 0; i < max_iter; i++) {
        n1 = dist(gen);
        n2 = dist(gen);
        while (n2 == n1) {
            n2 = dist(gen);
        }
        A1 = data_points[n1];
        A2 = data_points[n2];
        normal.x = -(A2.y - A1.y);
        normal.y = A2.x - A1.x;
        c = A1.x * A2.y - A2.x * A1.y;

        norm = cv::norm(normal);

        inLier.clear();

        for (const auto& pt : data_points) {
            double distance = std::abs((normal.x * pt.x + normal.y * pt.y + c) / norm);
            
            if (distance < max_dist) {
                inLier.push_back(pt);
            }
        }

        if (inLier.size() > most_inliers) {
            best_normal[0] = normal.x / norm;
            best_normal[1] = normal.y / norm;
            best_normal[2] = c / norm;
            best_inliers = inLier;
            most_inliers = inLier.size();
        }

        if (static_cast<double>(most_inliers) / N > 0.7) {
            //ROS_WARN("early stop!!!! with %d iter", i);
            break;
        }
        
    }

    result.first = best_normal;
    result.second = best_inliers;
    result.third = most_inliers;

    return result;
}

CM_pair mTracker::get_contrast(const cv::Point2f& x, const std::vector<dvs_msgs::Event>& ce, const std::vector<int>& x_length, const std::vector<int>& y_length) {
    int wsize = x_length.size();
    int N = ce.size();

    CM_pair result;
    std::vector<std::vector<double>> H(wsize, std::vector<double>(wsize, 0));

    double t_ref = ce.back().ts.toSec();
    double avg_H = 0;

    int x_min = x_length[0];
    int x_max = x_length[wsize - 1];
    int y_min = y_length[0];
    int y_max = y_length[wsize - 1];

    for (int i = 0; i < N; i++) {
        double wdt = ce[i].ts.toSec() - t_ref;
        double wx = ce[i].x - x.x*wdt;
        double wy = ce[i].y - x.y*wdt;

        if (wx <= x_max && wx >= x_min && wy <= y_max && wy >= y_min) {
            int idx_x = static_cast<int>(cvRound(wx - x_length[0]));
            int idx_y = static_cast<int>(cvRound(wy - y_length[0]));

            H[idx_y][idx_x] += 1;
            avg_H += 1;

        }
    }

    double f = 0;
    avg_H /= (wsize*wsize);
    
    for (int i = 0; i < wsize; i++) {
        for (int j = 0; j < wsize; j++) {
            double diff = H[i][j] - avg_H;
            f += diff * diff;
        }
    }
    f /= (wsize*wsize);

    result.first = f;
    result.second = H;

    return result;
}

double mTracker::line_search(const cv::Point2f& x, const cv::Point2f& dx, const cv::Point2f& df, const double& a, const double& b, 
                      const std::vector<dvs_msgs::Event>& ce, const std::vector<int>& x_length, const std::vector<int>& y_length) {
    double t = 1;
    cv::Point2f tmp;

    double initial_contrast = -get_contrast(x, ce, x_length, y_length).first;

    for (int i = 0; i < 50; i++) {
        tmp.x = x.x + t * dx.x;
        tmp.y = x.y + t * dx.y;

        double new_contrast = -get_contrast(tmp, ce, x_length, y_length).first;

        if (new_contrast < initial_contrast + a * t * (df.x * dx.x + df.y * dx.y)) {
            break;
        }
        t = b * t;
    }
    return t;
}

std::vector<std::vector<double>> mTracker::CM(const std::vector<dvs_msgs::Event>& events, const cv::Point2f& feature_pos, const cv::Point2f& flow_init, const int& wsize) {
    std::vector<std::vector<double>> H;
    std::vector<int> x_range, y_range;
    cv::Point2f x = flow_init;

    x_range.reserve(wsize);
    y_range.reserve(wsize);
    for (int i = -(wsize - 1)/2; i <= (wsize - 1)/2; i++) {
        x_range.push_back(i);
        y_range.push_back(i);
    }

    std::vector<dvs_msgs::Event> final_events;
    for (const auto& e : events) {
        dvs_msgs::Event centered_event = e;
        centered_event.x -= feature_pos.x;
        centered_event.y -= feature_pos.y;
        if (centered_event.x >= x_range[0] && centered_event.y >= y_range[0] && centered_event.x <= x_range[wsize-1] && centered_event.y <= y_range[wsize-1]) {
            final_events.push_back(centered_event);
        }
    }

    if (final_events.size() < 50) {
        return H;
    }

    int max_iter = 10;
    double jde = 30;
    int fade_rate = 1;
    std::vector<cv::Point2f> x_history;
    x_history.reserve(max_iter);

    cv::Point2f G, G_;
    double aa = 0.003;
    double bb = 0.77;

    for (int i = 0; i < max_iter; i++) {
        cv::Point2f x1, x0;
        double f1, f0, df_dx, df_dy;

        // calculate gradient in x direction
        x1 = x + cv::Point2f(jde, 0);
        x0 = x - cv::Point2f(jde, 0);

        f1 = -get_contrast(x1, final_events, x_range, y_range).first;
        f0 = -get_contrast(x0, final_events, x_range, y_range).first;
        df_dx = (f1 - f0)/(2 * jde);

        // calculate gradient in y direction
        x1 = x + cv::Point2f(0, jde);
        x0 = x - cv::Point2f(0, jde);

        f1 = -get_contrast(x1, final_events, x_range, y_range).first;
        f0 = -get_contrast(x0, final_events, x_range, y_range).first;
        df_dy = (f1 - f0)/(2 * jde);

        jde = jde * fade_rate;
        G = cv::Point2f(df_dx, df_dy);
        G_ = -G;

        double line_et = line_search(x, G_, G, aa, bb, final_events, x_range, y_range);

        x -= line_et * G;
        x_history.push_back(x);

        if ((i>1) && (cv::norm(G)<0.003 || cv::norm(x_history[i]-x_history[i-1])/cv::norm(x_history[i-1]) < 0.0005) && (cv::norm(x) > 4) && (cv::norm(x) < 1000)) {
            break;
        }
    }

    std::vector<std::vector<double>> H_ = get_contrast(x, final_events, x_range, y_range).second;

    // find max and normalize H_
    H.reserve(H_.size());

    double max_ = 0;
    for (const auto& row : H_) {
        double row_max = *std::max_element(row.begin(), row.end());
        if (row_max > max_) {
            max_ = row_max;
        }
    }

    for (const auto& row : H_) {
        std::vector<double> tmp;
        tmp.reserve(row.size());
        for (double value : row) {
            tmp.push_back(value / max_);
        }
        H.push_back(std::move(tmp));
    }

    return H;
}

std::vector<cv::Point2f> mTracker::sort_points(const std::vector<std::vector<double>>& sharp_frame, const int& wsize) {
    std::vector<std::pair<cv::Point2f, double>> H_pts_with_intensity;
    if (sharp_frame.empty()) {
        ROS_WARN("empty!!!!!!!!!!!!!!!");
    }
    for (int j = 0; j < wsize; j++) {
        for (int k = 0; k < wsize; k++) {
            if (sharp_frame[j][k] > PIX_THR_) {
                double col = static_cast<double>(k - (wsize - 1) / 2);
                double row = static_cast<double>(j - (wsize - 1) / 2);
                H_pts_with_intensity.push_back(std::make_pair(cv::Point2f(col, row), sharp_frame[j][k]));
            }
        }
    }

    // Sort H_pts based on intensity values in descending order
    std::sort(H_pts_with_intensity.begin(), H_pts_with_intensity.end(),
        [](const std::pair<cv::Point2f, double>& a, const std::pair<cv::Point2f, double>& b) {
            return a.second > b.second;
        });

    // Extract sorted points
    std::vector<cv::Point2f> H_pts;
    for (const auto& p : H_pts_with_intensity) {
        H_pts.push_back(p.first);
    }

    return H_pts;
}

std::vector<cv::Point2f> mTracker::proj_endpoints(const std::vector<cv::Point2f>& inliers, const V3d& normal) {
    cv::Point2f p_s0, p_e0;
    std::vector<double> line_pts_x, line_pts_y;
    
    // extract x and y coordinates of inliers
    
    for (int j = 0; j < inliers.size(); j++) {
        line_pts_x.push_back(inliers[j].x);
        line_pts_y.push_back(inliers[j].y);
    }

    // determine start and end point of line segment
    if (std::abs(normal[1]/normal[0]) >= 1) {
        // vertical line case
        auto max_it_x = std::max_element(line_pts_x.begin(), line_pts_x.end());
        auto min_it_x = std::min_element(line_pts_x.begin(), line_pts_x.end());
        p_s0 = inliers[std::distance(line_pts_x.begin(), max_it_x)];
        p_e0 = inliers[std::distance(line_pts_x.begin(), min_it_x)];
    }
    else {
        // horizontal line case
        auto max_it_y = std::max_element(line_pts_y.begin(), line_pts_y.end());
        auto min_it_y = std::min_element(line_pts_y.begin(), line_pts_y.end());
        p_s0 = inliers[std::distance(line_pts_y.begin(), max_it_y)];
        p_e0 = inliers[std::distance(line_pts_y.begin(), min_it_y)];
    }

    // compute the direction and normalize it
    double norm = cv::norm(cv::Point2f(normal[0], normal[1]));
    double dist_s_signed = (normal[0] * p_s0.x + normal[1] * p_s0.y + normal[2]) / norm;
    cv::Point2f p_s = p_s0 - (dist_s_signed/norm) * cv::Point2f(normal[0], normal[1]);

    double dist_e_signed = (normal[0] * p_e0.x + normal[1] * p_e0.y + normal[2]) / norm;
    cv::Point2f p_e = p_e0 - (dist_e_signed/norm) * cv::Point2f(normal[0], normal[1]);

    std::vector<cv::Point2f> res;
    res.push_back(p_s);
    res.push_back(p_e);

    return res;

}

std::tuple<cv::Point2f, cv::Point2f, int> mTracker::extend_endpoints(const cv::Point2f& p_s, const cv::Point2f& p_e, const cv::Mat &img_raw, 
                                                                    const std::vector<dvs_msgs::Event>& evs, const cv::Point2f& flow_init, const int& wsize, 
                                                                    const int& max_iter, const V3d& normal_) {
    
    int extend_cnt = 0;
    cv::Point2f p_s_temp = p_s;
    cv::Point2f p_e_temp = p_e;

    //std::cout << "p_s: " << p_s_temp << ", p_e: " << p_e_temp << std::endl;

    double intensity_threshold = 0.15;

    std::vector<std::vector<V3d>> extended_lines = extend_and_intensity_single(img_raw, p_s, p_e, intensity_threshold);

    //ROS_WARN("ex_line size: %zu, %zu", extended_lines.front().size(), extended_lines.back().size());

    if (extended_lines.front().size() >= 4) {
        cv::Point2f ps_new(extended_lines.front().back()[0], extended_lines.front().back()[1]);
        std::vector<std::vector<double>> extended_frame_s = CM(evs, ps_new, flow_init, wsize);
        if (!extended_frame_s.empty()) {
            std::vector<cv::Point2f> H_pts = sort_points(extended_frame_s, wsize);  

            if (H_pts.size() >= 10) {
                ransac_pair ransac_result = ransac_line(H_pts, max_iter, 1);
                V3d normal = ransac_result.first;
                std::vector<cv::Point2f> inliers = ransac_result.second;

                double inlier_ratio = static_cast<double>(ransac_result.third) / H_pts.size();
                //ROS_WARN("inlier_ratio: %f ", inlier_ratio);

                cv::Point2f v1(normal_[0], normal_[1]);
                cv::Point2f v2(normal[0], normal[1]);

                double temp = (v1.x * v2.x + v1.y * v2.y) / (cv::norm(v1) * cv::norm(v2));
                if (temp > 1.0) {
                                temp = 1.0;
                }
                if (temp < -1.0) {
                                temp = -1.0;
                }
                double theta = R2D * acos(temp);
                double min_theta = std::min(theta, 180 - theta);

                cv::Point2f p_s_new;
                if (inlier_ratio > 0.7 && min_theta < 10) {

                    std::vector<cv::Point2f> endpoints = proj_endpoints(inliers, normal);
                    p_s_new = endpoints.front();
                    p_s_temp = p_s_new + ps_new;
                    extend_cnt++;
                }
                else {
                    //p_s_new = cv::Point2f(0, 0);
                }
                // p_s_temp = p_s_new + ps_new;
                // extend_cnt++;
    
            }

        }
        
    }

    if (extended_lines.back().size() >= 4) {
        cv::Point2f pe_new(extended_lines.back().back()[0], extended_lines.back().back()[1]);
        std::vector<std::vector<double>> extended_frame_e = CM(evs, pe_new, flow_init, wsize);
        if (!extended_frame_e.empty()) {
            std::vector<cv::Point2f> H_pts = sort_points(extended_frame_e, wsize);

            if (H_pts.size() >= 10) {
                ransac_pair ransac_result = ransac_line(H_pts, max_iter, 1);

                V3d normal = ransac_result.first;
                std::vector<cv::Point2f> inliers = ransac_result.second;

                double inlier_ratio = static_cast<double>(ransac_result.third) / H_pts.size();
                //ROS_WARN("inlier_ratio: %f ", inlier_ratio);

                cv::Point2f v1(normal_[0], normal_[1]);
                cv::Point2f v2(normal[0], normal[1]);

                double temp = (v1.x * v2.x + v1.y * v2.y) / (cv::norm(v1) * cv::norm(v2));
                if (temp > 1.0) {
                                temp = 1.0;
                }
                if (temp < -1.0) {
                                temp = -1.0;
                }
                double theta = R2D * acos(temp);
                double min_theta = std::min(theta, 180 - theta);

                cv::Point2f p_e_new;
                if (inlier_ratio > 0.7 && min_theta < 10) {

                    std::vector<cv::Point2f> endpoints = proj_endpoints(inliers, normal);
                    p_e_new = endpoints.back();
                    p_e_temp = p_e_new + pe_new;
                    extend_cnt++;
                }
                else {
                    //p_e_new = cv::Point2f(0, 0);
                }
                // p_e_temp = p_e_new + pe_new;
                // extend_cnt++;

            }

        }
        
    }

    std::tuple<cv::Point2f, cv::Point2f, int> temp_tuple = std::make_tuple(p_s_temp, p_e_temp, extend_cnt);

    return temp_tuple;

}

std::vector<std::vector<V3d>> mTracker::extend_and_intensity_single(const cv::Mat &img_raw, const cv::Point2f& p_s, const cv::Point2f& p_e, const double& intensity_threshold) {
    cv::Point2f direction_ps = (p_s - p_e)/cv::norm(p_s - p_e);
    cv::Point2f direction_pe = (p_e - p_s)/cv::norm(p_e - p_s);

    std::vector<V3d> extended_line_ps = extend_line_single(img_raw, p_s, direction_ps, intensity_threshold);
    std::vector<V3d> extended_line_pe = extend_line_single(img_raw, p_e, direction_pe, intensity_threshold);


    std::vector<std::vector<V3d>> res;
    res.push_back(extended_line_ps);    
    res.push_back(extended_line_pe);  

    return res;
}

std::vector<V3d> mTracker::extend_line_single(const cv::Mat &img_raw, const cv::Point2f& point, const cv::Point2f& direction, const double& intensity_threshold) {
    cv::Point2f pointRounded(cvRound(point.x), cvRound(point.y));
    double intensity_init = get_intensity_mean(img_raw, pointRounded, direction);
    std::vector<V3d> extended_line;
    extended_line.push_back(V3d(point.x, point.y, intensity_init));

    cv::Point2f current_point = point;

    while (true) {
        cv::Point2f next_point = current_point + direction;

        if (is_outside_image(next_point, img_raw)) {
            break;
        }

        cv::Point2f pointRounded(cvRound(next_point.x), cvRound(next_point.y));
        double current_intensity = get_intensity_mean(img_raw, pointRounded, direction);

        if (std::abs(current_intensity - extended_line.back()[2]) > intensity_threshold) {
            break;
        }

        if (current_intensity == 0) {
            break;
        }

        extended_line.push_back(V3d(next_point.x, next_point.y, current_intensity));

        current_point = next_point;


    }

    return extended_line;

}

double mTracker::get_intensity_mean(const cv::Mat &img_raw, const cv::Point2f& pointRounded, const cv::Point2f& direction) {
    double in_1, in_2, in_3, in_4, in_5;
    if (std::abs(direction.x/direction.y) > 1) {
        in_1 = get_intensity(img_raw, pointRounded);
        in_2 = get_intensity(img_raw, pointRounded + cv::Point2f(0, 1));
        in_3 = get_intensity(img_raw, pointRounded + cv::Point2f(0, -1));
        in_4 = get_intensity(img_raw, pointRounded + cv::Point2f(0, 2));
        in_5 = get_intensity(img_raw, pointRounded + cv::Point2f(0, -2));
    }
    else {
        in_1 = get_intensity(img_raw, pointRounded);
        in_2 = get_intensity(img_raw, pointRounded + cv::Point2f(1, 0));
        in_3 = get_intensity(img_raw, pointRounded + cv::Point2f(-1, 0));
        in_4 = get_intensity(img_raw, pointRounded + cv::Point2f(2, 0));
        in_5 = get_intensity(img_raw, pointRounded + cv::Point2f(-2, 0));
    }
    double intensity = (in_1 + in_2 + in_3 + in_4 + in_5)/5;

    return intensity;

}

double mTracker::get_intensity(const cv::Mat &img_raw, const cv::Point2f& pixel) {
    int x = cvRound(pixel.x);
    int y = cvRound(pixel.y);

    int height = img_raw.rows;
    int width = img_raw.cols;

    double intensity = 0;

    if (x >= 0 && x < width && y >= 0 && y < height) {
        intensity = static_cast<double>(img_raw.at<uchar>(y, x));
        //ROS_WARN("inten: %f",intensity);
        intensity = intensity/255;
    }

    return intensity;

}

bool mTracker::is_outside_image(const cv::Point2f& point, const cv::Mat &img_raw) {
    int x = cvRound(point.x);
    int y = cvRound(point.y);

    int height = img_raw.rows;
    int width = img_raw.cols;

    bool is_outside = false;

    if (x < 0 || x >= width || y < 0 || y >= height) {
        is_outside = true;
    }

    return is_outside;
    
}


V4d mTracker::IMU_propa_q(const std::vector<sensor_msgs::ImuConstPtr>& IMUs, const double& before_t, const V3d& ba, const V3d& bg) {
    // Attitude update
    double last_imu_t = before_t;
    V4d q0(1, 0, 0, 0);
    V4d q_12;
    for (const auto &imu_msg : IMUs) {
        double now_imu_t = imu_msg->header.stamp.toSec();
        double imu_dt = now_imu_t - last_imu_t;

        V3d wib_b(imu_msg->angular_velocity.x, imu_msg->angular_velocity.y, imu_msg->angular_velocity.z);
        wib_b -= bg;

        V3d theta = wib_b*imu_dt;
        V4d del_quat = nesl::rvec2quat(theta);
        q_12 = nesl::quatMultiply(q0, del_quat);

        q0 = q_12;
        last_imu_t = now_imu_t;
    }

    return q_12;        // R_12

}

void mTracker::LoadParameters(ros::NodeHandle &n){
    // fu, fv, cu, cv
    n.getParam("/mono_tracker/cam0/intrinsics", cam_intrinsics_);

    // radtan: k1, k2, p1, p2, equidistant: k1, k2, k3, k4
    n.getParam("/mono_tracker/cam0/distortion_coeffs", cam_dcf_);
    std::vector<int> resolution(2); // width, height
    n.getParam("/mono_tracker/cam0/resolution", resolution);
    w_img_ = resolution[0];
    h_img_ = resolution[1];

    n.getParam("/mono_tracker/max_events", MAX_EVENTS_);
    n.getParam("/mono_tracker/life_pixel", LIFE_PIXEL_);
    n.getParam("/mono_tracker/max_cnt", MAX_CNT_);
    n.getParam("/mono_tracker/min_dist", MIN_DIST_);
    n.getParam("/mono_tracker/window_size", wsize_);
    n.getParam("/mono_tracker/cam0/distortion_model", DISTORTION_MODEL_);

    // line detection parameters
    n.getParam("/mono_tracker/eta", ETA_);
    int temp = cvRound(h_img_*ETA_);
    if (temp % 2 == 0) {
        PATCH_SIZE_ = temp + 1;
    }
    else {
        PATCH_SIZE_ = temp;
    }
    ROS_WARN("patch size : %d", PATCH_SIZE_);
    //n.getParam("/mono_tracker/patch_size", PATCH_SIZE_);
    n.getParam("/mono_tracker/pix_th", PIX_THR_);
    n.getParam("/mono_tracker/ransac_iter", RANSAC_ITER_);

    // line matching parameters
    n.getParam("/mono_tracker/dist_thr", DIST_THR_);
    n.getParam("/mono_tracker/theta_thr", THETA_THR_);


    // R^{vehicle}_{IMU}
    double roll_iv, pitch_iv, yaw_iv;
    n.param<double>("estimator_node/roll_imu_vehicle", roll_iv, 0.0);
    n.param<double>("estimator_node/pitch_imu_vehicle", pitch_iv, 0.0);
    n.param<double>("estimator_node/yaw_imu_vehicle", yaw_iv, 0.0);
    V3d euler_iv(D2R*roll_iv, D2R*pitch_iv, D2R*yaw_iv);
    M3d Rvi_ = nesl::euler2dcm(euler_iv);

    // read from yaml file
    // Cam/IMU extrinsics
    Eigen::Isometry3d T_cam_imu =
        msckf_vio::utils::getTransformEigen(n, "/estimator_node/cam0/T_cam_imu");
    M4d Tcb = M4d::Zero();
    Tcb.block<3,3>(0,0) = T_cam_imu.linear() * Rvi_.transpose();;
    Tcb.block<3,1>(0,3) = T_cam_imu.translation();
    Tcb(3,3) = 1.0;
    Tbc_ = Tcb.inverse();


    ROS_INFO("Frontend parameter loading finished !");
}


void mTracker::readMonoImage(const cv::Mat &img_raw, double cur_time){

    cur_time_ = cur_time;
    cv::Mat img; // current left image

    // Histogram equalization
    // using Constrast Limited Adaptive Histogram Equalization
    cv::Ptr<cv::CLAHE> clahe = cv::createCLAHE(3.0, cv::Size(8,8));
    clahe->apply(img_raw, img);

    if (cur_img_.empty()) {
        cur_img_ = prev_img_ = img;
    }
    else {
        cur_img_ = img;
    }

    cur_pts_.clear();

    

    // Temporal & static tracker
    if (prev_pts_.size() > 0){
        // Temporal feature tracker based on left images
        std::vector<uchar> status;
        std::vector<float> err;
        cv::calcOpticalFlowPyrLK(prev_img_, cur_img_,
                                 prev_pts_, cur_pts_,
                                 status, err, cv::Size(21,21), 3);

        for (int i = 0; i < int(cur_pts_.size()); i++) {
            if (status[i] && !inBorder(cur_pts_[i], w_img_, h_img_)) {
                status[i] = 0;
            }
        }

        reduceVector(prev_pts_, status);
        reduceVector(cur_pts_, status);
        reduceVector(ids_, status);
        reduceVector(track_cnt_, status);
        rejectWithF();

    }
    for (auto &n : track_cnt_)
        n++;

    setMask();

    

    // extract new features if necessary
    int n_max_cnt = MAX_CNT_ - static_cast<int>(cur_pts_.size());
    if (n_max_cnt > 0){
        if(mask_.empty())
            std::cout << "mask is empty " << std::endl;
        if (mask_.type() != CV_8UC1)
            std::cout << "mask type is wrong " << std::endl;
        if (mask_.size() != cur_img_.size())
            std::cout << "mask size is wrong " << std::endl;
        // cv::goodFeaturesToTrack(cur_img0, n_pts0,
        //                         MAX_CNT - cur_pts0.size(),
        //                         0.01, MIN_DIST, mask);
        cv::goodFeaturesToTrack(cur_img_, n_pts_,
                                MAX_CNT_ - cur_pts_.size(),
                                0.001, MIN_DIST_, mask_);
    }
    else {
        n_pts_.clear();
    }
    addPoints();

    // update member vars
    prev_img_ = cur_img_;
    prev_pts_ = cur_pts_;
    undistortedPoints();
    prev_time_ = cur_time_;

    
}



void mTracker::readEventImage(const cv::Mat &img_raw, double cur_time, const std::vector<dvs_msgs::Event>& evs, const V4d& q_12){

    ros::Time start_time = ros::Time::now();

    cur_time_e = cur_time;
    cv::Mat img; // current left image

    // Histogram equalization
    // using Constrast Limited Adaptive Histogram Equalization
    cv::Ptr<cv::CLAHE> clahe = cv::createCLAHE(3.0, cv::Size(8,8));
    clahe->apply(img_raw, img);
    //img = img_raw;

    if (cur_img_e.empty()) {
        cur_img_e = prev_img_e = img;
    }
    else {
        cur_img_e = img;
    }

    cur_pts_e.clear();
    cur_id_point_line.clear();
    final_lines.clear();

    if (several_windows) {
        prev_pts_e = after_pts_e;
    }

    // Temporal & static tracker
    if (prev_pts_e.size() > 0){

        // Temporal feature tracker based on left images
        std::vector<uchar> status;
        std::vector<float> err;
        cv::calcOpticalFlowPyrLK(prev_img_e, cur_img_e,
                                 prev_pts_e, cur_pts_e,
                                 status, err, cv::Size(21,21), 3);

        for (int i = 0; i < int(cur_pts_e.size()); i++) {
            if (status[i] && !inBorder(cur_pts_e[i], w_img_, h_img_)) {
                status[i] = 0;
            }
        }

        reduceVector(prev_pts_e, status);
        reduceVector(cur_pts_e, status);
        reduceVector(ids_e, status);
        reduceVector(track_cnt_e, status);
        rejectWithF_e();

    }

    for (auto &n : track_cnt_e)
        n++;

    double eventpoint_processing_time = (ros::Time::now()-start_time).toSec();
    eventpoint_sum += eventpoint_processing_time;

    // Event line extract
    if (prev_pts_e.size() > 0) {
        double dt = cur_time_e - prev_time_e;
        std::vector<cv::Point2f> finalCfts = cur_pts_e;
        std::vector<cv::Point2f> finalPfts = prev_pts_e;
        std::vector<V4d> linePts(finalCfts.size(), V4d(0, 0, 0, 0));

        ros::Time start_time = ros::Time::now();

        #pragma omp parallel for
        for (int i = 0; i < finalCfts.size(); i++) {
            int thread_num = omp_get_thread_num();
            cv::Point2f flow_init = (finalCfts[i] - finalPfts[i]) / dt;

            // Fine motion compensation (Contrast Maximization)
            std::vector<std::vector<double>> sharp_frame = CM(evs, finalCfts[i], flow_init, PATCH_SIZE_);

            if (sharp_frame.empty()) continue;

            // Informative points selection
            std::vector<cv::Point2f> H_pts = sort_points(sharp_frame, PATCH_SIZE_);

            if (H_pts.size() >= 10) {

                // RANSAC line fitting
                ransac_pair ransac_result = ransac_line(H_pts, RANSAC_ITER_, 1);

                V3d normal = ransac_result.first;
                std::vector<cv::Point2f> inliers = ransac_result.second;

                double inlier_ratio = static_cast<double>(ransac_result.third) / H_pts.size();

                if (inlier_ratio > 0.4) {       

                    std::vector<cv::Point2f> endpoints = proj_endpoints(inliers, normal);

                    cv::Point2f p_s = endpoints.front();
                    cv::Point2f p_e = endpoints.back();

                    p_s += finalCfts[i];
                    p_e += finalCfts[i];

                    linePts[i] = V4d(p_s.x, p_s.y, p_e.x, p_e.y);
                }                   
            }        
        }
        double line_extract_time = (ros::Time::now()-start_time).toSec();
        //ROS_WARN("line extraction t: %f", line_extract_time);


        cnt_line_ex++;
        line_extract_sum += line_extract_time;
        mean_line_extact = (double)line_extract_sum/cnt_line_ex;

        ROS_WARN("Fine M.C. & event-line detection: %f [ms]", mean_line_extact*1000);


    
        cur_lines_e = linePts;

        setMask_e();

        // ID & point & line set
        std::vector<std::tuple<int, cv::Point2f, V4d>> id_point_line;
        id_point_line.reserve(cur_pts_e.size()); 

        for (unsigned int i = 0; i < cur_pts_e.size(); i++) {
            id_point_line.emplace_back(ids_e[i], cur_pts_e[i], cur_lines_e[i]); 
        }

        cur_id_point_line = id_point_line;

        start_time = ros::Time::now();

        // Event line matching
        if (!prev_id_point_line.empty()) {
            M3d Kmtx = getK();
            M3d K_inv = Kmtx.inverse();
            M3d R_b2b1 = nesl::quat2dcm(q_12).transpose();
            M3d R_bc = Tbc_.block<3,3>(0,0);
            M3d R_c2c1 = R_bc.transpose()*R_b2b1*R_bc;

            std::vector<std::tuple<int, cv::Point2f, cv::Point2f, V4d, V4d>> synced_lines = synchronizeLines(cur_id_point_line, prev_id_point_line);

            # pragma omp parallel for
            for (int i = 0; i < synced_lines.size(); i++) {
                int id = std::get<0>(synced_lines[i]);
                cv::Point2f p_ft = std::get<1>(synced_lines[i]);
                cv::Point2f c_ft = std::get<2>(synced_lines[i]);
                V4d p_line = std::get<3>(synced_lines[i]);
                V4d c_line = std::get<4>(synced_lines[i]);

                if (isNonZero(p_line) && isNonZero(c_line)) {
                    cv::Point2f cur_ps = cv::Point2f(c_line[0], c_line[1]) - c_ft;
                    cv::Point2f cur_pe = cv::Point2f(c_line[2], c_line[3]) - c_ft;
                    cv::Point2f prev_ps = cv::Point2f(p_line[0], p_line[1]) - p_ft;
                    cv::Point2f prev_pe = cv::Point2f(p_line[2], p_line[3]) - p_ft;

                    //ROS_WARN("cfts: (%f, %f), pfts: (%f, %f)", c_ft.x, c_ft.y, p_ft.x, p_ft.y);

                    cv::Point2f mid_point = (cur_ps + cur_pe)/2;
                    cv::Point2f normal_p = cv::Point2f(-(prev_pe.y - prev_ps.y), (prev_pe.x - prev_ps.x));
                    double c = prev_ps.x*prev_pe.y - prev_pe.x*prev_ps.y;
                    double dist = std::abs(normal_p.x*mid_point.x + normal_p.y*mid_point.y + c)/cv::norm(normal_p);

                    V3d ps_temp(p_line[0], p_line[1], 1);
                    V3d pe_temp(p_line[2], p_line[3], 1);
                    V3d ps_c2 = R_c2c1*Kmtx.inverse()*ps_temp;
                    ps_c2 /= ps_c2[2];
                    ps_c2 = Kmtx*ps_c2;
                    V3d pe_c2 = R_c2c1*Kmtx.inverse()*pe_temp;
                    pe_c2 /= pe_c2[2];
                    pe_c2 = Kmtx*pe_c2;

                    cv::Point2f v1 = cv::Point2f(c_line[0], c_line[1]) - cv::Point2f(c_line[2], c_line[3]);
                    cv::Point2f v2 = cv::Point2f(p_line[0], p_line[1]) - cv::Point2f(p_line[2], p_line[3]);

                    double temp = (v1.x*v2.x + v1.y*v2.y)/(cv::norm(v1)*cv::norm(v2));
                    if (temp > 1.0) {
                        temp = 1.0;
                    }
                    if (temp < -1.0) {
                        temp = -1.0;
                    }
                    double theta = R2D*acos(temp);
                    double min_theta = std::min(theta, 180-theta);

                    //ROS_WARN("dist: %f, theta: %f", dist, min_theta);

                    if (dist < DIST_THR_ && min_theta < THETA_THR_) {
                        #pragma omp critical                        // access only one thread at once
                        final_lines.emplace_back(id, c_line);
                    }
                }
            }

        }

        else {
            for (int i = 0; i < cur_id_point_line.size(); i++) {
                int id = std::get<0>(cur_id_point_line[i]);
                V4d line = std::get<2>(cur_id_point_line[i]);
                if (isNonZero(line)) {
                    final_lines.emplace_back(id, line);
                }
            }
        }
        double line_matching_time = (ros::Time::now()-start_time).toSec();
        //ROS_WARN("line matching t: %f", line_matching_time);

        cnt_line_mat++;
        line_matching_sum += line_matching_time;
        mean_line_matching = (double)line_matching_sum/cnt_line_mat;

        ROS_WARN("Event-line tracking: %f [ms]", mean_line_matching*1000);

        prev_id_point_line = cur_id_point_line;
    }

    else {
        setMask_e_point();
    }

    start_time = ros::Time::now();

    // extract new features if necessary
    int n_max_cnt = MAX_CNT_ - static_cast<int>(cur_pts_e.size());

    if (n_max_cnt > 0){
        if(mask_e.empty())
            std::cout << "mask is empty " << std::endl;
        if (mask_e.type() != CV_8UC1)
            std::cout << "mask type is wrong " << std::endl;
        if (mask_e.size() != cur_img_e.size())
            std::cout << "mask size is wrong " << std::endl;
        // cv::goodFeaturesToTrack(cur_img0, n_pts0,
        //                         MAX_CNT - cur_pts0.size(),
        //                         0.01, MIN_DIST, mask);
        cv::goodFeaturesToTrack(cur_img_e, n_pts_e,
                                MAX_CNT_ - cur_pts_e.size(),
                                0.001, MIN_DIST_, mask_e);
    }
    else {
        n_pts_e.clear();
    }
    addPoints_e();

    // update member vars
    prev_img_e = cur_img_e;
    prev_pts_e = cur_pts_e;
    before_pts_e = cur_pts_e;
    undistortedPoints_e();
    undistortedLines_e();       // undistor and K project
    prev_time_e = cur_time_e;
    prev_time_fix = cur_time_e;  // for velociy calculation


    eventpoint_processing_time = (ros::Time::now()-start_time).toSec();
    eventpoint_sum += eventpoint_processing_time;

    cnt_eventpoint++;
    mean_eventpointprocess = (double)eventpoint_sum/cnt_eventpoint;

    ROS_WARN("Event-point detection & tracking: %f [ms]", mean_eventpointprocess*1000);

}


void mTracker::readEventImage_no_cnt(const cv::Mat &img_raw, double cur_time){
    cur_time_e = cur_time;
    cv::Mat img; // current left image

    // Histogram equalization
    // using Constrast Limited Adaptive Histogram Equalization

    cv::Ptr<cv::CLAHE> clahe = cv::createCLAHE(3.0, cv::Size(8,8));
    clahe->apply(img_raw, img);
    //img = img_raw;

    if (cur_img_e.empty()) {
        cur_img_e = prev_img_e = img;
    }
    else {
        cur_img_e = img;
    }

    // Temporal & static tracker
    if (before_pts_e.size() > 0){
        // Temporal feature tracker based on left images
        std::vector<uchar> status;
        std::vector<float> err;
        cv::calcOpticalFlowPyrLK(prev_img_e, cur_img_e,
                                 before_pts_e, after_pts_e,
                                 status, err, cv::Size(21,21), 3);

        for (int i = 0; i < int(after_pts_e.size()); i++) {
            if (status[i] && !inBorder(after_pts_e[i], w_img_, h_img_)) {
                status[i] = 0;
            }
        }

        reduceVector(before_pts_e, status);
        reduceVector(after_pts_e, status);
        reduceVector(ids_e, status);
        reduceVector(track_cnt_e, status);
        rejectWithF_e_no_cnt();
    }

    setMask_e_no_cnt();

    // extract new features if necessary
    int n_max_cnt = MAX_CNT_ - static_cast<int>(after_pts_e.size());

    if (n_max_cnt > 0){
        if(mask_e.empty())
            std::cout << "mask is empty " << std::endl;
        if (mask_e.type() != CV_8UC1)
            std::cout << "mask type is wrong " << std::endl;
        if (mask_e.size() != cur_img_e.size())
            std::cout << "mask size is wrong " << std::endl;
        // cv::goodFeaturesToTrack(cur_img0, n_pts0,
        //                         MAX_CNT - cur_pts0.size(),
        //                         0.01, MIN_DIST, mask);
        cv::goodFeaturesToTrack(cur_img_e, n_pts_e,
                                MAX_CNT_ - after_pts_e.size(),
                                0.001, MIN_DIST_, mask_e);
    }
    else {
        n_pts_e.clear();
    }
    addPoints_e_no_cnt();

    // update member vars
    prev_img_e = cur_img_e;
    before_pts_e = after_pts_e;
    prev_time_e = cur_time_e;
}


cv::Mat mTracker::makeEventFrame(const std::vector<dvs_msgs::Event>& evs, const std::vector<sensor_msgs::ImuConstPtr>& imus, const V3d& ba, const V3d& bg) {

    ros::Time start_time = ros::Time::now();

    M3d Kmtx = getK();
    V3d f_mean = V3d::Zero();
    V3d w_mean = V3d::Zero();

    for (unsigned int i = 0; i < imus.size(); i++) {
        V3d f_i(imus[i]->linear_acceleration.x,
            imus[i]->linear_acceleration.y,
            imus[i]->linear_acceleration.z);
        V3d w_i(imus[i]->angular_velocity.x,
            imus[i]->angular_velocity.y,
            imus[i]->angular_velocity.z);
        f_mean += f_i;
        w_mean += w_i; 
    }
    f_mean /= imus.size();
    w_mean /= imus.size();

    f_mean -= ba;
    w_mean -= bg;

    //ROS_INFO("w_mean: %f, %f, %f", w_mean[0], w_mean[1], w_mean[2]);

    //ROS_WARN("ba: %f, %f, %f", ba[0], ba[1], ba[2]);
    //ROS_WARN("bg: %f, %f, %f", bg[0], bg[1], bg[2]);

    double t0 = evs.back().ts.toSec();
    std::vector<cv::Point2f> warped_events;

    for (const auto& ev : evs) {
        double del_t = ev.ts.toSec() - t0;
        M3d R_bref_bi = nesl::rvec2dcm(w_mean * del_t);
        M3d R_bc = Tbc_.block<3,3>(0,0);
        M3d R_cref_ci = R_bc.transpose() * R_bref_bi * R_bc;

        V3d before_xy1(ev.x, ev.y, 1);
        V3d temp = R_cref_ci * Kmtx.inverse() * before_xy1;
        V3d temp_(temp[0]/temp[2], temp[1]/temp[2], temp[2]/temp[2]);

        V3d warped_xyz = Kmtx * temp_;

        double x_ = warped_xyz[0];
        double y_ = warped_xyz[1];

        if (x_ < w_img_ - 1 && x_ > 0 && y_ < h_img_ - 1 && y_ > 0) {
            warped_events.emplace_back(x_, y_);
        }
    }

    cv::Mat eventFrame_;
    eventFrame_ = cv::Mat(h_img_, w_img_, CV_8UC1);
    std::vector<std::vector<int>> H(h_img_, std::vector<int>(w_img_, 0));

    for (unsigned int i = 0; i < warped_events.size(); i++) {
        int img_x = cvRound(warped_events[i].x);
        int img_y = cvRound(warped_events[i].y);
        H[img_y][img_x] += 1;
    }
    int max_ = 0;
    for (const auto& row : H) {
        for (int value : row) {
            max_ = std::max(max_, value);
        }
    }
    for (auto& row : H) {
        for (int& value : row) {
            value = static_cast<int>((static_cast<double>(value) * 255) / max_);
        }
    }
    for (int i = 0; i < h_img_; ++i) {
        for (int j = 0; j < w_img_; ++j) {
            // Assign the pixel value from the matrix H to the corresponding position in the eventFrame
            eventFrame_.at<uchar>(i, j) = static_cast<uchar>(H[i][j]);
        }
    }

    //ROS_INFO("max of H: %d ", max_);
    //ROS_WARN("here");

    double eventframe_processing_time = (ros::Time::now()-start_time).toSec();
    eventframe_sum += eventframe_processing_time;

    cnt_eventframe++;
    mean_eventframe = (double)eventframe_sum/cnt_eventframe;

    ROS_WARN("Coarse M.C.: %f [ms]", mean_eventframe*1000);

    return eventFrame_;

}


void mTracker::setMask(){
    mask_ = cv::Mat(h_img_, w_img_, CV_8UC1, cv::Scalar(255));

    // prefer to keep features that are tracked for long time
    // first : left pts, second : right pts
    std::vector<std::pair<int, std::pair<cv::Point2f, int>>> cnt_pts_id;
    for (unsigned int i = 0; i < cur_pts_.size(); i ++){
        cnt_pts_id.push_back(std::make_pair(track_cnt_[i],
                                            std::make_pair(cur_pts_[i], ids_[i])));
    }
    std::sort(cnt_pts_id.begin(), cnt_pts_id.end(),
              [](const std::pair<int, std::pair<cv::Point2f, int>> &a,
                 const std::pair<int, std::pair<cv::Point2f, int>> &b)
    {
        return a.first > b.first;
    });

    cur_pts_.clear();
    ids_.clear();
    track_cnt_.clear();

    for (auto &it : cnt_pts_id){
        if (mask_.at<uchar>(it.second.first) == 255){
            cur_pts_.push_back(it.second.first);
            ids_.push_back(it.second.second);
            track_cnt_.push_back(it.first);
            // filled circle with MIN_DIST
            cv::circle(mask_, it.second.first, MIN_DIST_, 0, -1);
        }
    }

}

void mTracker::setMask_e(){
    mask_e = cv::Mat(h_img_, w_img_, CV_8UC1, cv::Scalar(255));

    // Prefer to keep features that are tracked for a long time
    std::vector<std::tuple<int, cv::Point2f, V4d, int>> cnt_pts_id;
    for (size_t i = 0; i < cur_pts_e.size(); ++i) {
        cnt_pts_id.emplace_back(track_cnt_e[i], cur_pts_e[i], cur_lines_e[i], ids_e[i]);
    }

    // Sort the features based on the tracking count (in descending order)
    std::sort(cnt_pts_id.begin(), cnt_pts_id.end(),
              [](const std::tuple<int, cv::Point2f, V4d, int> &a,
                 const std::tuple<int, cv::Point2f, V4d, int> &b) {
        return std::get<0>(a) > std::get<0>(b);
    });

    cur_pts_e.clear();
    cur_lines_e.clear();
    ids_e.clear();
    track_cnt_e.clear();

    for (const auto &it : cnt_pts_id) {
        cv::Point2f pt = std::get<1>(it);
        if (mask_e.at<uchar>(pt) == 255) {
            cur_pts_e.push_back(pt);
            cur_lines_e.push_back(std::get<2>(it));
            ids_e.push_back(std::get<3>(it));
            track_cnt_e.push_back(std::get<0>(it));
            // Filled circle with MIN_DIST_
            cv::circle(mask_e, pt, MIN_DIST_, 0, -1);
        }
    }

}


void mTracker::setMask_e_point() {
    mask_e = cv::Mat(h_img_, w_img_, CV_8UC1, cv::Scalar(255));

    // prefer to keep features that are tracked for long time
    // first : left pts, second : right pts
    std::vector<std::pair<int, std::pair<cv::Point2f, int>>> cnt_pts_id;
    for (unsigned int i = 0; i < cur_pts_e.size(); i ++){
        cnt_pts_id.push_back(std::make_pair(track_cnt_e[i],
                                            std::make_pair(cur_pts_e[i], ids_e[i])));
    }
    std::sort(cnt_pts_id.begin(), cnt_pts_id.end(),
              [](const std::pair<int, std::pair<cv::Point2f, int>> &a,
                 const std::pair<int, std::pair<cv::Point2f, int>> &b)
    {
        return a.first > b.first;
    });

    cur_pts_e.clear();
    ids_e.clear();
    track_cnt_e.clear();

    for (auto &it : cnt_pts_id){
        if (mask_e.at<uchar>(it.second.first) == 255){
            cur_pts_e.push_back(it.second.first);
            ids_e.push_back(it.second.second);
            track_cnt_e.push_back(it.first);
            // filled circle with MIN_DIST
            cv::circle(mask_e, it.second.first, MIN_DIST_, 0, -1);
        }
    }
}


void mTracker::setMask_e_no_cnt(){
    mask_e = cv::Mat(h_img_, w_img_, CV_8UC1, cv::Scalar(255));

    // prefer to keep features that are tracked for long time
    // first : left pts, second : right pts
    std::vector<std::pair<int, std::pair<cv::Point2f, int>>> cnt_pts_id;
    for (unsigned int i = 0; i < after_pts_e.size(); i ++){
        cnt_pts_id.push_back(std::make_pair(track_cnt_e[i],
                                            std::make_pair(after_pts_e[i], ids_e[i])));
    }
    std::sort(cnt_pts_id.begin(), cnt_pts_id.end(),
              [](const std::pair<int, std::pair<cv::Point2f, int>> &a,
                 const std::pair<int, std::pair<cv::Point2f, int>> &b)
    {
        return a.first > b.first;
    });

    after_pts_e.clear();
    ids_e.clear();
    track_cnt_e.clear();

    for (auto &it : cnt_pts_id){
        if (mask_e.at<uchar>(it.second.first) == 255){
            after_pts_e.push_back(it.second.first);
            ids_e.push_back(it.second.second);
            track_cnt_e.push_back(it.first);
            // filled circle with MIN_DIST
            cv::circle(mask_e, it.second.first, MIN_DIST_, 0, -1);
        }
    }
}

void mTracker::addPoints(){
    for (auto &p : n_pts_){
        cur_pts_.push_back(p);
        ids_.push_back(-1);
        track_cnt_.push_back(1);
    }
}

void mTracker::addPoints_e(){
    for (auto &p : n_pts_e){
        cur_pts_e.push_back(p);
        ids_e.push_back(-1);
        track_cnt_e.push_back(1);
    }
}

void mTracker::addPoints_e_no_cnt(){
    for (auto &p : n_pts_e){
        after_pts_e.push_back(p);
        ids_e.push_back(-1);
        track_cnt_e.push_back(0);
    }
}

bool mTracker::updateID(unsigned int i)
{
    if (i < ids_.size())
    {
        if (ids_[i] == -1)
            ids_[i] = n_id++;
        return true;
    }
    else
        return false;
}

bool mTracker::updateID_e(unsigned int i)
{
    if (i < ids_e.size())
    {
        if (ids_e[i] == -1)
            ids_e[i] = n_id_e++;
        return true;
    }
    else
        return false;
}

void mTracker::rejectWithF(){
    if (cur_pts_.size() >= 8){
        std::vector<cv::Point2f> prev_pts, cur_pts;
        double fu, fv, cu, cv;

        // left image ransac
        prev_pts = prev_pts_;
        cur_pts = cur_pts_;
        fu = cam_intrinsics_[0];
        fv = cam_intrinsics_[1];
        cu = cam_intrinsics_[2];
        cv = cam_intrinsics_[3];
        

        std::vector<cv::Point2f> un_prev_pts(prev_pts.size());
        std::vector<cv::Point2f> un_cur_pts(cur_pts.size());
        for (unsigned int i = 0; i < prev_pts.size(); i++){
            Eigen::Vector3d tmp_p;
            liftProjective(Eigen::Vector2d(prev_pts[i].x, prev_pts[i].y),
                           tmp_p);
            tmp_p.x() = fu * (tmp_p.x() / tmp_p.z()) + cu;
            tmp_p.y() = fv * (tmp_p.y() / tmp_p.z()) + cv;
            un_prev_pts[i] = cv::Point2f(tmp_p.x(), tmp_p.y());

            liftProjective(Eigen::Vector2d(cur_pts[i].x, cur_pts[i].y),
                           tmp_p);
            tmp_p.x() = fu * (tmp_p.x() / tmp_p.z()) + cu;
            tmp_p.y() = fv * (tmp_p.y() / tmp_p.z()) + cv;
            un_cur_pts[i] = cv::Point2f(tmp_p.x(), tmp_p.y());
        }
        std::vector<uchar> status;
        cv::findFundamentalMat(un_prev_pts, un_cur_pts,
                               cv::FM_RANSAC, 1.0, 0.99, status);
        reduceVector(prev_pts_, status);
        reduceVector(cur_pts_, status);
        reduceVector(ids_, status);
        reduceVector(track_cnt_, status);

    }
}

void mTracker::rejectWithF_e(){
    if (cur_pts_e.size() >= 8){
        std::vector<cv::Point2f> prev_pts, cur_pts;
        double fu, fv, cu, cv;

        // left image ransac
        prev_pts = prev_pts_e;
        cur_pts = cur_pts_e;
        fu = cam_intrinsics_[0];
        fv = cam_intrinsics_[1];
        cu = cam_intrinsics_[2];
        cv = cam_intrinsics_[3];
        

        std::vector<cv::Point2f> un_prev_pts(prev_pts.size());
        std::vector<cv::Point2f> un_cur_pts(cur_pts.size());
        for (unsigned int i = 0; i < prev_pts.size(); i++){
            Eigen::Vector3d tmp_p;
            liftProjective(Eigen::Vector2d(prev_pts[i].x, prev_pts[i].y),
                           tmp_p);
            tmp_p.x() = fu * (tmp_p.x() / tmp_p.z()) + cu;
            tmp_p.y() = fv * (tmp_p.y() / tmp_p.z()) + cv;
            un_prev_pts[i] = cv::Point2f(tmp_p.x(), tmp_p.y());

            liftProjective(Eigen::Vector2d(cur_pts[i].x, cur_pts[i].y),
                           tmp_p);
            tmp_p.x() = fu * (tmp_p.x() / tmp_p.z()) + cu;
            tmp_p.y() = fv * (tmp_p.y() / tmp_p.z()) + cv;
            un_cur_pts[i] = cv::Point2f(tmp_p.x(), tmp_p.y());
        }
        std::vector<uchar> status;
        cv::findFundamentalMat(un_prev_pts, un_cur_pts,
                               cv::FM_RANSAC, 1.0, 0.99, status);
        reduceVector(prev_pts_e, status);
        reduceVector(cur_pts_e, status);
        reduceVector(ids_e, status);
        reduceVector(track_cnt_e, status);
    }
}

void mTracker::rejectWithF_e_no_cnt(){
    if (after_pts_e.size() >= 8){
        std::vector<cv::Point2f> prev_pts, cur_pts;
        double fu, fv, cu, cv;

        // left image ransac
        prev_pts = before_pts_e;
        cur_pts = after_pts_e;
        fu = cam_intrinsics_[0];
        fv = cam_intrinsics_[1];
        cu = cam_intrinsics_[2];
        cv = cam_intrinsics_[3];
        

        std::vector<cv::Point2f> un_prev_pts(prev_pts.size());
        std::vector<cv::Point2f> un_cur_pts(cur_pts.size());
        for (unsigned int i = 0; i < prev_pts.size(); i++){
            Eigen::Vector3d tmp_p;
            liftProjective(Eigen::Vector2d(prev_pts[i].x, prev_pts[i].y),
                           tmp_p);
            tmp_p.x() = fu * (tmp_p.x() / tmp_p.z()) + cu;
            tmp_p.y() = fv * (tmp_p.y() / tmp_p.z()) + cv;
            un_prev_pts[i] = cv::Point2f(tmp_p.x(), tmp_p.y());

            liftProjective(Eigen::Vector2d(cur_pts[i].x, cur_pts[i].y),
                           tmp_p);
            tmp_p.x() = fu * (tmp_p.x() / tmp_p.z()) + cu;
            tmp_p.y() = fv * (tmp_p.y() / tmp_p.z()) + cv;
            un_cur_pts[i] = cv::Point2f(tmp_p.x(), tmp_p.y());
        }
        std::vector<uchar> status;
        cv::findFundamentalMat(un_prev_pts, un_cur_pts,
                               cv::FM_RANSAC, 1.0, 0.99, status);
        reduceVector(before_pts_e, status);
        reduceVector(after_pts_e, status);
        reduceVector(ids_e, status);
        reduceVector(track_cnt_e, status);
    }
}

void mTracker::liftProjective(const Eigen::Vector2d &p,
                              Eigen::Vector3d &P) const{
    // Intinrisc parameters
    double inv_K11, inv_K13, inv_K22, inv_K23;
    double fu, fv, cu, cv, k1, k2, k3, k4;

    inv_K11 = 1.0 / cam_intrinsics_[0];
    inv_K13 = -cam_intrinsics_[2] / cam_intrinsics_[0];
    inv_K22 = 1.0 / cam_intrinsics_[1];
    inv_K23 = -cam_intrinsics_[3] / cam_intrinsics_[1];
    k1 = cam_dcf_[0];
    k2 = cam_dcf_[1];
    k3 = cam_dcf_[2];
    k4 = cam_dcf_[3];
    fu = cam_intrinsics_[0];
    fv = cam_intrinsics_[1];
    cu = cam_intrinsics_[2];
    cv = cam_intrinsics_[3];
    

    double mx_d, my_d, mx_u, my_u;
    mx_d = inv_K11 * p(0) + inv_K13;
    my_d = inv_K22 * p(1) + inv_K23;

    if (DISTORTION_MODEL_.compare("radtan") == 0) {
        // Recursive distortion model
        int n = 8;
        Eigen::Vector2d d_u;
        distortion(Eigen::Vector2d(mx_d, my_d), d_u,
            Eigen::Vector4d(k1, k2, k3, k4));

        // Approximate value: 
        // subtract delta distortion (d_u) from current estimate
        mx_u = mx_d - d_u(0);
        my_u = my_d - d_u(1);
        for (int i = 1; i < n; ++i){
            distortion(Eigen::Vector2d(mx_u, my_u), d_u,
                Eigen::Vector4d(k1, k2, k3, k4));
            mx_u = mx_d - d_u(0);
            my_u = my_d - d_u(1);
        }
    } else if (DISTORTION_MODEL_.compare("equidistant") == 0) {
        // based on schneith's implementation
        // iteratively solve inverse equidistant model
        double d_d = sqrt(mx_d*mx_d + my_d*my_d);
        double r = d_d;
        double r2, r4, r6, r8, scaling;
        for (int i = 1; i < 20; ++i) {
            r2 = r*r;
            r4 = r2*r2;
            r6 = r4*r2;
            r8 = r4*r4;
            r = d_d / (1 + k1*r2 + k2*r4 + k3*r6 + k4*r8);
        }
        scaling = tan(r) / d_d;
        mx_u = mx_d*scaling;
        my_u = my_d*scaling;

        // // opencv implementation
        // cv::Mat mat(1, 2, CV_32F);
        // mat.at<float>(0, 0) = p(0);
        // mat.at<float>(0, 1) = p(1);
        // mat = mat.reshape(2); // Nx1, 1-channel
        // // Eigen mtx to cv mtx
        // cv::Matx33d camK;
        // camK(0, 0) = fu;
        // camK(0, 1) = 0;
        // camK(0, 2) = cu;
        // camK(1, 0) = 0;
        // camK(1, 1) = fv;
        // camK(1, 2) = cv;
        // camK(2, 0) = 0;
        // camK(2, 1) = 0;
        // camK(2, 2) = 1;
        // cv::Vec4d camD;
        // camD(0) = k1;
        // camD(1) = k2;
        // camD(2) = k3;
        // camD(3) = k4;
        // // Undistort it!
        // // in: pixel coordinate, out: normalized coordinate
        // cv::fisheye::undistortPoints(mat, mat, camK, camD);
        // // Constrcut our return vector
        // mx_u = mat.at<float>(0, 0);
        // my_u = mat.at<float>(0, 1);

    } else {
        ROS_ERROR("Please set appropriate distortion model !");
    }
    P << mx_u, my_u, 1.0;
}


void mTracker::distortion(const Eigen::Vector2d& p_u,
                          Eigen::Vector2d& d_u,
                          const Eigen::Vector4d& dist_coeff) const{
    double k1, k2, p1, p2;
    k1 = dist_coeff(0);
    k2 = dist_coeff(1);
    p1 = dist_coeff(2);
    p2 = dist_coeff(3);

    double mx2_u, my2_u, mxy_u, rho2_u, rad_dist_u;

    mx2_u = p_u(0) * p_u(0);
    my2_u = p_u(1) * p_u(1);
    mxy_u = p_u(0) * p_u(1);
    rho2_u = mx2_u + my2_u;
    rad_dist_u = k1 * rho2_u + k2 * rho2_u * rho2_u;
    d_u << p_u(0) * rad_dist_u + 2.0 * p1 * mxy_u + p2 * (rho2_u + 2.0 * mx2_u),
           p_u(1) * rad_dist_u + 2.0 * p2 * mxy_u + p1 * (rho2_u + 2.0 * my2_u);
}

cv::Point2f mTracker::undistort_a_point(const cv::Point2f& p) {

    Eigen::Vector2d a0(p.x, p.y);
    Eigen::Vector3d b0;
    liftProjective(a0, b0);
    V3d un_pt = getK()*b0;

    cv::Point2f un_p(un_pt[0], un_pt[1]);

    return un_p;
}

void mTracker::undistortedPoints(){
    cur_un_pts_.clear();
    cur_un_pts_map_.clear();
    for (unsigned int i = 0; i < cur_pts_.size(); i++)
    {
        Eigen::Vector2d a0(cur_pts_[i].x, cur_pts_[i].y);
        Eigen::Vector3d b0;
        liftProjective(a0, b0);
        cur_un_pts_.push_back(cv::Point2f(b0.x() / b0.z(), b0.y() / b0.z()));
        cur_un_pts_map_.insert(std::make_pair(ids_[i], cv::Point2f(b0.x() / b0.z(), b0.y() / b0.z())));

    }

    // caculate points velocity (left image)
    if (!prev_un_pts_map_.empty())
    {
        double dt = cur_time_ - prev_time_;
        pts_velocity_.clear();
        for (unsigned int i = 0; i < cur_un_pts_.size(); i++)
        {
            if (ids_[i] != -1)
            {
                std::map<int, cv::Point2f>::iterator it;
                it = prev_un_pts_map_.find(ids_[i]);
                if (it != prev_un_pts_map_.end())
                {
                    double v_x = (cur_un_pts_[i].x - it->second.x) / dt;
                    double v_y = (cur_un_pts_[i].y - it->second.y) / dt;
                    pts_velocity_.push_back(cv::Point2f(v_x, v_y));
                }
                else
                    pts_velocity_.push_back(cv::Point2f(0, 0));
            }
            else
            {
                pts_velocity_.push_back(cv::Point2f(0, 0));
            }
        }
    }
    else
    {
        for (unsigned int i = 0; i < cur_pts_.size(); i++)
        {
            pts_velocity_.push_back(cv::Point2f(0, 0));
        }
    }
    prev_un_pts_map_ = cur_un_pts_map_;
}

void mTracker::undistortedPoints_e(){
    cur_un_pts_e.clear();
    cur_un_pts_map_e.clear();
    for (unsigned int i = 0; i < cur_pts_e.size(); i++)
    {
        Eigen::Vector2d a0(cur_pts_e[i].x, cur_pts_e[i].y);
        Eigen::Vector3d b0;
        liftProjective(a0, b0);
        cur_un_pts_e.push_back(cv::Point2f(b0.x() / b0.z(), b0.y() / b0.z()));
        cur_un_pts_map_e.insert(std::make_pair(ids_e[i], cv::Point2f(b0.x() / b0.z(), b0.y() / b0.z())));
    }
    
    // caculate points velocity (left image)
    if (!prev_un_pts_map_e.empty())
    {
        double dt = cur_time_e - prev_time_fix;
        pts_velocity_e.clear();
        for (unsigned int i = 0; i < cur_un_pts_e.size(); i++)
        {
            if (ids_e[i] != -1)
            {
                
                std::map<int, cv::Point2f>::iterator it;
                it = prev_un_pts_map_e.find(ids_e[i]);
                if (it != prev_un_pts_map_e.end())
                {
                    double v_x = (cur_un_pts_e[i].x - it->second.x) / dt;
                    double v_y = (cur_un_pts_e[i].y - it->second.y) / dt;
                    pts_velocity_e.push_back(cv::Point2f(v_x, v_y));
                }
                else
                    pts_velocity_e.push_back(cv::Point2f(0, 0));
            }
            else
            {
                pts_velocity_e.push_back(cv::Point2f(0, 0));
            }
        }
        velocity_use = true;     
    }
    else
    {
        for (unsigned int i = 0; i < cur_pts_e.size(); i++)
        {
            pts_velocity_e.push_back(cv::Point2f(0, 0));
        }
    }
    prev_un_pts_map_e = cur_un_pts_map_e;
}


void mTracker::undistortedLines_e(){
    un_lines.clear();

    for (unsigned int i = 0; i < final_lines.size(); i++) {
        
        Eigen::Vector2d a0(final_lines[i].second[0], final_lines[i].second[1]);
        Eigen::Vector3d b0;
        liftProjective(a0, b0);
        V3d un_pt_s = getK()*b0;
        

        Eigen::Vector2d a1(final_lines[i].second[2], final_lines[i].second[3]);
        Eigen::Vector3d b1;
        liftProjective(a1, b1);
        V3d un_pt_e = getK()*b1;

        V4d un_line(un_pt_s[0], un_pt_s[1], un_pt_e[0], un_pt_e[1]);

        un_lines.emplace_back(final_lines[i].first, un_line);
    }
}

void mTracker::calTempWindow() {
    // set temporal window as N percentile
    std::vector<double> vlifetime(pts_velocity_e.size());
    for (size_t i = 0; i < vlifetime.size(); ++i) {
        double norm_flow = ((cam_intrinsics_[0]+cam_intrinsics_[1])/2)*cv::norm(pts_velocity_e[i]);
        vlifetime[i] = LIFE_PIXEL_ / norm_flow;
        //ROS_WARN("norm_flow: %f", norm_flow);
    }
    std::vector<double> valid_vlifetime;
    for (double lifetime : vlifetime) {
        // must not be inf
        if (lifetime < 1) {
            valid_vlifetime.push_back(lifetime);
        }
    }
    double temp_window;
    if (!valid_vlifetime.empty()) {
        std::sort(valid_vlifetime.begin(), valid_vlifetime.end());
        int index = static_cast<int>(valid_vlifetime.size() * 0.65); // 65th percentile
        temp_window = valid_vlifetime[index];
    } 
    else {
        temp_window = (t_window > 1e-4) ? t_window : 1.0;
    }
    t_window = temp_window;

}




