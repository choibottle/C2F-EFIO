/**
 * @file Attitude.h
 * @author Jaehyung Jung (lastflowers@snu.ac.kr)
 * @brief Attitude transform functions
 * @date 2021-03-30
 *
 * @copyright Copyright (c) 2020 Jaehyung Jung
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

#include <Eigen/Eigen>
#include "define.h"

#ifndef ATTITUDE_H
#define ATTITUDE_H

using namespace Eigen;

namespace nesl {

    M3d euler2dcm(const V3d &euler);
    V4d euler2quat(const V3d &euler);

    M3d quat2dcm(const V4d &quat);
    V3d quat2euler(const V4d &quat);
    V3d quat2rvec(const V4d &quat);
    M4d quatLeftComp(const V4d &quat);
    M4d quatRightComp(const V4d &quat);
    V4d quatInverse(const V4d &quat);
    V4d quatNormalize(const V4d &quat);
    V4d quatMultiply(const V4d &q, const V4d &p);

    V3d dcm2euler(const M3d &dcm);
    V4d dcm2quat(const M3d &dcm);
    V3d dcm2rvec(const M3d &dcm);

    M3d Vec2SkewMat(const V3d &vec);
    V3d SkewMat2Vec(const M3d &mat);
    M3d rvec2dcm(const V3d &vec);
    V4d rvec2quat(const V3d &rvec);

    M3d expmap(const V3d &phi);
    V3d logmap(const M3d &R);

    Eigen::MatrixXd matrixSqrt(const Eigen::MatrixXd& A);
    M3d validate(const M3d &R);
}

#endif // ATTITUDE_H
