/**
 * @file NavEarth.h
 * @author Jaehyung Jung (lastflowers@snu.ac.kr)
 *         Hanyeol Lee (han2110@snu.ac.kr)
 * @brief 
 * @date 2021-06-14
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


// Navigation on Earth related functions

#include "define.h"

#ifndef NAV_EARTH_H
#define NAV_EARTH_H

namespace NAV
{

    V3d GetWnie(V3d &Vn, V3d &Radpos);
    V3d GetWnen(V3d &Vn, V3d &Radpos);
    V2d GetRmRp(const double &lat); // Meridian radius, transverse radius
    V2d GetDelRmRp(const double &lat);
    double gravity_wgs84(const double &lat, const double &h);    
    V4d UpdateQuat(V3d &wbib, V4d &FormerQuat, const double &dt);
    void NormalizeQuat(V4d &quat);
    V3d llh2xyz(V3d &llh); // llh2ecef in WGS-84
    V3d xyz2llh(V3d &xyz);
    V3d xyz2ned(V3d &xyz, V3d &orgxyz);
}

#endif