/**
 * @file NavEarth.cpp
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


#include "NavEarth.h"

V3d NAV::GetWnie(V3d &Vn, V3d &Radpos)
{
    V2d RmRp;
    RmRp = NAV::GetRmRp(Radpos(0));
    double rm = RmRp(0);
    double rp = RmRp(1);
    double dlat, dlon, dh;
    dlat = Vn(0) / (rm + Radpos(2));
    dlon = Vn(1) / ((rp + Radpos(2)) * cos(Radpos(0)));
    dh = -Vn(2);

    double clat, slat;
    clat = cos(Radpos(0));
    slat = sin(Radpos(0));

    V3d Wnie;
    Wnie << ER * clat, 0, -ER * slat;

    return Wnie;
}

V3d NAV::GetWnen(V3d &Vn, V3d &Radpos)
{
    V2d RmRp;
    RmRp = NAV::GetRmRp(Radpos(0));
    double rm = RmRp(0);
    double rp = RmRp(1);
    double dlat, dlon, dh;
    dlat = Vn(0) / (rm + Radpos(2));
    dlon = Vn(1) / ((rp + Radpos(2)) * cos(Radpos(0)));
    dh = -Vn(2);

    double clat, slat;
    clat = cos(Radpos(0));
    slat = sin(Radpos(0));

    V3d Wnen;
    Wnen << dlon * clat, -dlat, -dlon * slat;

    return Wnen;
}

V2d NAV::GetRmRp(const double &lat)
{
    double e2 = ECC * ECC;
    double den = 1 - e2 * sin(lat) * sin(lat);
    V2d RmRp;
    RmRp << R0 * (1 - e2) / (pow(den, 1.5)), R0 / (sqrt(den));
    return RmRp;
}

V2d NAV::GetDelRmRp(const double &lat)
{
    double e2 = ECC * ECC;
    double den = 1 - e2 * sin(lat) * sin(lat);
    V2d RmmRtt;
    RmmRtt << (3 * R0 * (1 - e2) * e2 * sin(lat) * cos(lat)) / pow(den, 5.0 / 2.0),
        (R0 * e2 * sin(lat) * cos(lat)) / pow(den, 3.0 / 2.0);
    return RmmRtt;
}

double NAV::gravity_wgs84(const double &lat, const double &h)
{
    double GM = 3.986004418e14;  // earth's gravitational constant
    double ge = 9.7803253359;    // theoretical gravity at the equator
    double gp = 9.8321849378;    // theoretical gravity at the pole
    double k = 0.00193185265241; // theoretical gravity formula constant
    double m = 0.00344978650684; // m = ER^2 * a^2 * b / GM
    double f = 1 / 298.257223563;

    double e2, a2, slat, slat2, slat4, h2, g1, g2, g3, g4, g5, g6, g7, g8, gout;
    e2 = ECC * ECC;
    a2 = R0 * R0;
    slat = sin(lat);
    slat2 = slat * slat;
    slat4 = slat2 * slat2;
    h2 = h * h;
    g1 = k + e2 / 2.0;
    g2 = -2.0 * (1.0 + f + m) / R0;
    g3 = g1 * g2 + 4.0 * f / R0;
    g4 = 3.0 / a2;
    g5 = g1 * g4;
    g6 = e2 * (3.0 * e2 / 4.0 + k) / 2.0;
    g7 = 4.0 * g1 * f / R0 + g2 * g6;
    g8 = g4 * g6;

    gout = ge * (1.0 + g1 * slat2 + g2 * h + g3 * slat2 * h + g4 * h2 + g5 * slat2 * h2);
    return gout;
}


V4d NAV::UpdateQuat(V3d &wbib, V4d &FormerQuat, const double &dt)
{
    V3d dtheta = dt * wbib;
    double Wx, Wy, Wz;
    Wx = dtheta(0);
    Wy = dtheta(1);
    Wz = dtheta(2);
    double theta = sqrt(Wx * Wx + Wy * Wy + Wz * Wz);

    double ac = cos(theta / 2);
    double as;
    if (theta < 1e-15)
        as = 0;
    else
        as = sin(theta / 2) / theta;

    V4d quatrk;
    M4d Matrk;
    V4d quatnew;
    quatrk << ac, as * Wx, as * Wy, as * Wz;
    Matrk << quatrk(0), -quatrk(1), -quatrk(2), -quatrk(3),
     quatrk(1), quatrk(0), quatrk(3), -quatrk(2),
        quatrk(2), -quatrk(3), quatrk(0), quatrk(1), 
        quatrk(3), quatrk(2), -quatrk(1), quatrk(0);
    quatnew = Matrk * FormerQuat;

    return quatnew;
}

void NAV::NormalizeQuat(V4d &quat)
{
    double mag = quat.norm();
    double q_normal = 1.5 - 0.5 * mag * mag;
    quat = q_normal * quat;
}

V3d NAV::llh2xyz(V3d &llh)
{
    V2d RnRe = GetRmRp(llh(0));
    double x = (RnRe(1)+llh(2)) * cos(llh(0))  * cos(llh(1));
    double y = (RnRe(1)+llh(2)) * cos(llh(0))  * sin(llh(1));
    double z = (RnRe(1)*(1-ECC*ECC)+llh(2)) * sin(llh(0));

    V3d xyz;
    xyz << x,y,z;
    return xyz;
}

V3d NAV::xyz2llh(V3d &xyz)
{
    double x = xyz(0);
    double y = xyz(1);
    double z = xyz(2);  
    double a = 6378137.0000;
    double b = 6356752.3142;
    double e = pow(1-pow((b/a), 2), 0.5);
    double ep = e*(a/b);
    double r = pow(pow(x,2)+pow(y,2),0.5);
    double E2 = pow(a,2) - pow(b,2);
    double F = 54*pow(b,2)*pow(z,2);
    double G = pow(r,2) 
            + (1 - pow(e,2))*pow(z,2) - pow(e,2)*E2;
    double c = (pow(e,2)*pow(e,2)*F*pow(r,2))/pow(G,3);
    double s = pow(1+c+pow(pow(c,2)+2*c,0.5),1/3);
    double P = F/(3*pow((s + 1/s + 1),2)*pow(G,2));
    double Q = pow(1+2*pow(e,2)*pow(e,2)*P,0.5);
    double ro = -(P*pow(e,2)*r)/(1+Q) + pow((a*a/2)*(1+1/Q),0.5)
                -(P*(1-pow(e,2))*pow(z,2))/(Q*(1+Q))-P*pow(r,2)/2;
    double tmp = pow((r - pow(e,2)*ro),2);
    double U = pow(tmp + pow(z,2), 0.5);
    double V = pow(tmp + (1- pow(e,2))*pow(z,2),0.5);
    double zo = pow(b,2)*z/(a*V);
    double height = U * (1-pow(b,2)/(a*V));
    double lat = atan((z+ep*ep*zo)/r);
    double temp = atan(y/x);
    double longt;

    if(x>=0){
        longt = temp;
    }
    else if(x < 0 && y >= 0){
        longt = M_PI + temp;
    }
    else longt = temp - M_PI;
    V3d llh;
    llh(0) = lat;
    llh(1) = longt;
    llh(2) = height;
    return llh;
}

V3d NAV::xyz2ned(V3d &xyz, V3d &orgxyz)
{
    V3d tmpxyz = xyz;
    V3d tmporg = orgxyz;   
    V3d difxyz = tmpxyz - tmporg;
    V3d orgllh = NAV::xyz2llh(orgxyz);
    double phi = orgllh(0);
    double lam = orgllh(1);
    double sinphi = sin(phi);
    double cosphi = cos(phi);
    double sinlam = sin(lam);
    double coslam = cos(lam);
    M3d R;
    R << -sinlam,           coslam,      0,
        -sinphi*coslam,  -sinphi*sinlam, cosphi,
        cosphi*coslam,   cosphi*sinlam, sinphi;
    V3d enu = R * difxyz;
    V3d ned;
    ned(0) = enu(1);
    ned(1) = enu(0);
    ned(2) = -enu(2);
    return ned;
}
