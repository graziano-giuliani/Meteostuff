/***************************************************************************
 *   Copyright (C) 2008-2009 by Graziano Giuliani                          *
 *   graziano.giuliani at aquila.infn.it                                   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details. (see COPYING)            *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 *                                                                         *
 *   LIC: GPL                                                              *
 *                                                                         *
 ***************************************************************************/

#include <suntime.h>
#include <iostream>
#include <cmath>
#include <boost/date_time/local_time/local_time.hpp>

using namespace himet;
using namespace boost::local_time;

static const double RADEG = 180.0 / M_PI;
static const double DEGRAD = M_PI / 180.0;
static const double INV360 = 1.0 / 360.0;
static const double dlen = 15.04107;
static const double accuracy = 0.0001; // half second

static inline double revolution(double angle);
static inline double rev180(double angle);
static inline double sind(double angle);
static inline double cosd(double angle);
static inline double tand(double angle);
static inline double atand(double angle);
static inline double asind(double angle);
static inline double acosd(double angle);
static inline double atan2d(double a1, double a2);

inline double revolution(double angle)
{
  return (angle - 360.0 * floor( angle * INV360 ));
}
inline double rev180(double angle)
{
  return (angle - 360.0 * floor( angle * INV360 + 0.5 ));
}
inline double sind(double angle) { return sin(angle*DEGRAD); }
inline double cosd(double angle) { return cos(angle*DEGRAD); }
inline double tand(double angle) { return tan(angle*DEGRAD); }
inline double atand(double angle) { return RADEG * atan(angle); }
inline double asind(double angle) { return RADEG * asin(angle); }
inline double acosd(double angle) { return RADEG * acos(angle); }
inline double atan2d(double a1, double a2) { return RADEG*atan2(a1,a2); }

suntime::suntime( )
{
  upper_limb = 1;
  fxd = date(2000, Jan, 1);
}

void suntime::calculate(date d, double lon, double lat, suninfo &i)
{
  double sd = (d - fxd).days( ) + 0.5 - lon / 360.0;
  double r1, r2, r3, s1, s2, s3, xx;
  double alt = h.hgtdeg( );
  workhorse(sd, lon, lat, alt, r1, s1);
  r2 = 9;
  r3 = 0;
  while ( fabs(r2 - r3) > accuracy )
  {
    double d1 = sd + r1 / 24.0;
    workhorse(d1, lon, lat, alt, r2, xx);
    r1 = r3;
    double d2 = sd + r2 / 24.0;
    workhorse(d2, lon, lat, alt, r3, xx);
  }
  s2 = 9;
  s3 = 0;
  while ( fabs(s2 - s3) > accuracy )
  {
    double d1 = sd + s1 / 24.0;
    workhorse(d1, lon, lat, alt, xx, s2);
    s1 = s3;
    double d2 = sd + s2 / 24.0;
    workhorse(d2, lon, lat, alt, xx, s3);
  }

  i.sunrise = ptime(d)+seconds(r3*3600.0);
  i.sunset = ptime(d)+seconds(s3*3600.0);
  return;
}

void suntime::workhorse(double d, double lon, double lat,
                        double h, double &rise, double &set)
{
  sunpos p;
  calcpos(d, p);
  double sidtime = revolution(greenwich_mean_sideral_time(d) + 180.0 + lon);
  double tsouth  = 12.0 - rev180(sidtime - p.right_ascension) / dlen;
  double sradius = 0.2666 / p.right_ascension;
  if (upper_limb) h -= sradius;
  double cost = ( sind(h) - sind(lat) * sind(p.declination) ) /
                ( cosd(lat) * cosd(p.declination) );
  double t;
  if ( cost >= 1.0 ) t = 0.0;        // Sun always below h
  else if ( cost <= -1.0 ) t = 12.0; // Sun always above h
  else t = acosd(cost) / dlen;

  rise = tsouth - t;
  set =  tsouth + t;

  return;
}

inline double suntime::greenwich_mean_sideral_time(double d)
{
  return revolution( ( 180.0 + 356.0470 + 282.9404 ) +
                    ( 0.9856002585 + 4.70935E-5 ) * d );
}

void suntime::calcpos(double d, sunpos &a)
{
  double Mean_anomaly_of_sun = revolution(356.0470 + 0.9856002585*d);
  double Mean_longitude_of_perihelion = 282.9404 + 4.70935e-5 * d;
  double Eccentricity_of_Earth_orbit  = 0.016709 - 1.151e-9 * d;
  double Eccentric_anomaly = Mean_anomaly_of_sun +
         Eccentricity_of_Earth_orbit * RADEG * sind(Mean_anomaly_of_sun) *
         ( 1.0 + Eccentricity_of_Earth_orbit * cosd(Mean_anomaly_of_sun) );
  double x = cosd(Eccentric_anomaly) - Eccentricity_of_Earth_orbit;
  double y = sqrt( 1.0 - Eccentricity_of_Earth_orbit *
            Eccentricity_of_Earth_orbit ) * sind(Eccentric_anomaly);
  double True_anomaly = atan2d( y, x );
  a.distance = sqrt( x * x + y * y );
  a.ecliptic_longitude = True_anomaly + Mean_longitude_of_perihelion;
  if ( a.ecliptic_longitude >= 360.0 ) { a.ecliptic_longitude -= 360.0; }
  x = a.distance * cosd(a.ecliptic_longitude);
  y = a.distance * sind(a.ecliptic_longitude);
  double obl_ecl = 23.4393 - 3.563e-7 * d;
  double z = y * sind(obl_ecl);
  y = y * cosd(obl_ecl);
  a.right_ascension = atan2d( y, x );
  a.declination = atan2d( z, sqrt( x * x + y * y ) );
  return;
}

void suntime::localized(suninfo &utc, suninfo &local, std::string zn)
{
  time_zone_ptr zone(new posix_time_zone(zn));
  time_duration shift = -zone->base_utc_offset() + zone->dst_offset();
  local.sunrise = utc.sunrise + shift;
  local.sunset = utc.sunset + shift;
}

#ifdef TESTME

int main(int argc, char *argv[])
{
  suntime e;
  suninfo u, l;
  date d(2009, Oct, 15);

  for (int i = 0; i < 8000; i ++)
  {
    e.calculate(d, 13.37, 42.37, u);
    e.localized(u, l, "CET-1CEST,M3.5.0,M10.5.0/3");
    std::cout << "Localtime " << l << std::endl;
  }
  return 0;
}

#endif
