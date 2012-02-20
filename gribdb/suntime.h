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

#ifndef __SUNRISE__H__
#define __SUNRISE__H__

#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <iostream>
#include <cmath>

namespace himet {

  using namespace boost::gregorian;
  using namespace boost::posix_time;

  typedef enum {
    SUNHGT_CENTER            = 0,
    SUNHGT_UPPER             = 1,
    SUNHGT_CENTER_REFRACTION = 2,
    SUNHGT_UPPER_REFRACTION  = 3,
    SUNHGT_CIVIL             = 4,
    SUNHGT_ASTRO_AMATEUR     = 5,
    SUNHGT_ASTRONOMICAL      = 6
  } t_enum_refhgt;

  static double refhgt_degval[8] = {
        0, -0.25, -0.583, -0.833, -6, -12, -15, -18 };
  static char refhgt_str[8][72]  = {
    "Center of Sun's disk touches a mathematical horizon",
    "Sun's upper limb touches a mathematical horizon",
    "Center of Sun's disk touches the horizon with atmospheric refraction",
    "Sun's upper limb touches the horizon with atmospheric refraction",
    "Civil twilight (no reading outside without artificial illumination)",
    "Amateur astronomical twilight (sky dark for astronomical observation)",
    "Astronomical twilight (sky completely dark)"
  };

  class refhgt {
    public:
      refhgt( ) { indx = 3; }
      refhgt(t_enum_refhgt type) { indx = (int) type; }
      void change(t_enum_refhgt type) { indx = (int) type; }
      double hgtdeg( ) { return refhgt_degval[indx]; }
      friend std::ostream& operator<< (std::ostream &out, refhgt &type)
      {
         out << refhgt_str[type.indx];
         return out;
      }
    private:
      int indx;
  };

  class suninfo {
    public:
      ptime sunrise;
      ptime sunset;
      friend std::ostream& operator<< (std::ostream &out, suninfo &i)
      {
        out << "Sun rises at " << i.sunrise << ", sets at " << i.sunset;
        return out;
      }
  };

  typedef struct {
    double right_ascension;
    double declination;
    double ecliptic_longitude;
    double distance;
  } sunpos;

  class suntime {
    public:
      suntime( );
      inline void setref(t_enum_refhgt t) { h.change(t); }
      void calculate(date d, double lon, double lat, suninfo &output);
      void localized(suninfo &utc, suninfo &local, std::string zn);
    private:
      void workhorse(double d, double lon, double lat,
                     double h, double &rise, double &set);
      inline double greenwich_mean_sideral_time(double d);
      void calcpos(double d, sunpos &a);
      int upper_limb;
      refhgt h;
      date fxd;
  };

}

#endif
