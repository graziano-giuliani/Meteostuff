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

#ifndef __WEATHER_H__
#define __WEATHER_H__

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <boost/date_time/posix_time/posix_time.hpp>

namespace himet {

  using namespace boost::posix_time;

  typedef enum {
    SKY_clear          = 0,
    SKY_mostly_clear   = 1,
    SKY_partly_cloudy  = 2,
    SKY_cloudy         = 3,
    SKY_covered        = 4,
    SKY_light_rain     = 5,
    SKY_medium_rain    = 6,
    SKY_heavy_rain     = 7,
    SKY_thunderstorm   = 8,
    SKY_snowing        = 9,
    SKY_fog            = 10,
    SKY_veiled         = 11,
    SKY_variable_cloud = 12,
    SKY_variable_rain  = 13,
    SKY_variable_storm = 14,
    SKY_variable_snow  = 15
  } t_enum_sky_condition;

  typedef enum {
    TT_stationary = 0,
    TT_increase   = 1,
    TT_decrease   = 2
  } t_enum_tendency;

  typedef enum {
    WDIR_absent   = 0,
    WDIR_N        = 1,
    WDIR_NNE      = 2,
    WDIR_NE       = 3,
    WDIR_ENE      = 4,
    WDIR_E        = 5,
    WDIR_ESE      = 6,
    WDIR_SE       = 7,
    WDIR_SSE      = 8,
    WDIR_S        = 9,
    WDIR_SSW      = 10,
    WDIR_SW       = 11,
    WDIR_WSW      = 12,
    WDIR_W        = 13,
    WDIR_WNW      = 14,
    WDIR_NW       = 15,
    WDIR_NNW      = 16,
    WDIR_variable = 17
  } t_enum_wind_direction;

  typedef enum {
    WSPD_absent      = 0,
    WSPD_calm        = 1,
    WSPD_light       = 2,
    WSPD_moderate    = 3,
    WSPD_strong      = 4,
    WSPD_very_strong = 5,
    WSPD_storm_wind  = 6
  } t_enum_wind_intensity;

  typedef enum {
    VIS_optimal    = 0,
    VIS_good       = 1,
    VIS_fair       = 2,
    VIS_light_haze = 3,
    VIS_haze       = 4,
    VIS_fog        = 5
  } t_enum_visibility;

  typedef enum {
    KIND_absent       = 0,
    KIND_intermittent = 1,
    KIND_persistent   = 2
  } t_enum_phenomena_kind;

  typedef enum {
    RAIN_absent       = 0,
    RAIN_light        = 1,
    RAIN_moderate     = 2,
    RAIN_heavy        = 3,
    RAIN_shower       = 4,
    RAIN_heavy_shower = 5,
    RAIN_storm        = 6,
    RAIN_heavy_storm  = 7
  } t_enum_rain_kind;

  typedef enum {
    SNOW_absent     = 0,
    SNOW_light      = 1,
    SNOW_moderate   = 2,
    SNOW_heavy      = 3,
    SNOW_very_heavy = 4,
    SNOW_extreme    = 5
  } t_enum_snow_kind;

  typedef enum {
    SLOT_night     = 0,
    SLOT_morning   = 1,
    SLOT_afternoon = 2,
    SLOT_evening   = 3
  } t_enum_dayslot;

  class weather {
    public:
      weather( );
      std::string &skycond(int code);
      std::string &tendency(int code);
      std::string &wdir(int code);
      std::string &wname(int code);
      std::string &wspd(int code);
      std::string &visibility(int code);
      std::string &phenkind(int code);
      std::string &rain(int code);
      std::string &snow(int code);
      std::string &dayslot(int code);
      t_enum_dayslot UTCtime2UTCslot(ptime utc);
      t_enum_dayslot UTCtime2LOCslot(ptime utc, std::string zonedef);
      void uv2sd(const std::vector <float> &u, const std::vector <float> &v,
                 std::vector <float> &s, std::vector <float> &d);
      void sd2desc(const std::vector <float> &s, const std::vector <float> &d,
                   std::vector <t_enum_wind_direction> &dir,
                   std::vector <t_enum_wind_intensity> &spd);
      void perceived(const std::vector <float> &tc,
                     const std::vector <float> &rh,
                     const std::vector <float> &spd,
                     std::vector <float> &pt);
      void calcvis(const std::vector <float> &tc,
                   const std::vector <float> &prs,
                   const std::vector <float> &tdc,
                   const std::vector <float> &q,
                   const std::vector <float> &clw,
                   const std::vector <float> &rnw,
                   const std::vector <float> &snow,
                   const std::vector <float> &ice,
                   std::vector <float> &vis);
       void visdesc(const std::vector <float> &vis,
                    std::vector <t_enum_visibility> &evis);
    private:
      std::map<int, std::string> sky_condition_map;
      std::map<int, std::string> tendency_map;
      std::map<int, std::string> wind_direction_map;
      std::map<int, std::string> wind_direction_name_map;
      std::map<int, std::string> wind_intensity_map;
      std::map<int, std::string> visibility_map;
      std::map<int, std::string> phenomena_kind_map;
      std::map<int, std::string> rain_kind_map;
      std::map<int, std::string> snow_kind_map;
      std::map<int, std::string> dayslot_map;
  };

}

#endif
