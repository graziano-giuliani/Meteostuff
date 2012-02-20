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

#include <weather.h>
#include "boost/date_time/local_time/local_time.hpp"

using namespace himet;
using namespace boost::local_time;

static const int max_name_len = 64;

static const int num_sky_cond = 16;
static const int num_tendency = 3;
static const int num_wdirs    = 18;
static const int num_wdirname = 18;
static const int num_wspds    = 7;
static const int num_viscond  = 6;
static const int num_kinds    = 3;
static const int num_rainkind = 8;
static const int num_snowkind = 6;
static const int num_dayslot  = 4;

static char sky_condition_names[num_sky_cond][max_name_len] = {
  "cielo sereno",
  "cielo prevalentemente sereno",
  "cielo parzialmente nuvoloso",
  "cielo molto nuvoloso",
  "cielo coperto",
  "pioggia",
  "pioggia",
  "pioggia",
  "temporale",
  "neve",
  "nebbia",
  "cielo velato",
  "nuvolosit&agrave; variabile",
  "nuvolosit&agrave; variabile con pioggia",
  "nuvolosit&agrave; variabile con temporali",
  "nuvolosit&agrave; variabile con neve"
};

static char tendency_names[num_tendency][max_name_len] = {
  "stazionaria", "in aumento", "in diminuzione"
};

static char wind_direction_names[num_wdirs][max_name_len] = {
  "assente",
  "N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
  "S", "SSO", "SO", "OSO", "O", "ONO", "NO", "NNO",
  "da direzioni variabili"
};

static char wind_direction_name_names[num_wdirname][max_name_len] = {
  "assente",
  "Tramontana", "Tramontana", "Grecale", "Grecale",
  "Levante", "Levante", "Scirocco", "Scirocco",
  "Mezzogiorno", "Mezzogiorno", "Libeccio", "Libeccio",
  "Ponente", "Ponente", "Maestrale", "Maestrale",
  "da direzioni variabili"
};

static char wind_intensity_names[num_wspds][max_name_len] = {
  "vento assente", "vento calmo", "vento debole",
  "vento moderato", "vento forte", "vento molto forte",
  "tempesta di vento"
};

static char visibility_names[num_viscond][max_name_len] = {
  "visibilit&agrave; ottima", "visibilit&agrave; buona",
  "visibilit&agrave; discreta", "leggera foschia",
  "foschia", "nebbia"
};

static char phenomena_kind_names[num_kinds][max_name_len] = {
  "assente",
  "a tratti",
  "persistente"
};

static char rain_kind_names[num_rainkind][max_name_len] = {
  "assente",
  "pioggia debole",
  "pioggia moderata",
  "pioggia forte",
  "rovesci di pioggia",
  "nubifragio",
  "a carattere temporalesco",
  "a violento carattere temporalesco"
};

static char snow_kind_names[num_snowkind][max_name_len] = {
  "assente", "deboli nevicate", "nevicate moderate",
  "forti nevicate", "nevicate molto forti", "nevicate estreme"
};

static char dayslot_names[num_dayslot][max_name_len] = {
  "notte", "mattina", "pomeriggio", "sera"
};

static const double lowmark_precip = 0;
static const double lowmark_convective_precip = 10;
static const double R2D = 180.0/M_PI;
static const double D2R = M_PI/180.0;
static const double windlow = 0.001;
static inline float kelvin_to_celsius(float k);
static inline float kelvin_to_fahrenheit(float k);
static inline float celsius_to_fahrenheit(float c);
static inline float fahrenheit_to_celsius(float f);
static inline float Pa_to_mb(float pa);
static inline float kmh2knt(float spd);

inline float kelvin_to_celsius(float k) { return k-273.15; }
inline float kelvin_to_fahrenheit(float k) { return 1.8*k-459.67; }
inline float celsius_to_fahrenheit(float c) { return 1.8*c+32.0; }
inline float fahrenheit_to_celsius(float f) { return 5.0/9.0*(f-32.0); }
inline float Pa_to_mb(float pa) { return pa/100.0; }
inline float kmh2knt(float spd) { return spd/1.852; }

static std::string &search_and_report(std::map<int, std::string> &where,
                                      int what);

weather::weather( )
{
  for (int i = 0; i < num_sky_cond; i ++)
    sky_condition_map.insert(
          std::pair<int, std::string>(i, sky_condition_names[i]));
  for (int i = 0; i < num_tendency; i ++)
    tendency_map.insert(
          std::pair<int, std::string>(i, tendency_names[i]));
  for (int i = 0; i < num_wdirs; i ++)
    wind_direction_map.insert(
          std::pair<int, std::string>(i, wind_direction_names[i]));
  for (int i = 0; i < num_wdirname; i ++)
    wind_direction_name_map.insert(
          std::pair<int, std::string>(i, wind_direction_name_names[i]));
  for (int i = 0; i < num_wspds; i ++)
    wind_intensity_map.insert(
          std::pair<int, std::string>(i, wind_intensity_names[i]));
  for (int i = 0; i < num_viscond; i ++)
    visibility_map.insert(
          std::pair<int, std::string>(i, visibility_names[i]));
  for (int i = 0; i < num_kinds; i ++)
    phenomena_kind_map.insert(
          std::pair<int, std::string>(i, phenomena_kind_names[i]));
  for (int i = 0; i < num_rainkind; i ++)
    rain_kind_map.insert(
          std::pair<int, std::string>(i, rain_kind_names[i]));
  for (int i = 0; i < num_snowkind; i ++)
    snow_kind_map.insert(
          std::pair<int, std::string>(i, snow_kind_names[i]));
  for (int i = 0; i < num_dayslot; i ++)
    dayslot_map.insert(
          std::pair<int, std::string>(i, dayslot_names[i]));
}

std::string &search_and_report(std::map<int, std::string> &where, int what)
{
  std::map<int, std::string>::iterator iter = where.find(what);
  if (iter != where.end())
    return iter->second;
  else
    throw "Item not found";
}

std::string &weather::skycond(int code)
{
  return search_and_report(sky_condition_map, code);
}

std::string &weather::tendency(int code)
{
  return search_and_report(tendency_map, code);
}

std::string &weather::wdir(int code)
{
  return search_and_report(wind_direction_map, code);
}

std::string &weather::wspd(int code)
{
  return search_and_report(wind_intensity_map, code);
}

std::string &weather::wname(int code)
{
  return search_and_report(wind_direction_name_map, code);
}

std::string &weather::visibility(int code)
{
  return search_and_report(visibility_map, code);
}

std::string &weather::phenkind(int code)
{
  return search_and_report(phenomena_kind_map, code);
}

std::string &weather::rain(int code)
{
  return search_and_report(rain_kind_map, code);
}

std::string &weather::snow(int code)
{
  return search_and_report(snow_kind_map, code);
}

std::string &weather::dayslot(int code)
{
  return search_and_report(dayslot_map, code);
}

t_enum_dayslot weather::UTCtime2UTCslot(ptime utc)
{
  time_duration td = hours(1);
  if (utc.time_of_day() < td*6)  return SLOT_night;
  if (utc.time_of_day() < td*12) return SLOT_morning;
  if (utc.time_of_day() < td*18) return SLOT_afternoon;
  return SLOT_evening;
}

t_enum_dayslot weather::UTCtime2LOCslot(ptime utc, std::string zn)
{
  time_zone_ptr zone(new posix_time_zone(zn));
  time_duration shift = -zone->base_utc_offset() + zone->dst_offset();
  ptime ltime = utc + shift;
  return UTCtime2UTCslot(ltime);
}

void weather::uv2sd(const std::vector <float> &u, const std::vector <float> &v,
                    std::vector <float> &s, std::vector <float> &d)
{
  if ( u.size( ) != v.size( ) ||
       s.size( ) != d.size( ) ||
       s.size( ) != u.size( ) )
    throw "Different size of vectors in uv2sd";

  double ux, vx;
  for (int i = 0; i < (int) u.size( ); i ++)
  {
    ux = u[i];
    vx = v[i];
    s[i] = 3.6 * sqrt(ux*ux+vx*vx);
    d[i] = atan2(ux, vx)*R2D+180.0;
  }
  return;
}

void weather::sd2desc(const std::vector <float> &spd,
                      const std::vector <float> &dir,
                      std::vector <t_enum_wind_direction> &edir,
                      std::vector <t_enum_wind_intensity> &espd)
{
  if ( spd.size( ) !=  dir.size( ) ||
       spd.size( ) != edir.size( ) ||
       spd.size( ) != espd.size( ) )
    throw "Different size of vectors in sd2desc";

  for (int i = 0; i < (int) spd.size( ); i ++)
  {
    if (spd[i] < windlow)
    {
      edir[i] = WDIR_absent;
      espd[i] = WSPD_absent;
      continue;
    }

    if      (dir[i] >=   0 && dir[i] <    5) edir[i] = WDIR_N;
    else if (dir[i] >=   5 && dir[i] <=  40) edir[i] = WDIR_NNE;
    else if (dir[i] >   40 && dir[i] <   50) edir[i] = WDIR_NE;
    else if (dir[i] >=  50 && dir[i] <=  85) edir[i] = WDIR_ENE;
    else if (dir[i] >   85 && dir[i] <   95) edir[i] = WDIR_E;
    else if (dir[i] >=  95 && dir[i] <= 130) edir[i] = WDIR_ESE;
    else if (dir[i] >  130 && dir[i] <  140) edir[i] = WDIR_SE;
    else if (dir[i] >= 140 && dir[i] <= 175) edir[i] = WDIR_SSE;
    else if (dir[i] >  175 && dir[i] <  185) edir[i] = WDIR_S;
    else if (dir[i] >= 185 && dir[i] <= 220) edir[i] = WDIR_SSW;
    else if (dir[i] >  220 && dir[i] <  230) edir[i] = WDIR_SW;
    else if (dir[i] >= 230 && dir[i] <= 265) edir[i] = WDIR_WSW;
    else if (dir[i] >  265 && dir[i] <  275) edir[i] = WDIR_W;
    else if (dir[i] >= 275 && dir[i] <= 310) edir[i] = WDIR_WNW;
    else if (dir[i] >  310 && dir[i] <  320) edir[i] = WDIR_NW;
    else if (dir[i] >= 320 && dir[i] <= 355) edir[i] = WDIR_NNW;
    else edir[i] = WDIR_N;

    if      (spd[i] >  0 && spd[i] <=  5) espd[i] = WSPD_calm;
    else if (spd[i] >  5 && spd[i] <= 19) espd[i] = WSPD_light;
    else if (spd[i] > 19 && spd[i] <= 38) espd[i] = WSPD_moderate;
    else if (spd[i] > 38 && spd[i] <= 61) espd[i] = WSPD_strong;
    else if (spd[i] > 61 && spd[i] <= 88) espd[i] = WSPD_very_strong;
    else espd[i] = WSPD_storm_wind;
  }

  return;
}

void weather::perceived(const std::vector <float> &tc,
                        const std::vector <float> &rh,
                        const std::vector <float> &spd,
                        std::vector <float> &pt)
{
  if ( tc.size( ) != rh.size( )  ||
       tc.size( ) != spd.size( ) ||
       tc.size( ) != pt.size( ) )
    throw "Different size of vectors in perceived";

  double cels, kmh, hum, rft;
  for (int i = 0; i < (int) tc.size( ); i ++)
  {
    cels = tc[i] < -50.0 ? -50.0 : tc[i];
    kmh  = spd[i] > 100.0 ? 100.0 : spd[i];
    hum  = rh[i] < 0 ? 0.0 : rh[i];

    if (cels <= 10.0)
    {
      if (kmh >= 5.0)
      {
        rft = 13.12 + 0.6215*cels - 11.37*pow(kmh, 0.16)
                    + 0.3965*cels*pow(kmh, 0.16);
      }  
      else
        rft = cels;
    }
    else if (cels >= 27.0)
    {
      if (hum >= 40.0 && kmh < 20.0)
      {
        float far = celsius_to_fahrenheit(cels);
        rft = -42.379 + 2.04901523*far + 10.14333127*hum
                      - 0.22475541*far*hum - 0.00683783*(far*far)
                      - 0.05481717*(hum*hum) + 0.00122874*(far*far)*hum
                      + 0.00085282*far*(hum*hum)
                      - 0.00000199*(far*far)*(hum*hum);
        rft = fahrenheit_to_celsius(rft);
      }
      else
        rft = cels;
    }
    else
      rft = cels;

    pt[i] = nearbyint(rft);
  }
}

void weather::calcvis(const std::vector <float> &tc,
                      const std::vector <float> &prs,
                      const std::vector <float> &tdc,
                      const std::vector <float> &q,
                      const std::vector <float> &clw,
                      const std::vector <float> &rnw,
                      const std::vector <float> &snow,
                      const std::vector <float> &ice,
                      std::vector <float> &vis)
{
  if ( tc.size( ) != prs.size( )  ||
       tc.size( ) != tdc.size( )  ||
       tc.size( ) != q.size( )    ||
       tc.size( ) != clw.size( )  ||
       tc.size( ) != rnw.size( )  ||
       tc.size( ) != snow.size( ) ||
       tc.size( ) != ice.size( )  ||
       tc.size( ) != vis.size( )   )
    throw "Different size of vectors in calcvis";

  double beta, rhoair, tv, qrpc, qcld;
  double conclc, conclp, concfc, concfp, vovermd;
  for (int i = 0; i < (int) tc.size( ); i ++)
  {
    tv = tc[i] * (1.0 + 0.608 * q[i]);
    rhoair = prs[i] / (287.04 * tv);
    qrpc = rnw[i] + snow[i];
    qcld = clw[i] + ice[i];
    vovermd = (1.0 + q[i])       / rhoair +
              (clw[i] + rnw[i])  / 1000.0 +
              (ice[i] + snow[i]) / 917.0; 
    conclc = 1000.0 * clw[i] / vovermd;
    conclp = 1000.0 * rnw[i] / vovermd;
    concfc = 1000.0 * ice[i] / vovermd;
    concfp = 1000.0 * snow[i] / vovermd;

    beta = 327.80 * concfc +
            10.36 * pow(concfp, 0.7776) +
           144.70 * pow(conclc, 0.8800) +
             2.24 * pow(conclp, 0.7500) + 1.e-10;
    vis[i] = 1000.0 * (-log(0.02) / beta);
    if (vis[i] > 25000.0) vis[i] = 25000.0;
  }
  return;
}

void weather::visdesc(const std::vector <float> &vis,
                      std::vector <t_enum_visibility> &evis)
{
  if ( vis.size( ) != evis.size( ) )
    throw "Different size of vectors in visdesc";

  for (int i = 0; i < (int) vis.size( ); i ++)
  {
    if      (vis[i] > 10000.0)                      evis[i] = VIS_optimal;
    else if (vis[i] >= 4000.0 && vis[i] <= 10000.0) evis[i] = VIS_good;
    else if (vis[i] >= 2000.0 && vis[i] <= 4000.0)  evis[i] = VIS_fair;
    else if (vis[i] >= 1000.0 && vis[i] <= 2000.0)  evis[i] = VIS_light_haze;
    else if (vis[i] >= 100.0  && vis[i] <= 1000.0)  evis[i] = VIS_haze;
    else evis[i] = VIS_fog;
  }
  return;
}

#ifdef TESTME

#include <boost/date_time/gregorian/gregorian.hpp>
using namespace boost::gregorian;

int main(int argc, char *argv[])
{
  weather w;
  std::cout << w.skycond(SKY_covered) << std::endl;
  ptime T(date(2009, Oct, 18), hours(17)+minutes(11));
  std::cout << w.dayslot(w.UTCtime2UTCslot(T)) << std::endl;
  std::string zn = "CET-1CEST,M3.5.0,M10.5.0/3";
  std::cout << w.dayslot(w.UTCtime2LOCslot(T, zn)) << std::endl;
  std::vector <float> u(1), v(1), t(1), rh(1), prs(1), tdc(1);
  std::vector <float> q(1), clw(1), rnw(1), snow(1), ice(1);
  std::vector <float> s(1), d(1), pt(1), vis(1);
  std::vector <t_enum_wind_direction> ed(1);
  std::vector <t_enum_wind_intensity> ei(1);
  std::vector <t_enum_visibility> ev(1);
  u[0]    =   40.0;
  v[0]    =  -10.0;
  t[0]    =   30.0;
  rh[0]   =   60.0;
  prs[0]  = 1000.0;
  tdc[0]  =   25.0;
  q[0]    =   0.01;
  clw[0]  =   0.0001;
  rnw[0]  =   0.001;
  snow[0] =   0.003;
  ice[0]  =   0.00001;
  w.uv2sd(u, v, s, d);
  w.sd2desc(s, d, ed, ei);
  std::cout << w.wdir(ed[0]) << " (" << w.wname(ed[0]) << ")" << std::endl;
  std::cout << w.wspd(ei[0]) << std::endl;
  w.perceived(t,rh,s, pt);
  std::cout << pt[0] << std::endl;
  t[0] =  8.0;
  w.perceived(t,rh,s, pt);
  std::cout << pt[0] << std::endl;
  w.calcvis(t, prs, tdc, q, clw, rnw, snow, ice, vis);
  std::cout << "vis : " << vis[0] << std::endl;
  w.visdesc(vis, ev);
  std::cout << w.visibility(ev[0]) << std::endl;
  return 0;
}

#endif
