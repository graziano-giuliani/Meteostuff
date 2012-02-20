/***************************************************************************
 *   Copyright (C) 2007 by Graziano Giuliani                               *
 *   graziano.giuliani at poste.it                                         *
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
// $Id: $

#ifndef METEOSALUTE_H_
#define METEOSALUTE_H_

#include <string>
#include <cmath>

namespace Meteo
{
  struct vectorcomp_field
  {
      float *u;
      float *v;
  };

  struct dirspeed_field
  {
      float *direction;
      float *speed;
  };

  static const float d2r = 3.14159265358979323846 / 180.0;
  static const float r2d = 180.0 / 3.14159265358979323846;
  static const float TWO_M_PI= 3.14159265358979323846 * 2.0;

  class Meteosalute
  {
    public:
      Meteosalute(int vector_size);

      int vector_size;

      dirspeed_field *dirspeed(vectorcomp_field &uv); // Components to module and direction

      float *tempktoC(float *tempk);
      
      float *pressPatomb(float *pa);
      
      float *pressvap(float *temp);
      
      float *defsat(float *dewpoint, float *pvap);
      
      float *temprad(float *tsoil, // Soil temperature in K
                     float *rshort, // Short wave radiation in W/m^2
                     float *rdiffuse, // Diffuse radiation in W/m^2
                     float *sunelev, // Sun elevation
                     float *wvtens); // Water vapour tension

      float *sunelev(float *lat, // Latitude in degrees North
                     float *lon, // Longitude in degrees East
                     float model_time); // Time in seconds (Unix Time)

      float *radteoric(float *lat, // Latitude in degrees North
                       float *topography, // Topographic HGT in m
                       float model_time); // Time in seconds (Unix Time)

      float *rdiffuse(float *radteoric, // Teoric SW radiation
                      float *rshort); // Actual SW radiation down

      float *tempradsoil(float *tsoil, // Soil temperature in K
                         float *rshort, // Short wave rad. down in W/m^2
                         float *rlongup, // Short wave rad. up in W/m^2
                         float *sunelev); // Sun elevation

      float *pmv_hoppe_iso(float *t, // Air temperature in °C
                           float *rh, // Relative humidity in %
                           float *wind, // Wind speed in m/s
                           float *mtrad, // Mean radiant temperature K
                           float *iclo); // Clothing index

      float pmv_hoppe_iso(float t, // Air temperature in °C
                          float rh, // Relative humidity in %
                          float wind, // Wind speed in m/s
                          float mtrad, // Mean radiant temperature K
                          float iclo); // Clothing index

      float *ppd(float *pmv); // Predicted Mean Vote

      float *clomin(float *t, // Air temperature in °C
                    float *rh, // Relative humidity in %
                    float *wind, // Wind speed in m/s
                    float *mtrad); // Mean radiant temperature K

      float *tapparent(float *t, float *rh, float *wind);

      float *wellness(float *t, float *rh, float *wind);

      float *poda(float *t, float *p, float *pvap);

      float *hi(float *t,  // Air temperature in °C
                float *rh ); // Relative humidity in %

      float *humidex(float *t,  // Air temperature in °C
                     float *rh); // Relative humidity in %

      float *net(float *t,  // Air temperature in °C
                 float *rh, // Relative humidity in %
                 float *wind); // Wind speed in m/s

      float *ssi(float *t,  // Air temperature in °C
                 float *rh ); //  Relative humidity in %

      float *thom(float *t,  // Air temperature in °C
                  float *p_hPa); // Pressure in hPa (mb)

      float *wbgt(float *t,  // Air temperature in °C
                  float *rh, // Relative humidity in %
                  float *wind); // Wind speed in m/s

      float *steadman_indoor(float *t,  // Air temperature in °C
                             float *rh, // Relative humidity in %
                             float *p_hPa); // Pressure in hPa (mb)

      float *steadman_outdoor_sun(float *t,  // Air temperature in °C
                                  float *rh, // Relative humidity in %
                                  float *wind, // Wind speed in m/s
                                  float *rshort,  // SW rad. down in W/m^2
                                  float *sunelev);  // Sun elevation

      float *steadman_outdoor_shade(float *t,  // Air temperature in °C
                                    float *rh, // Relative humidity in %
                                    float *wind); // Wind speed in m/s

 
      float *utci(float *t,    // Temperature Celsius
                  float *tmrt, // Mean radiant temperature
                  float *wind, // Wind m/s
                  float *rh);  // Relative humidity %

      float *new_variable(float *t);

      // Utility functions

      std::string compass(float direction);
      float inchesHg2hPa(float p_inchesHg);
      float hPa2inchesHg(float p_hPa);
      float mmHg2hPa(float p_mmHg);
      float hPa2mmHg(float p_hPa);
      float mm2inches(float mm);
      float inches2mm(float inches);
      float FtoC(float degF);
      float CtoF(float degC);
      float mphtoms(float mph);
      float mstomph(float ms);
      float kmhtomph(float kmh);
      float mphtokmh(float mph);
      float knotstomph(float knots);
      float mphtoknots(float mph);
      float knotstoms(float knots);
      float mstoknots(float ms);
      float knotstokmh(float knots);
      float kmhtoknots(float kmh);
      float windchill(float t, float wind);
      float heatindex(float t, float rh);
      float hi_wind(float t, float w);
      float hi_class(float t, float rh, float wind);
      float sharlau_class(float t, float rh);
      float windchill_class(float t, float wind);
  };
}

#endif /*METEOSALUTE_H_*/
