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

#include <cmath>
#include <string>

#include <meteosalute.h>
#include <psychrometric.h>

using namespace std;
using namespace Meteo;

string Meteosalute::compass(float direction)
{
  const static char dir[16][4] = { "N", "NNE", "NE", "ENE", "E", "ESE", "SE",
                                   "SSE", "S", "SSW", "SW", "WSW", "W", "WNW",
                                   "NW", "NNW" };
  return dir[((int) ((direction+11.25)/22.5))%16];
}

float Meteosalute::inchesHg2hPa(float p_inchesHg)
{
  return (p_inchesHg * (1.0 / 33.8638864));
}

float Meteosalute::hPa2inchesHg(float p_hPa)
{
  return (p_hPa / 33.8638864);
}

float Meteosalute::mmHg2hPa(float p_mmHg)
{
  return (p_mmHg * 1.33322);
}

float Meteosalute::hPa2mmHg(float p_hPa)
{
  return (p_hPa * 0.750062);
}

float Meteosalute::mm2inches(float mm)
{
  return mm * 0.0393701;
}

float Meteosalute::inches2mm(float inches)
{
  return inches * 25.4;
}

float Meteosalute::FtoC(float degF)
{
  return degF * 0.555556- 17.7778;
}

float Meteosalute::CtoF(float degC)
{
  return degC * 1.8+ 32.0;
}

float Meteosalute::mphtoms(float mph)
{
  return mph * 0.44704;
}

float Meteosalute::mstomph(float ms)
{
  return ms * 2.23694;
}

float Meteosalute::kmhtomph(float kmh)
{
  return kmh * 0.621371;
}

float Meteosalute::mphtokmh(float mph)
{
  return mph * 1.60934;
}

float Meteosalute::knotstomph(float knots)
{
  return knots * 1.15078;
}

float Meteosalute::mphtoknots(float mph)
{
  return mph * 0.868976;
}

float Meteosalute::knotstoms(float knots)
{
  return knots * 0.514444;
}

float Meteosalute::mstoknots(float ms)
{
  return ms * 1.94384;
}

float Meteosalute::knotstokmh(float knots)
{
  return knots * 1.852;
}

float Meteosalute::kmhtoknots(float kmh)
{
  return kmh * 0.539957;
}

float Meteosalute::windchill(float t, float wind)
{
  float wkmh = wind / 0.27778;
  return (13.12 + (0.6215 * t)- (11.37 * powf(wkmh, 0.16))+(0.3965 * t)
      * powf(wkmh, 0.16));
}

inline bool inside_interval(float val, float c, float hr)
{
  return (val >= c - hr && val <= c + hr);
}

inline float veratemp(float temperature, float relhum)
{
  float irh = rintf(relhum);

  if (temperature <= 20.5)
    return rintf(-0.0003 * powf(irh, 2.0) + 0.0807 * irh + 15.944);
  else if (temperature > 20.5&& temperature <= 21.5)
    return rintf(-0.0001 * powf(irh, 2.0) + 0.0571 * irh + 17.734);
  else if (temperature > 21.5&& temperature <= 22.5)
    return rintf(0.05 * irh + 18.773);
  else if (temperature > 22.5&& temperature <= 23.5)
    return rintf(-0.0002 * powf(irh, 2.0) + 0.0698 * irh + 19.748);
  else if (temperature > 23.5&& temperature <= 24.5)
    return rintf(0.0001 * powf(irh, 2.0) + 0.0468 * irh + 21.203);
  else if (temperature > 24.5&& temperature <= 25.5)
    return rintf(0.0596 * irh + 21.986);
  else if (temperature > 25.5&& temperature <= 26.5)
    return rintf(0.0002 * powf(irh, 2.0) + 0.0452 * irh + 23.119);
  else if (temperature > 26.5&& temperature <= 27.5)
    return 27.0;
  else if (temperature > 27.5&& temperature <= 28.5)
    return 28.0;
  else
    return 29.0;

  return -9999.0;
}

float Meteosalute::heatindex(float t, float rh)
{
  float tf, tf2, ur2, hif;
  const float epsilon_t = 0.1;

  if (t < 27.0)
    return veratemp(t, rh);
  else if (inside_interval(t, 27.0, epsilon_t) && rh <= 40.0)
    return 0.0007 * (rh * rh)+ 0.0604 * rh + 23.493;
  else if ( (inside_interval(t, 28.0, epsilon_t) && rh <= 10.0)
      ||(inside_interval(t, 29.0, epsilon_t) && rh <= 6.0))
    return 26.0;
  else
  {
    tf = t * (9.0 / 5.0)+ 32.0;
    tf2 = powf(tf, 2.0);
    ur2 = powf(rh, 2.0);

    hif = -42.379 + 2.04901523 * tf + 10.1433127 * rh - 0.22475541 *tf * rh
        - 6.83783 * 0.001* tf2 - 5.481717 * 0.01* ur2 +1.22874 * 0.001* tf2* rh
        + 8.5282 * 0.0001* tf * ur2 -1.99 * 0.000001* tf2* ur2;

    return ((5.0 / 9.0) * (hif - 32.0));
  }
}

float Meteosalute::hi_wind(float t, float w)
{
  int windc;
  int wint = (int) w;

  if (wint <= 2)
    windc = 2;
  else if (wint > 2&& wint <= 6)
    windc = 6;
  else if (wint > 6&& wint <= 10)
    windc = 10;
  else if (wint > 10&& wint <= 14)
    windc = 14;
  else
    windc = 18;

  const static float hi_6[31] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                  -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
                                  -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
                                  -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
                                  -1.0 };
  const static float hi_10[31] = { -3.0, -3.0, -2.0, -2.0, -2.0, -2.0, -2.0,
                                   -2.0, -2.0, -1.0, -1.0, -1.0, -1.0, 0.0,
                                   0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0,
                                   1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0 };
  const static float hi_14[31] = { -4.0, -4.0, -3.0, -3.0, -3.0, -3.0, -3.0,
                                   -3.0, -3.0, -2.0, -2.0, -2.0, -1.0, -1.0,
                                   0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0,
                                   2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.0 };
  const static float hi_18[31] = { -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -3.0,
                                   -3.0, -3.0, -3.0, -2.0, -2.0, -1.0, -1.0,
                                   0.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0,
                                   3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 2.0 };

  if (windc == 2)
    return 0.0;
  if (t < 20.0)
    return 0.0;
  if (t > 50.0)
    return 2.0;

  if (windc == 6)
  {
    return hi_6[(int) rintf(t-20.0)];
  }
  else if (windc == 10)
  {
    return hi_10[(int) rintf(t-20.0)];
  }
  else if (windc == 14)
  {
    return hi_14[(int) rintf(t-20.0)];
  }
  return hi_18[(int) rintf(t-20.0)];
}

float Meteosalute::windchill_class(float t, float wind)
{
  float wc = windchill(t, wind);

  if (wc > 10.0)
    return 5.0;
  if (wc > 5.0&& wc <= 10.0)
    return 4.0;
  if (wc > 0.0&& wc <= 5.0)
    return 3.0;
  if (wc <= 0.0&& wc > -5.0)
    return 2.0;
  return 1.0;
}

float Meteosalute::sharlau_class(float t, float rh)
{
  float t_critic = -0.0003 * (rh*rh)+ 0.1497 * rh - 7.7133;
  float dt = (t -t_critic);

  if (dt <= -5.0)
    return 1.0;
  if (dt <= -3.0)
    return 2.0;
  if (dt <= -1.0)
    return 3.0;
  if (dt < 0.0)
    return 4.0;
  return 5.0;
}

float Meteosalute::hi_class(float t, float rh, float wind)
{
  float hi = heatindex(t, rh) + hi_wind(t, wind);

  // if (hi < 18.0) return 5.0;
  // if (hi < 21.0) return 5.0;
  // if (hi < 24.0) return 5.0;
  if (hi < 27.0)
    return 5.0;
  if (hi < 30.0)
    return 6.0;
  if (hi < 35.0)
    return 7.0;
  if (hi < 40.0)
    return 8.0;
  return 9.0;
}

float *Meteosalute::tempktoC(float *tempk)
{
  if (vector_size <= 0)
    return 0;
  float *tempc = new float[vector_size];

  for (int i = 0; i < vector_size; i ++)
  {
    tempc[i] = tempk[i] - 273.16;
  }
  return tempc;
}

float *Meteosalute::pressPatomb(float *pa)
{
  if (vector_size <= 0)
    return 0;
  float *pmb = new float[vector_size];

  for (int i = 0; i < vector_size; i ++)
  {
    pmb[i] = pa[i] * 0.01;
  }
  return pmb;
}

float *Meteosalute::pressvap(float *temp)
{
  if (vector_size <= 0)
    return 0;
  float *pv = new float[vector_size];

  for (int i = 0; i < vector_size; i ++)
  {
    pv[i] = func_t2es(temp[i]) * 10.0;
  }
  return pv;
}

float *Meteosalute::defsat(float *dewpoint, float *pvap)
{
  if (vector_size <= 0)
    return 0;
  float *df = new float[vector_size];

  for (int i = 0; i < vector_size; i ++)
  {
    df[i] = func_t2es(dewpoint[i] - 273.15) * 10.0 - pvap[i];
  }
  return df;
}

#ifdef TEST_ME

#include <iostream>

int main(int argc, char *argv[])
{
  std::cout << "0.0 : " << compass(0.0) << std::endl;
  std::cout << "90.5 : " << compass(90.5) << std::endl;
  std::cout << "185.5 : " << compass(185.5) << std::endl;
  std::cout << "220.0 : " << compass(220.0) << std::endl;
  return 0;
}

#endif

