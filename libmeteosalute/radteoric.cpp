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
#include <ctime>

#include <meteosalute.h>

using namespace std;
using namespace Meteo;

float *Meteosalute::radteoric(float *lat, float *topography, float model_time)
{
  if (vector_size <= 0)
    return 0;
  float *radteoric = new float[vector_size];
  struct tm *tpnt;
  time_t itime;
  float jday, hday, dr, soldec, factor, rt;

  itime = (time_t) model_time;
  tpnt = gmtime(&itime);
  jday = tpnt->tm_yday;
  hday = tpnt->tm_hour;
  dr = 1.0 + 0.033 * cosf(TWO_M_PI * (jday / 365.0));
  soldec = 23.45 * d2r * cos(TWO_M_PI / 365.0* (172 - fmodf(jday, 365.0)));

  for (int i = 0; i < vector_size; i ++)
  {
    factor = sinf(lat[i] * d2r) * sinf(soldec) +cosf(lat[i] * d2r)
        * cosf(soldec)* cosf((hday * 15.0) * d2r);
    rt = (1360.0 * dr) * factor * -1.0;
    if (rt < 0.0)
      rt = 0.0;
    radteoric[i] = (0.75 + 2.0 * powf(10.0, -5.0) * topography[i]) * rt;
  }
  return radteoric;
}

#ifdef TEST_ME

#include <iostream>

int main(int argc, char *argv[])
{
  const int V_SIZE = 7;
  Meteosalute m(V_SIZE);
  float utime = 1191680523.0; // Sat Oct  6 14:22:03 GMT 2007
  float lat[] =
  { -90.0, -45.0, -10.0, 0.0, 10.0, 45.0, 90.0};
  float topo[] =
  { 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0};
  float *se = m.radteoric(lat, topo, utime);

  float *res = se;
  for (int i = 0; i < V_SIZE; i ++)
  {
    std::cout << lat[i] << ", " << topo[i] << ", "
    << (long) utime << " : " << res[i] << std::endl;
  }
  float vlat[] =
  { 45.0, 45.0, 45.0, 45.0, 45.0, 45.0, 45.0};
  float vtopo[] =
  { 10.0, 100.0, 150.0, 250.0, 500.0, 1000.0, 2000.0};
  delete [] se;
  se = m.radteoric(vlat, vtopo, utime);
  res = se;
  for (int i = 0; i < V_SIZE; i ++)
  {
    std::cout << lat[i] << ", " << topo[i] << ", "
    << (long) utime << " : " << res[i] << std::endl;
  }

  return 0;
}

#endif
