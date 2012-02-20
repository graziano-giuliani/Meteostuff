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

#include <meteosalute.h>

using namespace std;
using namespace Meteo;

float *Meteosalute::rdiffuse(float *radteoric, float *rshort)
{
  if (vector_size <= 0)
    return 0;
  float *rdiffuse = new float[vector_size];
  float kg, kd;

  for (int i = 0; i < vector_size; i ++)
  {
    if (radteoric[i] <= 0.0)
      rdiffuse[i] = 0.0;
    else
    {
      kg = rshort[i] / radteoric[i];
      kd = (0.365 - 0.25 * kg) * sinf(M_PI * kg);
      rdiffuse[i] = (kd * radteoric[i]) > 0.0 ? (kd * radteoric[i]) : 0.0;
      rdiffuse[i] = (rdiffuse[i] > rshort[i]) ? rshort[i] : rdiffuse[i];
    }
  }
  return rdiffuse;
}

#ifdef TEST_ME

#include <iostream>

int main(int argc, char *argv[])
{
  const int V_SIZE = 1;
  Meteosalute m(V_SIZE);
  float rteoric[V_SIZE] =
  { 800.0};
  float rshort[V_SIZE] =
  { 758.0};
  float *rd = m.rdiffuse(rteoric, rshort);

  float *res = rd;
  for (int i = 0; i < V_SIZE; i ++)
  {
    std::cout << res[i] << std::endl;
  }
  return 0;
}

#endif

