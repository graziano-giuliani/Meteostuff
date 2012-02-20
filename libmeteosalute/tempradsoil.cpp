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

static inline float proj(float sunelev)
{
  static const float d2r= M_PI / 180.0;
  if (sunelev < 0.0)
    return 0.0;
  return 0.308 * cosf(d2r * (sunelev* (0.998- (powf(sunelev, 2.0) / 50000.0))));
}

float *Meteosalute::tempradsoil(float *tsoil, float *rshort, float *rshortup,
                                float *sunelev)
{
  if (vector_size <= 0)
    return 0;
  float *tempradsoil = new float[vector_size];
  float pj, tsk;
  const float sig = 5.67e-8;

  for (int i = 0; i < vector_size; i ++)
  {
    tsk = tsoil[i] + 273.12;
    pj = proj(sunelev[i]);
    tempradsoil[i] = powf(
                          273.0 + (0.97 * powf(tsk, 4.0))+ (1.0 - 0.33)
                              * rshortup[i]/ (sig* 0.97)+ (1 - 0.33) * pj
                              * ((rshort[i])/(sig* 0.97)), 0.25) - 273.16;
  }
  return tempradsoil;
}

#ifdef TEST_ME

#include <iostream>

int main(int argc, char *argv[])
{
  const int V_SIZE = 1;
  Meteosalute m(V_SIZE);
  float tsoil[V_SIZE] =
  { 285.0};
  float rshort[V_SIZE] =
  { 758.0};
  float rshortup[V_SIZE] =
  { 600.0};
  float sunelev[V_SIZE] =
  { 0.84};
  float *tr = m.tempradsoil(tsoil, rshort, rshortup, sunelev);
  std::cout << "Proj :" << proj(0.84) << std::endl;

  float *res = tr;
  for (int i = 0; i < V_SIZE; i ++)
  {
    std::cout << res[i] << std::endl;
  }
  return 0;
}

#endif
