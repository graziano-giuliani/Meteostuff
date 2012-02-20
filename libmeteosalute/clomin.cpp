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

float *Meteosalute::clomin(float *t, float *rh, float *wind, float *mtrad)
{
  if (vector_size <= 0)
    return 0;
  float *clomin = new float[vector_size];
  const int MAX_ITER = 20;
  const float PMV_GOOD = 0.5;
  float pmv = -1.0;

  for (int i = 0; i < vector_size; i ++)
  {
    clomin[i] = 0.1;
    for (int j = 0; j < MAX_ITER; j ++)
    {
      pmv = pmv_hoppe_iso(t[i], rh[i], wind[i], mtrad[i], clomin[i]);
      if (pmv > PMV_GOOD)
        break;
      clomin[i] += 0.2;
    }
  }
  return clomin;
}

#ifdef TEST_ME

#include <iostream>

int main(int argc, char *argv[])
{
  const int V_SIZE = 1;
  Meteosalute m(V_SIZE);
  float t[V_SIZE] =
  { 10.0};
  float rh[V_SIZE] =
  { 85.0};
  float wind[V_SIZE] =
  { 5.0};
  float mtrad[V_SIZE] =
  { 0.0};
  float *clm = m.clomin(t, rh, wind, mtrad);

  float *res = clm;
  for (int i = 0; i < V_SIZE; i ++)
  {
    std::cout << res[i] << std::endl;
  }
  return 0;
}

#endif

