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

#include <meteosalute.h>
#include <cmath>

using namespace std;
using namespace Meteo;

float *Meteosalute::hi(float *t, float *rh)
{
  if (vector_size <= 0) return 0;
  float *hi = new float[vector_size];

  for (int i = 0; i < vector_size; i ++)
  {
    hi[i] = 999.9;
    if (rh[i] > 100.1 || rh[i] < 0.0)
      continue;
    else if (t[i] > 100.0 || t[i] < -100.0)
      continue;
    else
      hi[i] = -8.784695+(1.61139411*t[i])+
                (2.338549*rh[i])-
                (0.14611605*t[i]*rh[i])-(1.2308094*pow(10,-2)*pow(t[i],2))
               -(1.6424828*pow(10,-2)*pow(rh[i],2))
               +(2.211732*pow(10,-3)*pow(t[i],2)*rh[i])
               +(7.2546*pow(10,-4)*t[i]*pow(rh[i],2))
               -(3.582*pow(10,-6)*pow(rh[i],2));
  }
  return hi;
}

#ifdef TEST_ME

#include <iostream>

int main(int argc, char *argv[])
{
  const int V_SIZE = 1;
  Meteosalute m(V_SIZE);
  float t[V_SIZE] = { 30.0};
  float rh[V_SIZE] = { 75.0};
  float *we = m.hi(t, rh);
  float *res = we;
  for (int i = 0; i < V_SIZE; i ++)
  {
    std::cout << res[i] << std::endl;
  }
  return 0;
}

#endif
