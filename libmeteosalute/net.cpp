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

float *Meteosalute::net(float *t, float *rh, float *wind)
{
  if (vector_size <= 0) return 0;
  float *net = new float[vector_size];

  for (int i = 0; i < vector_size; i ++)
  {
    net[i] = 999.9;
    if (rh[i] > 100.1 || rh[i] < 0.0)
      continue;
    else if (wind[i] > 130.0 || wind[i] < 0.0)
      continue;
    else if (t[i] > 100.0 || t[i] < -100.0)
      continue;
    else
      net[i] = 37.0-((37.0-t[i])/(0.68-(0.0014*rh[i])) + 
        (1.0/(1.76+(1.4*pow(wind[i],0.75)))))-(0.29*t[i]*(1.0-(0.01*rh[i])));
  }
  return net;
}

#ifdef TEST_ME

#include <iostream>

int main(int argc, char *argv[])
{
  const int V_SIZE = 1;
  Meteosalute m(V_SIZE);
  float t[V_SIZE] = { 30.0};
  float rh[V_SIZE] = { 75.0};
  float wind[V_SIZE] = { 12.0};
  float *we = m.net(t, rh, wind);

  float *res = we;
  for (int i = 0; i < V_SIZE; i ++)
  {
    std::cout << res[i] << std::endl;
  }
  return 0;
}

#endif
