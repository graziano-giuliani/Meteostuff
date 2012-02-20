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

dirspeed_field *Meteosalute::dirspeed(vectorcomp_field &uv)
{
  if (vector_size <= 0)
    return 0;
  dirspeed_field *ds = new dirspeed_field;
  float *direction = new float[vector_size];
  float *speed = new float[vector_size];
  float *u = uv.u;
  float *v = uv.v;

  ds->speed = speed;
  ds->direction = direction;

  for (int i = 0; i < vector_size; i ++)
  {
    speed[i] = sqrtf(u[i] * u[i]+ v[i] * v[i]);
    direction[i] = fmodf((r2d * (atanf(u[i] / v[i]) - M_PI/2.0)- 180.0), 180.0);
  }

  return ds;
}
