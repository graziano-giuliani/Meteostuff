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

#include <radargeom.h>
#include <cmath>

using namespace cetemps;

const double FourThirdEarthRadius = 8494667.0;
const double FourThirdEarthRadiusSquared = 72159367440889.00;
const double DegToRad =  0.01745329251994329576;
const double RadToDeg = 57.29577951308232087684;

radargeom::radargeom(float antenna_h_in_meter_from_ground)
{
  antenna_h = antenna_h_in_meter_from_ground;
}
 
float radargeom::beam_height(float distance_along_ray, float elevation_angle)
{
  return (antenna_h - FourThirdEarthRadius +
      sqrt(distance_along_ray*distance_along_ray+FourThirdEarthRadiusSquared+
           2*distance_along_ray*FourThirdEarthRadius*
           sin(elevation_angle*DegToRad)));
}

float radargeom::earth_distance(float distance_along_ray, float elevation_angle)
{
  float hgt = beam_height(distance_along_ray, elevation_angle);
  return sqrt(distance_along_ray*distance_along_ray - hgt*hgt);
}

float radargeom::elevation_angle(float distance, float height)
{
  double rayd = distance_along_ray(distance, height);
  double sint = (pow(height-antenna_h+FourThirdEarthRadius, 2.0)-
    rayd*rayd-FourThirdEarthRadiusSquared)/(2*rayd*FourThirdEarthRadius);
  return (asin(sint)*RadToDeg);
}

float radargeom::distance_along_ray(float distance, float height)
{
  return sqrt(distance*distance+height*height);
}

#ifdef TESTME

#include <iostream>

int main (int argc, char *argv[])
{
  radargeom g(50);

  std::cout << "Height from Earth after 25 km at elevation 2: "
    << g.beam_height(25000, 2) << " meters." << std::endl;
  std::cout << "Height from Earth after 100 km at elevation 1: "
    << g.beam_height(100000, 1) << " meters." << std::endl;
  std::cout << "Distance on Earth after 100 km on ray at elevation 1: "
    << g.earth_distance(100000, 1) << " meters." << std::endl;
  std::cout << "Distance along ray if see at 4000 m after 100 km: "
    << g.distance_along_ray(100000, 4000) << std::endl;
  std::cout << "Elevation Angle if see at 4000 m after 100 km: "
    << g.elevation_angle(100000, 4000) << std::endl;
}

#endif
