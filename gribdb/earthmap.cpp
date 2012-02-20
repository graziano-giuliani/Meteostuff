/***************************************************************************
 *   Copyright (C) 2008-2010 by Graziano Giuliani                          *
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

#include <cmath>
#include <earthmap.h>

#ifdef TESTME
#include <iostream>
#endif

using namespace himet;

// Constants using Libm Pi

static const double rad2deg        = 180.0/M_PI;
static const double deg2rad        = M_PI/180.0;
static const double earth_radius_m = 6371200.0;
static const double fill           = 1.0e+21;

// Create object

earthmap::earthmap( ) 
{
  this->init = false;
}

// Set all map projection variables

void earthmap::set( t_proj_latlon_parameters &parms )
{
  this->init = false;

  // Argument check

  if (fabs(parms.lat1) > 90.0001 ||
      (parms.nx < 0 || parms.ny < 0))
    return;

  // Grid dimension
  nx   = parms.nx;
  ny   = parms.ny;

  // Latitude longitude of (1,1) and (NX,NY)
  lat1 = parms.lat1;
  lon1 = parms.lon1;
  lat2 = parms.lat2;
  lon2 = parms.lon2;

  // Parse scan flag (GRIB1 && 2)
  iscan = parms.scanflag / 128 % 2;
  jscan = parms.scanflag /  64 % 2;
  nscan = parms.scanflag /  32 % 2;
  oscan = parms.scanflag /  16 % 2;
 
  // Multiplier for increments
  hi = pow(-1.0, iscan); // Default is +i dir (WE)
  hj = pow(-1.0, jscan); // Default is -j dir (NS)

  // Increments on axes in latitude and longitude
  double temp = hi * (lon2 - lon1) + 3599;
  dlon = hi * (fmod(temp, 360.0) + 1) / (nx - 1);
  dlat = hj * (lat2 - lat1) / (ny - 1);

  // Extrema of indexes
  xmin = 0.0;
  temp = 360.0 / fabs(dlon);
  if (nx == (int) rintf(temp))
    xmax = (double) (nx + 2);
  else
    xmax = (double) (nx + 1);
  ymin = 0.0;
  ymax = (double) (ny + 1);
 
  // Set up callbacks
  llij = &himet::earthmap::llij_latlon;
  ijll = &himet::earthmap::ijll_latlon;
  gluv = &himet::earthmap::grid2lluv_latlon;
  lguv = &himet::earthmap::ll2griduv_latlon;

  // All ok
  this->init = true;
  return;
}

void earthmap::set( t_proj_lambert_parameters &parms )
{
  this->init = false;

  // Argument check

  if (fabs(parms.lat1)     > 90.0001 ||
      fabs(parms.truelat1) > 90.0    ||
      fabs(parms.truelat2) > 90.0    ||
      fabs(parms.dx)      <= 0.0     ||
      fabs(parms.dy)      <= 0.0     ||
      (parms.nx < 0 || parms.ny < 0))
    return;

  // Grid dimension
  nx       = parms.nx;
  ny       = parms.ny;

  // Projection parameters
  lat1     = parms.lat1;
  lon1     = parms.lon1;
  orient   = parms.orient;
  truelat1 = parms.truelat1;
  truelat2 = parms.truelat2;
  iproj    = parms.projflag / 128 % 2; // Hemisphere flag

  irot     = parms.resflag / 8 % 2; // wind earth or grid relative ?

  // Resolution of grid in meters
  dx       = parms.dx;
  dy       = parms.dy;

  // Parse scan flag (GRIB1 && 2)
  iscan = parms.scanflag / 128 % 2;
  jscan = parms.scanflag /  64 % 2;
  nscan = parms.scanflag /  32 % 2;
  oscan = parms.scanflag /  16 % 2;

  // Multiplier for increments: Expected default (+i,+j) == scanflag bit 2 set
  hi = pow(-1.0, iscan);
  hj = pow(-1.0, (1 - jscan)); // Bit 2 should be set

  // Pole position (1 North pole, 0 South pole)
  pp = pow(-1.0, iproj);
 
  // Increments with direction
  dxs = dx * hi;
  dys = dy * hj;

  // Calculate cone factor
  if (truelat1 == truelat2)
    cone = sin(pp * truelat1 * deg2rad);
  else
    cone = log(cos(truelat1 * deg2rad) / cos(truelat2 * deg2rad)) /
         log(tan((pp * 90.0 - truelat1) * 0.5 * deg2rad) /
         tan((pp * 90.0 - truelat2) * 0.5 * deg2rad));

  // Calculate earth relative distances
  double d1 = (double) tan((pp * truelat1 + 90.0) * 0.5 * deg2rad);
  de = cos(truelat1 * deg2rad) * earth_radius_m * pow(d1, cone) / cone;

  // Special case for pole
  if (pp * lat1 == 90.0)
  {
    xp = 1.0;
    yp = 1.0;
  }
  else
  {
    d1 = (double) tan((pp * lat1 + 90) * 0.5 * deg2rad);
    double dr = de / pow(d1, cone);
    double temp1 = lon1 - orient + 3780.0;
    double difflon = fmod(temp1, 360.0) - 180.0;
    xp = -pp * sin(cone * difflon * deg2rad) * dr / dxs;
    yp = cos(cone * difflon * deg2rad) * dr / dys;
  }

  // Utility factors for calculation
  invcone2 = 1.0 / (cone * 2.0);
  de2 = de * de;

  // Extrema
  xmin = 0.0;
  xmax = (double) (nx + 1);
  ymin = 0.0;
  ymax = (double) (ny + 1);

  // Setup callbacks
  llij = &himet::earthmap::llij_lc;
  ijll = &himet::earthmap::ijll_lc;
  gluv = &himet::earthmap::grid2lluv_lc;
  lguv = &himet::earthmap::ll2griduv_lc;

  // All Ok
  this->init = true;
  return;
}

void earthmap::set( t_proj_polar_parameters &parms )
{
  this->init = false;

  // Argument check

  if (fabs(parms.lat1)    > 90.0001 ||
      fabs(parms.dx)     <= 0.0     ||
      fabs(parms.dy)     <= 0.0     ||
      (parms.nx < 0 || parms.ny < 0))
    return;

  // Grid dimension
  nx       = parms.nx;
  ny       = parms.ny;

  // Projection parameters
  lat1     = parms.lat1;
  lon1     = parms.lon1;
  orient   = parms.orient;
  iproj    = parms.projflag / 128 % 2; // Hemisphere flag

  irot     = parms.resflag / 8 % 2; // wind earth or grid relative ?

  // Grid resolution in meters
  dx       = parms.dx;
  dy       = parms.dy;

  // Parse scan flag (GRIB1 && 2)
  iscan = parms.scanflag / 128 % 2;
  jscan = parms.scanflag /  64 % 2;
  nscan = parms.scanflag /  32 % 2;
  oscan = parms.scanflag /  16 % 2;

  // Multiplier for increments: Expected default (+i,+j) == scanflag bit 2 set
  hi = pow(-1.0, iscan);
  hj = pow(-1.0, (1 - jscan));

  // Pole position (1 North pole, 0 South pole)
  pp = pow(-1.0, iproj);
 
  // Increments with direction
  dxs = dx * hi;
  dys = dy * hj;

  // Calculate earth relative distances and pole position
  de = (sin(60.0 * deg2rad) + 1.0) * earth_radius_m;
  double dr = de * cos(lat1 * deg2rad) / (pp * sin(lat1 * deg2rad) + 1.0);
  xp = 0.0 - pp * sin((lon1 - orient) * deg2rad) * dr / dxs;
  yp = 0.0 + cos((lon1 - orient) * deg2rad) * dr / dys;

  // Utility constants in calculation
  de2 = de * de;

  // Extrema for indexes
  xmin = 0.0;
  xmax = (double) (nx + 1);
  ymin = 0.0;
  ymax = (double) (ny + 1);

  // Setup callbacks
  llij = &himet::earthmap::llij_ps;
  ijll = &himet::earthmap::ijll_ps;
  gluv = &himet::earthmap::grid2lluv_ps;
  lguv = &himet::earthmap::ll2griduv_ps;

  // All Ok
  this->init = true;
  return;
}

void earthmap::set( t_proj_mercator_parameters &parms )
{
  this->init = false;

  // Argument check

  if (fabs(parms.lat1)   > 90.0001 ||
      fabs(parms.reflat) > 90.0    ||
      fabs(parms.dx)    <= 0.0     ||
      fabs(parms.dy)    <= 0.0     ||
      (parms.nx < 0 || parms.ny < 0))
    return;

  // Grid dimension
  nx       = parms.nx;
  ny       = parms.ny;

  // Projection parameters
  lat1     = parms.lat1;
  lon1     = parms.lon1;
  lat2     = parms.lat2;
  lon2     = parms.lon2;
  truelat1 = parms.reflat;

  // Grid resolution in meters
  dx       = parms.dx;
  dy       = parms.dy;

  // Parse scan flag (GRIB1 && 2)
  iscan = parms.scanflag / 128 % 2;
  jscan = parms.scanflag /  64 % 2;
  nscan = parms.scanflag /  32 % 2;
  oscan = parms.scanflag /  16 % 2;
 
  // Multiplier for increments: Expected default (+i,+j) == scanflag bit 2 set
  hi = pow(-1.0, iscan);
  hj = pow(-1.0, (1 - jscan));

  // Setup increments on lat/lon
  double temp1 = hi * (lon2 - lon1) + 3599;
  dlon = hi * (fmod(temp1, 360.0) + 1.0) / (nx - 1.0);
  dlat = hj * dy / (cos(truelat1 * deg2rad) * earth_radius_m);

  // Mercator factor
  ye = -log(tan((lat1 + 90.0) * 0.5 * deg2rad)) / dlat;

  // Extrema for indexes
  xmin = 0.0;
  temp1 = 360.0 / fabs(dlon);
  if (nx == (int) rintf(temp1))
    xmax = (double) (nx + 2);
  else
    xmax = (double) (nx + 1);
  ymin = 0.0;
  ymax = (double) (ny + 1);

  // Setup callbacks
  llij = &himet::earthmap::llij_merc;
  ijll = &himet::earthmap::ijll_merc;
  gluv = &himet::earthmap::grid2lluv_merc;
  lguv = &himet::earthmap::ll2griduv_merc;

  // All Ok
  this->init = true;
  return;
}

void earthmap::set( t_proj_rotlat_parameters &parms )
{
  this->init = false;

  // Argument check

  if (fabs(parms.lat1)    > 90.0001 ||
      fabs(parms.polelat) > 90.0    ||
      (parms.nx < 0 || parms.ny < 0))
    return;

  // Grid dimension
  nx       = parms.nx;
  ny       = parms.ny;

  // Projection parameters
  lat1     = parms.lat1;
  lon1     = parms.lon1;
  lat2     = parms.lat2;
  lon2     = parms.lon2;

  // Derived from parameters and used in calculations
  lat00    = ( 90.0 + parms.polelat) * deg2rad;
  lon00    = parms.polelon  * deg2rad;
  rotation = parms.rotation * deg2rad;

  // Parse scan flag (GRIB1 && 2)
  iscan = parms.scanflag / 128 % 2;
  jscan = parms.scanflag /  64 % 2;
  nscan = parms.scanflag /  32 % 2;
  oscan = parms.scanflag /  16 % 2;

  // Multiplier for increments: Expected default (+i,+j) == scanflag bit 2 set
  hi = pow(-1.0, iscan);
  hj = pow(-1.0, (1 - jscan));

  // Increments on rotated grid
  dlon = hi * (lon2 - lon1) / (nx - 1);
  dlat = hj * (lat2 - lat1) / (ny - 1);

  // Extrema for indexes
  xmin = 0.0;
  double temp = 360.0 / fabs(dlon);
  if (nx == (int) rintf(temp))
    xmax = (double) (nx + 2);
  else
    xmax = (double) (nx + 1);
  ymin = 0.0;
  ymax = (double) (ny + 1);

  // Setup callbacks
  llij = &himet::earthmap::llij_rl;
  ijll = &himet::earthmap::ijll_rl;
  gluv = &himet::earthmap::grid2lluv_rl;
  lguv = &himet::earthmap::ll2griduv_rl;

  // All Ok
  this->init = true;
  return;
}

void earthmap::set( t_proj_rotwrf_parameters &parms )
{
  this->init = false;

  // Argument check

  if (fabs(parms.lat1) > 90.0001 ||
      fabs(parms.lat0) > 90.0    ||
      (parms.nx < 0 || parms.ny < 0))
    return;

  // Grid dimension
  nx       = parms.nx;
  ny       = parms.ny;

  // Projection parameters
  lat0     = parms.lat0;
  lon0     = parms.lon0;
  lat00    = lat0 * deg2rad;
  lon00    = lon0 * deg2rad;

  // Calculate first point coordinates on true grid
  double slat1 = sin(parms.lat1*deg2rad);
  double clat1 = cos(parms.lat1*deg2rad);
  double slat0 = sin(lat00);
  double clat0 = cos(lat00);
  double sign = copysign(1.0, fmod(parms.lon1-lon0+180+3600,360.0)-180);
  double clon1 = cos((parms.lon1-lon0) * deg2rad);
  double slatr = clat0*slat1-slat0*clat1*clon1;
  double clatr = sqrt(1.0-slatr*slatr);
  double clonr = (clat0*clat1*clon1+slat0*slat1)/clatr;
  lat1 = asin(slatr) * rad2deg;
  lon1 = sign * acos(clonr) * rad2deg;

  // Parse scan flag (GRIB1 && 2)
  kscan = parms.scanflag / 256 % 2;  // Documented in NCEP for staggering
  iscan = parms.scanflag / 128 % 2;
  jscan = parms.scanflag /  64 % 2;
  nscan = parms.scanflag /  32 % 2;
  oscan = parms.scanflag /  16 % 2;

  // Multiplier for increments: Expected default (+i,+j) == scanflag bit 2 set
  hi = pow(-1.0, iscan);
  hj = pow(-1.0, (1 - jscan));

  // Increments with direction on rotated grid
  dy = hi * lat1/(-(ny-1)/3*0.5);
  dx = hj * lon1/(-(nx-1)*0.5);

  // Extrema for indexes
  xmin = 0.0;
  double xdlon = 360.0 / fabs(dlons);
  if (nx == rintf(xdlon))
    xmax = (double) (nx + 2);
  else
    xmax = (double) (nx + 1);
  ymin = 0.0;
  ymax = (double) (ny + 1);

  // Setup callbacks
  llij = &himet::earthmap::llij_rwrf;
  ijll = &himet::earthmap::ijll_rwrf;
  gluv = &himet::earthmap::grid2lluv_rwrf;
  lguv = &himet::earthmap::ll2griduv_rwrf;

  // All Ok
  this->init = true;
  return;
}

// Converts input lat/lon values to the cartesian (i,j) value

void earthmap::latlon_to_ij( double lat, double lon, double *i, double *j )
{
  if (! this->init) throw;
  (*this.*llij) (lat, lon, i, j);
  if (fabs(*i) < 1e-12) *i = 0.0;
  if (fabs(*j) < 1e-12) *j = 0.0;
  return;
}

// Computes geographical latitude and longitude for a given (i,j) point

void earthmap::ij_to_latlon( double i, double j, double *lat, double *lon )
{
  if (! this->init) throw;
  (*this.*ijll) (i, j, lat, lon);
  return;
}

//------------------------------------
// Actual calculations
//------------------------------------

// Compute the i/j point from latitude and longitude for Polar Stereo Proj

void earthmap::llij_ps(double lat, double lon, double *i, double *j)
{
  if (fabs(lat) <= 90.0 && fabs(lon) <= 360.0 && pp * lat != -90.0)
  {
    double dr = de * tan((90.0 - pp * lat) * 0.5 * deg2rad);
    double ri, rj;
    ri = xp + pp * sin((lon - orient) * deg2rad) * dr / dxs;
    rj = yp - cos((lon - orient) * deg2rad) * dr / dys;
    if (ri < 0 || ri >= nx || rj < 0 || rj >= ny)
    {
      ri = fill;
      rj = fill;
    }
    *i = ri;
    *j = rj;
  }
  else
  {
    *i = fill;
    *j = fill;
  }
  return;
}

// Compute the latitude and longitude of an i/j point for Polar Stereo Proj

void earthmap::ijll_ps(double i, double j, double *lat, double *lon)
{
  if (i >= xmin && i <= xmax && j >= ymin && j <= ymax)
  {
    double di = (i - xp) * dxs;
    double dj = (j - yp) * dys;
    double dr2 = di*di + dj*dj;
    if (dr2 < de2 * 1e-6)
    {
      *lon = 0.0;
      *lat = pp * 90.0;
    }
    else
    {
      double temp1 = orient + pp * rad2deg * atan2(di, -dj) + 3600;
      *lon = fmod(temp1, 360.0);
      *lat = pp * rad2deg * asin((de2 - dr2) / (de2 + dr2));
    }
  }
  else
  {
    *lon = fill;
    *lat = fill;
  }
  return;
}

// Compute the geographical latitude and longitude values
// to the cartesian x/y on a Lambert Conformal projection.

void earthmap::llij_lc(double lat, double lon, double *i, double *j)
{
  if (fabs(lat) <= 90.0 && fabs(lon) <= 360.0 && pp * lat != -90.0)
  {
    double d1 = tan((90.0 - pp * lat) * 0.5 * deg2rad);
    double dr = de * pow(d1, cone);
    double temp1 = lon - orient + 3780;
    double difflon = fmod(temp1, 360.0) - 180.0;
    double ri, rj;
    ri = xp + pp * sin(cone * difflon * deg2rad) * dr / dxs;
    rj = yp - cos(cone * difflon * deg2rad) * dr / dys;
    if (iproj == 1)
    {
      ri = 2 - ri;
      rj = 2 - rj;
    }
    if (ri < 0 || ri >= nx || rj < 0 || rj >= ny)
    {
      ri = fill;
      rj = fill;
    }
    *i = ri;
    *j = rj;
  }
  else
  {
    *i = fill;
    *j = fill;
  }

  return;
}

// Convert from the (i,j) cartesian coordinate to the
// geographical latitude and longitude for a Lambert Conformal projection.

void earthmap::ijll_lc(double i, double j, double *lat, double *lon)
{
  if (i >= xmin && i <= xmax && j >= ymin && j <= ymax)
  {
    double xpts, ypts;
    if (iproj == 1)
    {
      xpts = 2 - i;
      ypts = 2 - j;
    }
    else
    {
      xpts = i;
      ypts = j;
    }
    double di, dj, dr2;
    di = (xpts - xp) * dxs;
    dj = (ypts - yp) * dys;
    dr2 = di*di + dj*dj;
    if (dr2 < de2 * 1.0e-6)
    {
      *lon = 0.0;
      *lat = pp * 90.0;
    }
    else
    {
      double temp1 = orient + pp / cone * rad2deg * atan2(di, -dj) + 3600.0;
      *lon = fmod(temp1, 360.0);
      temp1 = (de2 / dr2);
      *lat = pp * (atan(pow(temp1, invcone2)) * 2.0 * rad2deg - 90.0);
    }
  }
  else
  {
    *lon = fill;
    *lat = fill;
  }
  return;
}

// Compute i/j coordinate from lat lon for mercator projection

void earthmap::llij_merc(double lat, double lon, double *i, double *j)
{
  if (fabs(lat) <= 90.0 && fabs(lon) <= 360.0)
  {
    double temp1 = hi * (lon - lon1) + 3600;
    double ri, rj;
    ri = hi * fmod(temp1, 360.0) / dlon;
    rj = ye + log(tan((lat + 90.0) * 0.5 * deg2rad)) / dlat;
    if (ri < 0 || ri >= nx || rj < 0 || rj >= ny)
    {
      ri = fill;
      rj = fill;
    }
    *i = ri;
    *j = rj;
  }
  else
  {
    *i = fill;
    *j = fill;
  }
  return;
}

// Compute the lat/lon from i/j for mercator projection

void earthmap::ijll_merc(double i, double j, double *lat, double *lon)
{
  if (i >= xmin && i <= xmax && j >= ymin && j <= ymax)
  {
    double temp1 = lon1 + dlon * i + 3600;
    *lon = fmod(temp1, 360.0);
    *lat = atan(exp(dlat * (j - ye))) * 2.0 * rad2deg - 90.0;
  }
  else
  {
    *lon = fill;
    *lat = fill;
  }
  return;
}

// Compute the i/j location of a lat/lon on a LATLON grid.

void earthmap::llij_latlon(double lat, double lon, double *i, double *j)
{
  if (fabs(lat) <= 90.0 && fabs(lon) <= 360.0)
  {
    double temp = hi * (lon - lon1) + 3600;
    double ri, rj;
    ri = hi * fmod(temp, 360.0) / dlon;
    rj = (lat - lat1) / dlat;
    if (ri < 0 || ri >= nx || rj < 0 || rj >= ny)
    {
      ri = fill;
      rj = fill;
    }
    *i = ri;
    *j = rj;
  }
  else
  {
    *i = fill;
    *j = fill;
  }
  return;
}

// Compute the lat/lon location of a i/j on a LATLON grid.

void earthmap::ijll_latlon(double i, double j, double *lat, double *lon)
{
  if (i >= xmin && i <= xmax && j >= ymin && j <= ymax)
  {
    double temp1, temp2;

    temp1 = lon1 + dlon * i + 3600;
    *lon = fmod(temp1, 360.0);
    temp2 = lat1 + dlat * j;
    temp1 = temp2 > -90.0 ? temp2 : -90.0;
    *lat = temp1 < 90.0 ? temp1 : 90.0;
  }
  else
  {
    *lon = fill;
    *lat = fill;
  }
  return;
}

// Compute the i/j location of a lat/lon on a Rotated LATLON grid.

void earthmap::llij_rl(double lat, double lon, double *i, double *j)
{
  if (fabs(lat) <= 90.0 && fabs(lon) <= 360.0)
  {
    double p1, g1, p, g, rlat, rlon;

    p1 = deg2rad * lat;
    g1 = -deg2rad * lon;

    p = asin( sin(lat00) * sin(p1) + cos(lat00) * cos(p1) * cos (g1 + lon00) );
    g = rotation + atan2(cos(p1) * sin (g1 + lon00),
          sin(lat00) * cos(p1) * cos (g1 + lon00) - cos(lat00) * sin(p1) );

    rlat = rad2deg * asin( -cos(p) * cos (g) );
    rlon = -rad2deg * atan2(cos(p) * sin (g), sin(p) );

    double ri, rj;
    ri = (rlon - lon1) / dlon;
    rj = (rlat - lat1) / dlat;
    if (ri < 0 || ri >= nx || rj < 0 || rj >= ny)
    {
      ri = fill;
      rj = fill;
    }
    *i = ri;
    *j = rj;
  }
  else
  {
    *i = fill;
    *j = fill;
  }
  return;
}

// Compute the lat/lon location of a i/j on a Rotated LATLON grid.

void earthmap::ijll_rl(double i, double j, double *lat, double *lon)
{
  if (i >= xmin && i <= xmax && j >= ymin && j <= ymax)
  {
    double rlat, rlon;
    double pr, gr, pm, gm;

    rlon = lon1 + i * dlon;
    rlat = lat1 + j * dlat;

    pr = deg2rad * rlat;
    gr = -deg2rad * rlon;
    pm = asin( cos(pr) * cos (gr) );
    gm = atan2(cos(pr) * sin (gr), -sin(pr) );

    *lat = rad2deg * asin( sin(lat00) * sin(pm) -
           cos(lat00) * cos(pm) * cos (gm - rotation) );
    *lon = -rad2deg * (-lon00 + atan2(cos(pm) * sin (gm - rotation),
           sin(lat00) * cos(pm) * cos (gm - rotation) + cos(lat00) * sin(pm)));
  }
  else
  {
    *lon = fill;
    *lat = fill;
  }
  return;
}

// Compute the i/j location of a lat/lon on a WRF Rotated LATLON grid.

void earthmap::llij_rwrf(double lat, double lon, double *i, double *j)
{
  if (fabs(lat) <= 90.0 && fabs(lon) <= 360.0)
  {
    double p1, g1, p, g;
    double talon, talat, lon360;

    // To account for issues around the dateline, convert the incoming
    // longitudes to be 0->360.
    if (lon < 0.0) lon360 = lon + 360.0;
    else lon360 = lon;

    p1 = deg2rad * lat;
    g1 = -deg2rad * lon360;

    p = asin( sin(lat00) * sin(p1) + cos(lat00) * cos(p1) * cos(g1 + lon00) );
    g = atan2(cos(p1) * sin(g1 + lon00),
        sin(lat00) * cos(p1) * cos(g1 + lon00) - cos(lat00) * sin(p1) );

    talat = rad2deg * asin( -cos(p) * cos(g) );
    talon = -rad2deg * atan2(cos(p) * sin(g), sin(p) );

    double ri, rj;
    ri = (talon - lon1) / dx;
    rj = (talat - lat1) / dy;
    if (ri < 0 || ri >= nx || rj < 0 || rj >= ny)
    {
      ri = fill;
      rj = fill;
    }
    *i = ri;
    *j = rj;
  }
  else
  {
    *i = fill;
    *j = fill;
  }
  return;
}

// Compute the lat/lon location of a i/j on a WRF Rotated LATLON grid.

void earthmap::ijll_rwrf(double i, double j, double *lat, double *lon)
{
  if (i >= xmin && i <= xmax && j >= ymin && j <= ymax)
  {
    double talat, talon;
    double pr, gr, pm, gm;

    talat = lat1 + j*dy;
    talon = lon1 + i*dx;

    pr = deg2rad * talat;
    gr = -deg2rad * talon;

    pm = asin( cos(pr) * cos (gr) );
    gm = atan2(cos(pr) * sin (gr), -sin(pr) );

    *lat = rad2deg * asin( sin(lat00) * sin(pm) - cos(lat00) *
          cos(pm) * cos (gm) );
    *lon = -rad2deg * (-lon00 + atan2(cos(pm) * sin (gm),
           sin(lat00) * cos(pm) * cos (gm) + cos(lat00) * sin(pm)));
    if ((fabs(*lat) > 90.0) || (fabs(i*dx) > 360.0))
    {
      *lat = fill;
      *lon = fill;
    }
    else
    {
      *lon = *lon + 360.0;
      *lon = fmod(*lon, 360.0);
      if (*lon > 180.0) *lon = *lon - 360.0;
    }
  }
  else
  {
    *lon = fill;
    *lat = fill;
  }
  return;
}

// Convert a wind from grid north to true north.

void earthmap::gridwind_to_truewind(double i, double j,
                                    double ugrid, double vgrid,
                                    double *utrue, double *vtrue)
{
  if (! this->init) throw;
  double crot, srot;
  (*this.*gluv) (i, j, &crot, &srot);
  *utrue = crot * ugrid + srot * vgrid;
  *vtrue = -srot * ugrid + crot * vgrid;
  return;
}

// Compute grid-relative u/v wind components from the earth-relative
// values for a given projection.

void earthmap::truewind_to_gridwind(double lat, double lon,
                                    double utrue, double vtrue,
                                    double *ugrid, double *vgrid)
{
  if (! this->init) throw;
  double crot, srot;
  (*this.*lguv) (lat, lon, &crot, &srot);
  *ugrid = crot * utrue + srot * vtrue;
  *vgrid = -srot * utrue + crot * vtrue;
  return;
}

void earthmap::grid2lluv_latlon(double i, double j,
                                double *crot, double *srot)
{
  *crot = 1.0;
  *srot = 0.0;
  return;
}

void earthmap::ll2griduv_latlon(double lat, double lon,
                                double *crot, double *srot)
{
  *crot = 1.0;
  *srot = 0.0;
  return;
}

void earthmap::grid2lluv_lc(double i, double j,
                            double *crot, double *srot)
{
  if (irot == 1)
  {
    double lat, lon;
    ijll_lc(i, j, &lat, &lon);
    double temp1 = lon - orient + 3780;
    double difflon = fmod(temp1, 360.0) - 180.0;
    double alpha = cone * difflon * deg2rad;
    *crot = pp * cos(alpha);
    *srot = sin(alpha);
  }
  else
  {
    *crot = 1.0;
    *srot = 0.0;
  }
  return;
}

void earthmap::ll2griduv_lc(double lat, double lon,
                            double *crot, double *srot)
{
  if (irot == 1)
  {
    double temp1 = orient - lon + 3780;
    double difflon = fmod(temp1, 360.0) - 180.0;
    double alpha = cone * difflon * deg2rad;
    *crot = pp * cos(alpha);
    *srot = sin(alpha);
  }
  else
  {
    *crot = 1.0;
    *srot = 0.0;
  }
  return;
}

void earthmap::grid2lluv_ps(double i, double j,
                            double *crot, double *srot)
{
  if (irot == 1)
  {
    double lat, lon;
    ijll_ps(i, j, &lat, &lon);
    double temp1 = lon - orient + 3780;
    double difflon = fmod(temp1, 360.0) - 180.0;
    double alpha = cone * difflon * deg2rad;
    *crot = pp * cos(alpha);
    *srot = sin(alpha);
  }
  else
  {
    *crot = 1.0;
    *srot = 0.0;
  }
  return;
}

void earthmap::ll2griduv_ps(double lat, double lon,
                            double *crot, double *srot)
{
  if (irot == 1)
  {
    double temp1 = orient - lon + 3780;
    double difflon = fmod(temp1, 360.0) - 180.0;
    double alpha = cone * difflon * deg2rad;
    *crot = pp * cos(alpha);
    *srot = sin(alpha);
  }
  else
  {
    *crot = 1.0;
    *srot = 0.0;
  }
  return;
}

void earthmap::grid2lluv_merc(double i, double j,
                              double *crot, double *srot)
{
  *crot = 1.0;
  *srot = 0.0;
  return;
}

void earthmap::ll2griduv_merc(double lat, double lon,
                              double *crot, double *srot)
{
  *crot = 1.0;
  *srot = 0.0;
  return;
}

void earthmap::grid2lluv_rl(double i, double j,
                            double *crot, double *srot)
{
  if (irot == 1)
  {
    double lat, lon;
    ijll_rl(i, j, &lat, &lon);
    double phi0 = lat1 * deg2rad;
    double temp1 = lon0 - lon + 3780;
    double difflon = fmod(temp1, 360.0) - 180.0;
    double relm = difflon * deg2rad;
    double cos_phi0 = cos(phi0);
    double sin_phi0 = sin(phi0);
    double rlat = lat * deg2rad;
    double bigd = cos(asin(cos_phi0 * sin(rlat) -
                    sin_phi0 * cos(rlat) * cos(relm)));
    *crot = (cos_phi0 * cos(rlat) + sin_phi0 * sin(rlat) * cos(relm)) / bigd;
    *srot = sin_phi0 * sin(relm) / bigd;
  }
  else
  {
    *crot = 1.0;
    *srot = 0.0;
  }
  return;
}

void earthmap::ll2griduv_rl(double lat, double lon,
                            double *crot, double *srot)
{
  if (irot == 1)
  {
    double phi0 = lat1 * deg2rad;
    double temp1 = lon - lon0 + 3780;
    double difflon = fmod(temp1, 360.0) - 180.0;
    double relm = difflon * deg2rad;
    double cos_phi0 = cos(phi0);
    double sin_phi0 = sin(phi0);
    double rlat = lat * deg2rad;
    double bigd = cos(asin(cos_phi0 * sin(rlat) -
                    sin_phi0 * cos(rlat) * cos(relm)));
    *crot = (cos_phi0 * cos(rlat) + sin_phi0 * sin(rlat) * cos(relm)) / bigd;
    *srot = sin_phi0 * sin(relm) / bigd;
  }
  else
  {
    *crot = 1.0;
    *srot = 0.0;
  }
  return;
}

void earthmap::grid2lluv_rwrf(double i, double j,
                              double *crot, double *srot)
{
  if (irot == 1)
  {
    double lat, lon;
    ijll_rwrf(i, j, &lat, &lon);
    double phi0 = lat1 * deg2rad;
    double temp1 = lon0 - lon + 3780;
    double difflon = fmod(temp1, 360.0) - 180.0;
    double relm = difflon * deg2rad;
    double cos_phi0 = cos(phi0);
    double sin_phi0 = sin(phi0);
    double rlat = lat * deg2rad;
    double bigd = cos(asin(cos_phi0 * sin(rlat) -
                    sin_phi0 * cos(rlat) * cos(relm)));
    *crot = (cos_phi0 * cos(rlat) + sin_phi0 * sin(rlat) * cos(relm)) / bigd;
    *srot = sin_phi0 * sin(relm) / bigd;
  }
  else
  {
    *crot = 1.0;
    *srot = 0.0;
  }
  return;
}

void earthmap::ll2griduv_rwrf(double lat, double lon,
                              double *crot, double *srot)
{
  if (irot == 1)
  {
    double phi0 = lat1 * deg2rad;
    double temp1 = lon - lon0 + 3780;
    double difflon = fmod(temp1, 360.0) - 180.0;
    double relm = difflon * deg2rad;
    double cos_phi0 = cos(phi0);
    double sin_phi0 = sin(phi0);
    double rlat = lat * deg2rad;
    double bigd = cos(asin(cos_phi0 * sin(rlat) -
                    sin_phi0 * cos(rlat) * cos(relm)));
    *crot = (cos_phi0 * cos(rlat) + sin_phi0 * sin(rlat) * cos(relm)) / bigd;
    *srot = sin_phi0 * sin(relm) / bigd;
  }
  else
  {
    *crot = 1.0;
    *srot = 0.0;
  }
  return;
}

void earthmap::rotuv_ll2gr(double lat, double lon, double *crot, double *srot)
{
  if (! this->init) throw;
  (*this.*lguv) (lat, lon, crot, srot);
  return;
}

void earthmap::rotuv_gr2ll(double i, double j, double *crot, double *srot)
{
  if (! this->init) throw;
  (*this.*gluv) (i, j, crot, srot);
  return;
}

void earthmap::ij2pl(double i, double j, double *pixel, double *line)
{
  if (! this->init) throw;
  if (i < 0 || j < 0 || i >= nx || j >= ny)
  {
    *pixel = fill;
    *line = fill;
  }
  *pixel = i;
  *line = j;
  if (oscan)
    *pixel = nx - 1 - i;
  if (nscan)
  {
    double tmp = *pixel;
    *pixel = *line;
    *line = tmp;
  }
  return;
}

int earthmap::ij2p(int i, int j)
{
  if (! this->init) throw;
  if (i < 0 || j < 0 || i >= nx || j >= ny) return -1;
  // j = jscan ? j : ny-1 - j; // In ll2ij
  // i = (iscan) ? nx-1 - i : i; // in ll2ij
  // We will take here care ONLY of oscan and nscan flags
  // Adjacent points in i/j direction are consecutive (nscan)
  // Adjacent rows scan in the opposite direction (oscan, only seen in GRIB V2)
  i = (oscan && (j % 2 == 1)) ?  nx - 1 - i : i;
  return (nscan ? j + i*ny : i + nx*j);
}

// Utility function

double earthmap::great_circle(double xlat1, double xlon1,
                              double xlat2, double xlon2)
{
  double a = (90.0-xlat1)*deg2rad;
  double b = (90.0-xlat2)*deg2rad;
  double theta = (xlon2-xlon1)*deg2rad;
  double c = acos((cos(a)*cos(b)) + (sin(a)*sin(b)*cos(theta)));
  return earth_radius_m*c;
}

// Test function main

#ifdef TESTME

using namespace std;

int main(int argc, char *argv[])
{
  int i;
  double ri, rj;
  double lat, lon;
  double ugrid, vgrid;
  double utrue, vtrue;

  earthmap emap;
  t_proj_latlon_parameters llp;
  // ECMWF grid
  //llp.nx   = 201;
  //llp.ny   = 121;
  //llp.lat1 = 60.0;
  //llp.lon1 = -15.0;
  //llp.lat2 = 30.0;
  //llp.lon2 = 35.0;
  //llp.scanflag = 0;
  // GFS Grid
  llp.nx   = 720;
  llp.ny   = 361;
  llp.lat1 = 90;
  llp.lon1 = 0.0;
  llp.lat2 = -90.0;
  llp.lon2 = 359.50;
  llp.scanflag = 0;
  emap.set(llp);
  cout << "Regular ll" << std::endl;
  if (! emap.is_set( )) cerr << "Error setting up the map" << endl;
  for (i = 0; i < 10; i ++)
  {
    ri = rj = i;
    emap.ij_to_latlon(ri, rj, &lat, &lon);
    cout << ri << "," << rj << " : " << lat << "," << lon << endl;
    emap.latlon_to_ij(lat, lon, &ri, &rj);
    cout << ri << "," << rj << " : " << lat << "," << lon << endl;
  }

  emap.truewind_to_gridwind(lat, lon, 1.0, 1.0, &ugrid, &vgrid);
  cout << ugrid << "    " << vgrid << endl;
  emap.gridwind_to_truewind(ri, rj, ugrid, vgrid, &utrue, &vtrue);
  cout << utrue << "    " << vtrue << endl;

  t_proj_mercator_parameters mp;
  mp.nx = 201;
  mp.ny = 121;
  mp.lat1 =  60.0;
  mp.lon1 = -15.0;
  mp.lat2 =  30.0;
  mp.lon2 =  35.0;
  mp.reflat = 12.0;
  mp.dx   = 18000.0;
  mp.dy   = 18000.0;
  mp.scanflag = 0;
  emap.set(mp);
  if (! emap.is_set( )) cerr << "Error setting up the map" << endl;
  cout << "Mercator" << std::endl;
  for (i = 0; i < 10; i ++)
  {
    ri = rj = i;
    emap.ij_to_latlon(ri, rj, &lat, &lon);
    cout << ri << "," << rj << " : " << lat << "," << lon << endl;
    emap.latlon_to_ij(lat, lon, &ri, &rj);
    cout << ri << "," << rj << " : " << lat << "," << lon << endl;
  }

  emap.truewind_to_gridwind(lat, lon, 1.0, 1.0, &ugrid, &vgrid);
  cout << ugrid << "    " << vgrid << endl;
  emap.gridwind_to_truewind(ri, rj, ugrid, vgrid, &utrue, &vtrue);
  cout << utrue << "    " << vtrue << endl;

  t_proj_lambert_parameters lcp;
  lcp.lat1 =  34.654;
  lcp.lon1 =  5.139;
  lcp.resflag = 136;
  lcp.projflag = 0;
  lcp.scanflag = 64;
  lcp.dx   = 18000.0;
  lcp.dy   = 18000.0;
  lcp.orient   = 12.0;
  lcp.truelat1   = 38.0;
  lcp.truelat2   = 45.0;
  lcp.nx = 71;
  lcp.ny = 89;
  emap.set(lcp);
  if (! emap.is_set( )) cerr << "Error setting up the map" << endl;
  cout << "Lambert" << std::endl;
  for (i = 0; i < 10; i ++)
  {
    ri = rj = i;
    emap.ij_to_latlon(ri, rj, &lat, &lon);
    cout << ri << "," << rj << " : " << lat << "," << lon << endl;
    emap.latlon_to_ij(lat, lon, &ri, &rj);
    cout << ri << "," << rj << " : " << lat << "," << lon << endl;
  }

  emap.truewind_to_gridwind(lat, lon, 1.0, 1.0, &ugrid, &vgrid);
  cout << ugrid << "    " << vgrid << endl;
  emap.gridwind_to_truewind(ri, rj, ugrid, vgrid, &utrue, &vtrue);
  cout << utrue << "    " << vtrue << endl;

  t_proj_polar_parameters pp;
  pp.nx = 135;
  pp.ny = 95;
  pp.lat1 = 27.203000;
  pp.lon1 = -135.213000;
  pp.orient = 249.000000;
  pp.dx = 60000.0;
  pp.dy = 60000.0;
  pp.resflag = 136;
  pp.projflag = 0;
  pp.scanflag = 64;
  emap.set(pp);
  if (! emap.is_set( )) cerr << "Error setting up the map" << endl;
  cout << "Polar" << std::endl;
  for (i = 0; i < 10; i ++)
  {
    ri = rj = (double) i;
    emap.ij_to_latlon(ri, rj, &lat, &lon);
    cout << ri << "," << rj << " : " << lat << "," << lon << endl;
    emap.latlon_to_ij(lat, lon, &ri, &rj);
    cout << ri << "," << rj << " : " << lat << "," << lon << endl;
  }

  emap.truewind_to_gridwind(lat, lon, 1.0, 1.0, &ugrid, &vgrid);
  cout << ugrid << "    " << vgrid << endl;
  emap.gridwind_to_truewind(ri, rj, ugrid, vgrid, &utrue, &vtrue);
  cout << utrue << "    " << vtrue << endl;

  t_proj_rotwrf_parameters wp;
  wp.nx = 127;
  wp.ny = 231;
  wp.lat1 = 33.142;
  wp.lon1 = -2.328;
  wp.resflag = 136;
  wp.lat0 = 42.6;
  wp.lon0 = 9.2;
  wp.scanflag = 64;
  wp.resflag = 136;
  emap.set(wp);
  if (! emap.is_set( )) cerr << "Error setting up the map" << endl;
  cout << "Arakawa E on rotated" << std::endl;
  for (i = 0; i < 10; i ++)
  {
    ri = rj = i;
    emap.ij_to_latlon(ri, rj, &lat, &lon);
    cout << ri << "," << rj << " : " << lat << "," << lon << endl;
    emap.latlon_to_ij(lat, lon, &ri, &rj);
    cout << ri << "," << rj << " : " << lat << "," << lon << endl;
  }

  emap.truewind_to_gridwind(lat, lon, 1.0, 1.0, &ugrid, &vgrid);
  cout << ugrid << "    " << vgrid << endl;
  emap.gridwind_to_truewind(ri, rj, ugrid, vgrid, &utrue, &vtrue);
  cout << utrue << "    " << vtrue << endl;

  t_proj_rotlat_parameters rlp;
  rlp.nx = 297;
  rlp.ny = 313;
  rlp.lat1 = -25.0;
  rlp.lon1 = -8.5;
  rlp.lat2 = -6.5;
  rlp.lon2 = 11.0;
  rlp.polelat = -32.5;
  rlp.polelon = 10.0;
  rlp.rotation = 0.0;
  emap.set(rlp);
  if (! emap.is_set( )) cerr << "Error setting up the map" << endl;
  cout << "Rotated LL" << std::endl;
  for (i = 0; i < 10; i ++)
  {
    ri = rj = i;
    emap.ij_to_latlon(ri, rj, &lat, &lon);
    cout << ri << "," << rj << " : " << lat << "," << lon << endl;
    emap.latlon_to_ij(lat, lon, &ri, &rj);
    cout << ri << "," << rj << " : " << lat << "," << lon << endl;
  }

  emap.truewind_to_gridwind(lat, lon, 1.0, 1.0, &ugrid, &vgrid);
  cout << ugrid << "    " << vgrid << endl;
  emap.gridwind_to_truewind(ri, rj, ugrid, vgrid, &utrue, &vtrue);
  cout << utrue << "    " << vtrue << endl;

  return 0;
}

#endif
