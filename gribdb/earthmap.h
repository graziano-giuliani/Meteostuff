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

#ifndef __EARTHMAP_H__
#define __EARTHMAP_H__

#include <cmath>

namespace himet {

  typedef enum {
    PROJ_G1_LATLON = 0,
    PROJ_G1_MERC   = 1,
    PROJ_G1_LC     = 3,
    PROJ_G1_PS     = 5,
    PROJ_G1_ROTLAT = 10,
    PROJ_G1_ROTWRF = 203
  } t_enum_g1_projection;

  typedef enum {
    PROJ_G2_LATLON = 0,
    PROJ_G2_ROTLAT = 1,
    PROJ_G2_MERC   = 10,
    PROJ_G2_PS     = 20,
    PROJ_G2_LC     = 30,
    PROJ_G2_ROTWRF = 32768
  } t_enum_g2_projection;

  typedef struct {
    long nx;
    long ny;
    double lat1;
    double lon1;
    double lat2;
    double lon2;
    int scanflag;
  } t_proj_latlon_parameters;

  typedef struct {
    long nx;
    long ny;
    double lat1;
    double lon1;
    double orient;
    double dx;
    double dy;
    double truelat1;
    double truelat2;
    int resflag;
    int projflag;
    int scanflag;
  } t_proj_lambert_parameters;

  typedef struct {
    long nx;
    long ny;
    double lat1;
    double lon1;
    double orient;
    double dx;
    double dy;
    int resflag;
    int projflag;
    int scanflag;
  } t_proj_polar_parameters;

  typedef struct {
    long nx;
    long ny;
    double lat1;
    double lon1;
    double lat2;
    double lon2;
    double dx;
    double dy;
    double reflat;
    int scanflag;
  } t_proj_mercator_parameters;

  typedef struct {
    long nx;
    long ny;
    double lat1;
    double lon1;
    double lat2;
    double lon2;
    double polelat;
    double polelon;
    double rotation;
    int scanflag;
  } t_proj_rotlat_parameters;

  typedef struct {
    long nx;
    long ny;
    double lat1;
    double lon1;
    double lat0;
    double lon0;
    int resflag;
    int scanflag;
  } t_proj_rotwrf_parameters;

  class earthmap {

    public:

      earthmap( );

      long get_size( ) { return (this->nx*this->ny); }
      long get_nx( ) { return (this->nx); }
      long get_ny( ) { return (this->ny); }

      bool is_set( ) { return this->init; }

      void set( t_proj_latlon_parameters &parms );
      void set( t_proj_lambert_parameters &parms );
      void set( t_proj_polar_parameters &parms );
      void set( t_proj_mercator_parameters &parms );
      void set( t_proj_rotlat_parameters &parms );
      void set( t_proj_rotwrf_parameters &parms );

      void latlon_to_ij( double lat, double lon, double *i, double *j );
      void ij_to_latlon( double i, double j, double *lat, double *lon );

      void gridwind_to_truewind(double i, double j,
                                double ugrid, double vgrid,
                                double *utrue, double *vtrue);
      void truewind_to_gridwind(double lat, double lon,
                                double utrue, double vtrue,
                                double *ugrid, double *vgrid);

      double great_circle(double lat1, double lon1, double lat2, double lon2);

      int ij2p(int i, int j);
      void ij2pl(double i, double j, double *pixel, double *line);

      void rotuv_ll2gr(double lat, double lon, double *crot, double *srot);
      void rotuv_gr2ll(double i, double j, double *crot, double *srot);

    private:

      // parameters

      double lat1;     // Latitude (1,1) in degrees (-90->90N)
      double lon1;     // Longitude (1,1) in degrees (-180->180E)
      double lat0;     // Center latitude (NX/2, NY/2) in degrees (-90->90N)
      double lon0;     // Center longitude (NX/2, NY/2) in degrees (-180->180E)
      double lat2;     // Latitude (NX, NY) in degrees (-90->90N)
      double lon2;     // Longitude (NX, NY) in degrees (-180->180E)
      int iscan;       // I direction point scan order
      int jscan;       // J direction point scan order
      int kscan;       // Reserved in GDS Octet 28, used NCEP 203 as stagger
      int nscan;       // Indexing along lines or columns
      int oscan;       // rows scanning in opposite direction
      int irot;        // Wind rotated on grid or easterly and northerly
      int iproj;       // Projection flag for N/S pole and tang/symmetric
      double dlat;     // Lat increment for lat/lon grids
      double dlon;     // Lon increment for lat/lon grids
      double dx;       // Grid spacing in meters at truelats (ps, lc, merc proj)
      double dy;       // Grid spacing in meters at truelats (ps, lc, merc proj)
      double orient;   // The orientation of the grid; i.e., the east longitude
                       // value of the meridian which is parallel to the y-axis
                       // (or columns of the grid) along which latitude
                       // increases as the y-coordinate increases.
                       // (Note: The orientation longitude may, or may not,
                       // appear within a particular grid.)
      double truelat1; // First true latitude (all projections)
      double truelat2; // Second true lat (LC only)
      double rotation; // Clockwise rotation of rotated globe in degrees

      bool  init;     // Flag to indicate if this class is ready

      long nx;
      long ny;

      double xmin, xmax, ymin, ymax;

      double hi, hj; // Scanning flag on i,j plane: default is WE:NS (scan=0)

      double pp, dxs, dys, xp, yp, cone,
             de, de2, invcone2;                   // LC-PS internal params
      double ye;                                  // Mercator internal params
      double zsinpol, zcospol, lat00, lon00;      // ROTLL internal parms
      double dlats, dlons, clat0, slat0, is1;     // ROTWRF internal parms

      void (himet::earthmap::*llij)(double, double, double*, double*);
      void (himet::earthmap::*ijll)(double, double, double*, double*);
      void (himet::earthmap::*gluv)(double, double, double*, double*);
      void (himet::earthmap::*lguv)(double, double, double*, double*);

      void llij_latlon(double lat, double lon, double *i, double *j);
      void ijll_latlon(double i, double j, double *lat, double *lon);
      void llij_lc(double lat, double lon, double *i, double *j);
      void ijll_lc(double i, double j, double *lat, double *lon);
      void llij_ps(double lat, double lon, double *i, double *j);
      void ijll_ps(double i, double j, double *lat, double *lon);
      void llij_merc(double lat, double lon, double *i, double *j);
      void ijll_merc(double i, double j, double *lat, double *lon);
      void llij_rl(double lat, double lon, double *i, double *j);
      void ijll_rl(double i, double j, double *lat, double *lon);
      void llij_rwrf(double lat, double lon, double *i, double *j);
      void ijll_rwrf(double i, double j, double *lat, double *lon);

      void grid2lluv_latlon(double i, double j, double *crot, double *srot);
      void ll2griduv_latlon(double lat, double lon, double *crot, double *srot);
      void grid2lluv_lc(double i, double j, double *crot, double *srot);
      void ll2griduv_lc(double lat, double lon, double *crot, double *srot);
      void grid2lluv_ps(double i, double j, double *crot, double *srot);
      void ll2griduv_ps(double lat, double lon, double *crot, double *srot);
      void grid2lluv_merc(double i, double j, double *crot, double *srot);
      void ll2griduv_merc(double lat, double lon, double *crot, double *srot);
      void grid2lluv_rl(double i, double j, double *crot, double *srot);
      void ll2griduv_rl(double lat, double lon, double *crot, double *srot);
      void grid2lluv_rwrf(double i, double j, double *crot, double *srot);
      void ll2griduv_rwrf(double lat, double lon, double *crot, double *srot);
  };

} // namespace himet

#endif
