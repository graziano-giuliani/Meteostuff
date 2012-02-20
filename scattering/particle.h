/***************************************************************************
 *   Copyright (C) 2008-2009 Graziano Giuliani                             *
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

#ifndef __CETEMPS__SCATTERING__PARTICLE__
#define __CETEMPS__SCATTERING__PARTICLE__

#include <iostream>
#include <simplecomplex.h>
#include <boost/numeric/ublas/vector.hpp>

using namespace boost::numeric::ublas;

namespace cetemps {

  static const int NPOINT_DROP_CHEBYSHEV = 11;
  static const double drop_chebyshev[NPOINT_DROP_CHEBYSHEV] =
      { -0.0481, 0.0359, -0.1263, 0.0244, 0.0091, -0.0099,
         0.0015, 0.0025, -0.0016, -0.0002, 0.0010 };

  typedef enum {
    radius_equal_volume_sphere       = 1,  // calculate ratio
    radius_equal_surface_area_sphere = 2   // given ratio
  } t_sphere_radius_type;
  
  typedef enum {
    particle_raindrop  = -3,
    particle_cylinder  = -2,
    particle_spheroid  = -1,
    particle_chebyshev =  0
  } t_particle_shape;

  class scattering_particle {
    public:
      t_particle_shape shape;

      double radius; // particle radius in microns
      double ratio;  // effective axis ratio
      double rev;    // radius * ratio

      int nmax(double wavelenght);

      dcplex refractive_index;

      virtual void rsp(vector <double> &x,
                       vector <double> &r, vector <double> &dr) {}

      friend std::ostream& operator<< (std::ostream& os,
          scattering_particle &p);
  };

  class spheroid : public scattering_particle {
    public:
      spheroid(double radius, double axis_ratio, t_sphere_radius_type type,
               dcplex refractive_index);
      void rsp(vector <double> &x, vector <double> &r, vector <double> &dr);
      double get_eps( ) { return eps; }
    private:
      double eps;
  };

  class chebyshev : public scattering_particle {
    public:
      chebyshev(double radius, double deformation,
                int polynomial_degree, t_sphere_radius_type type,
                dcplex refractive_index);
      void rsp(vector <double> &x, vector <double> &r, vector <double> &dr);
      int get_polynomial_degree( ) { return np; };
    private:
      double eps;
      int np;
  };

  class cylinder : public scattering_particle {
    public:
      cylinder(double radius, double dim_ratio, t_sphere_radius_type type,
               dcplex refractive_index);
      void rsp(vector <double> &x, vector <double> &r, vector <double> &dr);
      double get_hv_dimension_ratio( ) { return eps; }
    private:
      double eps;
  };

  class raindrop : public scattering_particle {
    public:
      raindrop(double radius, double axis_ratio,
               dcplex refractive_index);
      void rsp(vector <double> &x, vector <double> &r, vector <double> &dr);
      double get_rov( ) { return rov; }
    private:
      double rov;
  };

}

#endif
