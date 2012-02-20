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

#include <cmath>
#include <particle.h>
#include <gauss.h>
#include <iostream>

namespace cetemps {

  int scattering_particle::nmax(double wavelenght)
  {
    const double one_third =  1.0/3.0;
    double xev = 2.0 * M_PI * (rev / wavelenght);
    int imax = (int) (xev + 4.05 * pow(xev, one_third));
    imax = imax < 4 ? 4 : imax;
    return imax;
  }

  spheroid::spheroid(double radius,
                     double axis_ratio,
                     t_sphere_radius_type type,
                     dcplex refractive_index)
  {
    shape = particle_spheroid;
    this->radius = radius;
    this->refractive_index = refractive_index;
    eps = axis_ratio;

    const double one_third =  1.0/3.0;
    const double two_third =  2.0/3.0;
    const double four_third = 4.0/3.0;
    switch(type)
    {
      case radius_equal_surface_area_sphere:
        if (axis_ratio < 1.0)  // Prolate spheroid
        {
          double e = sqrt(1.0 - (axis_ratio*axis_ratio));
          double r = sqrt(0.50 * (pow(axis_ratio, two_third) +
            pow(axis_ratio, -one_third) * (asin(e)/e)));
          ratio = 1.0/r;
        }
        else                   // Oblate spheroid
        {
          double e = sqrt(1.0 - (1.0 / (axis_ratio*axis_ratio)));
          double r = sqrt(0.25 * (2.0 * pow(axis_ratio, two_third) +
            pow(axis_ratio, -four_third) * log((1.0+e)/(1.0-e))/e));
          ratio = 1.0/r;
        }
        rev = radius*ratio;
        break;
      case radius_equal_volume_sphere:
        ratio = axis_ratio; // accept input value
        rev = radius;
        break;
      default:
        throw "Particle size not correctly specified";
    }
  }

  void spheroid::rsp(vector <double> &x,
                     vector <double> &r, vector <double> &dr)
  {
    int ng = r.size( );
    int ngauss = ng / 2;
    BOOST_UBLAS_CHECK(ng == dr.size( ), external_logic());
    BOOST_UBLAS_CHECK(x.size() >= ngauss, external_logic());

    double a, c, s, aa, cc, ee, rr, ss, ee1;

    a = rev * pow(eps, 1.0/3.0);
    aa = a * a;
    ee = eps * eps;
    ee1 = ee - 1.0;
    for (int i = 0; i < ngauss; ++i)
    {
	    c = x(i);
	    cc = c * c;
	    ss = 1.0 - cc;
	    s = sqrt(ss);
	    rr = 1.0 / (ss + ee * cc);
	    r(i) = aa * rr;
	    r(ng - i - 1) = r(i);
	    dr(i) = rr * c * s * ee1;
	    dr(ng - i - 1) = -dr(i);
    }
    return;
  }

  chebyshev::chebyshev(double radius,
                       double deformation, 
                       int polynomial_degree,
                       t_sphere_radius_type type,
                       dcplex refractive_index)
  {
    const double one_third =  1.0/3.0;
    const double three_fourth = 3.0/4.0;
    static const int ng = 60;

    shape = particle_chebyshev;
    this->radius = radius;
    this->refractive_index = refractive_index;
    eps = deformation;
    np = polynomial_degree;
    switch(type)
    {
      case radius_equal_surface_area_sphere:
        {
          vector <double> x(ng);
          vector <double> w(ng);
          double a, s, v, a2, dn, en;
          double ds, dx, rs, rv, dcn, dsn, ens, dxn;

          dn = (double) (polynomial_degree);
          en = deformation * dn;
          gauss(ng, x, w);
          s = 0.;
          v = 0.;
          for (int i = 0; i < ng; i++)
          {
          	dx = acos(x(i));
    	      dxn = dn * dx;
        	  ds = sin(dx);
          	dsn = sin(dxn);
          	dcn = cos(dxn);
          	a = deformation * dcn + 1.;
          	a2 = a * a;
          	ens = en * dsn;
          	s += w(i) * a * sqrt(a2 + ens * ens);
          	v += w(i) * (ds * a + x(i) * ens) * ds * a2;
          }
          rs = sqrt(s * 0.5);
          rv = pow((v * three_fourth), one_third);
          ratio = rv / rs;
        }
        rev = radius*ratio;
        break;
      case radius_equal_volume_sphere:
        ratio = deformation; // accept input value
        rev = radius;
        break;
      default:
        throw "Particle size not correctly specified";
        break;
    }
  }

  void chebyshev::rsp(vector <double> &x,
                      vector <double> &r, vector <double> &dr)
  {
    int ng = r.size( );
    BOOST_UBLAS_CHECK(ng == dr.size( ), external_logic());
    BOOST_UBLAS_CHECK(ng == x.size( ), external_logic());

    double a, r0, dn, ep, ri, xi, dn4, dnp;

    dnp = (double) (np);
    dn = dnp * dnp;
    dn4 = dn * 4.0;
    ep = eps * eps;
    a = 1.0 + ep * 1.5 * ((dn4 - 2.0) / (dn4 - 1.0));
    int i = (int) ((dnp + 0.1) * 0.5);
    i *= 2;
    if (i == np) {
	    a = a - eps * 3.0 * (ep * 0.25 + 1.0) / (dn - 1.0) -
          ep * 0.25 * eps / (dn * 9.0 - 1.0);
    }
    r0 = rev * pow(a, (-1.0/3.0));
    for (i = 0; i < ng; i++)
    {
    	xi = acos(x(i)) * dnp;
	    ri = r0 * (eps * cos(xi) + 1.0);
	    r(i) = ri * ri;
	    dr(i) = -r0 * eps * dnp * (sin(xi) / ri);
    }
    return;
  }

  cylinder::cylinder(double radius,
                     double dim_ratio,
                     t_sphere_radius_type type,
                     dcplex refractive_index)
  {
    const double one_third =  1.0/3.0;

    shape = particle_cylinder;
    this->radius = radius;
    this->refractive_index = refractive_index;
    eps = dim_ratio;
    switch(type)
    {
      case radius_equal_surface_area_sphere:
        ratio = pow((1.5 / dim_ratio), one_third);
        ratio = ratio / sqrt((dim_ratio+2.0)/(2.0*dim_ratio));
        rev = radius*ratio;
        break;
      case radius_equal_volume_sphere:
        ratio = dim_ratio; // accept input value
        rev = radius;
        break;
      default:
        throw "Particle size not correctly specified";
        break;
    }
  }

  void cylinder::rsp(vector <double> &x,
                     vector <double> &r, vector <double> &dr)
  {
    int ng = r.size( );
    int ngauss = ng / 2;
    BOOST_UBLAS_CHECK(ng == dr.size( ), external_logic());
    BOOST_UBLAS_CHECK(x.size() >= ngauss, external_logic());

    double a, h;
    double co, si, rad, rthet;

    h = rev * pow((2.0 / (eps * 3.0 * eps)), (1.0/3.0));
    a = h * eps;
    for (int i = 0; i < ngauss; i++)
    {
	    co = -x(i);
    	si = sqrt(1.0 - co * co);
	    if (si / co > a / h)
      {
    	  rad = a / si;
	      rthet = -a * co / (si * si);
	    }
      else
      {
	      rad = h / co;
	      rthet = h * si / (co * co);
      }
	    r(i) = rad * rad;
	    r(ng - i - 1) = r(i);
	    dr(i) = -rthet / rad;
	    dr(ng - i - 1) = -dr(i);
    }
    return;
  }

  raindrop::raindrop(double radius, double dim_ratio,
                     dcplex refractive_index)
  {
    const double one_third =  1.0/3.0;
    static const int ng = 60;

    shape = particle_raindrop;
    this->radius = radius;
    this->refractive_index = refractive_index;
    this->ratio  = dim_ratio;
    vector <double> x(ng);
    vector <double> w(ng);
    double s, v, ci, ri, si, wi;
    double xi, rs, rv, dri, xin, risi;
    gauss(ng, x, w);
    s = 0.;
    v = 0.;
    for (int i = 0; i < ng; i++)
    {
    	 xi = acos(x(i));
       wi = w(i);
       ri = drop_chebyshev[0] + 1.;
    	 dri = 0.0;
    	 for (int n = 1; n < NPOINT_DROP_CHEBYSHEV; ++n)
       {
    	   xin = xi * n;
    	   ri += drop_chebyshev[n] * cos(xin);
    	   dri -= drop_chebyshev[n] * n * sin(xin);
    	 }
    	 si = sin(xi);
    	 ci = x(i);
       risi = ri * si;
       s += wi * ri * sqrt(ri * ri + dri * dri);
    	 v += wi * ri * risi * (risi - dri * ci);
    }
    rs = sqrt(s * 0.5);
    rv = pow((v * 3.0 * 0.25), one_third);
    if ((fabs(ratio - 1.0)) > 1e-8) {
    	ratio = rv / rs;
    }
    rov = 1.0 / rv;
    rev = radius*ratio;
  }

  void raindrop::rsp(vector <double> &x,
                     vector <double> &r, vector <double> &dr)
  {
    int ng = r.size( );
    BOOST_UBLAS_CHECK(ng == dr.size( ), external_logic());
    BOOST_UBLAS_CHECK(ng == x.size( ), external_logic());

    double r0, ri, xi, dri, xin;

    r0 = rev * rov;
    for (int i = 0; i < ng; i++)
    {
	    xi = acos(x(i));
	    ri = drop_chebyshev[0] + 1.0;
	    dri = 0.0;
	    for (int n = 1; n < NPOINT_DROP_CHEBYSHEV; n++)
      {
	      xin = xi * n;
	      ri += drop_chebyshev[n] * cos(xin);
	      dri -= drop_chebyshev[n] * n * sin(xin);
      }
	    ri *= r0;
	    dri *= r0;
	    r(i) = ri * ri;
	    dr(i) = dri / ri;
	  }
    return;
  }

  std::ostream& operator<< (std::ostream& os, scattering_particle &p)
  {
    os << "Particle modeled as ";

    switch (p.shape)
    {
      case particle_spheroid:
        {
          spheroid *s = (spheroid *) &p;
          if (s->get_eps( ) > 0.0)
            os << "Oblate ";
          else
            os << "Prolate ";
          os << "Spheroid.";
        }
        break;
      case particle_cylinder:
        {
          cylinder *c = (cylinder *) &p;
          os << "Cylinder with dimension ratio a/b "
             << c->get_hv_dimension_ratio( );
        }
        break;
      case particle_chebyshev:
        {
          chebyshev *y = (chebyshev *) &p;
          os << "Generic Chebyshev with polynomial degree " <<
            y->get_polynomial_degree( ) << ".";
        }
        break;
      case particle_raindrop:
        os << "Raindrop model (Chebyshev 11 deg).";
        break;
      default:
        os << "???????????";
        break;
    }
    os << std::endl;
    os << "Radius              : " << p.radius << " in units of wavelenght"
       << std::endl
       << "Axis ratio          : " << p.ratio << std::endl
       << "Radius/ratio        : " << p.rev << std::endl
       << "Refractive index is : " << p.refractive_index;
    return os;
  }

}

#ifdef TESTME

#include <iostream>
#include <boost/numeric/ublas/io.hpp>

using namespace cetemps;
using namespace boost::numeric;
using namespace std;

int main(int argc, char *argv[])
{
  spheroid a(10.0, 0.5, radius_equal_volume_sphere, (1.5,0.2));

  cout << a << endl;
  cout << a.nmax(31.2500000000000) << endl;

  ublas::vector <double> x(144);
  ublas::vector <double> r(144);
  ublas::vector <double> dr(144);

  x[0] = -0.99945;
  x[1] = -0.997103;
  x[2] = -0.992885;
  x[3] = -0.986803;
  x[4] = -0.978869;
  x[5] = -0.969097;
  x[6] = -0.957505;
  x[7] = -0.944116;
  x[8] = -0.928955;
  x[9] = -0.912049;
  x[10] = -0.893431;
  x[11] = -0.873136;
  x[12] = -0.851202;
  x[13] = -0.827669;
  x[14] = -0.802583;
  x[15] = -0.77599;
  x[16] = -0.74794;
  x[17] = -0.718486;
  x[18] = -0.687683;
  x[19] = -0.655589;
  x[20] = -0.622264;
  x[21] = -0.587771;
  x[22] = -0.552175;
  x[23] = -0.515542;
  x[24] = -0.477941;
  x[25] = -0.439442;
  x[26] = -0.400119;
  x[27] = -0.360044;
  x[28] = -0.319294;
  x[29] = -0.277944;
  x[30] = -0.236072;
  x[31] = -0.193757;
  x[32] = -0.151079;
  x[33] = -0.108116;
  x[34] = -0.0649512;
  x[35] = -0.0216639;
  x[36] = 0.0216639;
  x[37] = 0.0649512;
  x[38] = 0.108116;
  x[39] = 0.151079;
  x[40] = 0.193757;
  x[41] = 0.236072;
  x[42] = 0.277944;
  x[43] = 0.319294;
  x[44] = 0.360044;
  x[45] = 0.400119;
  x[46] = 0.439442;
  x[47] = 0.477941;
  x[48] = 0.515542;
  x[49] = 0.552175;
  x[50] = 0.587771;
  x[51] = 0.622264;
  x[52] = 0.655589;
  x[53] = 0.687683;
  x[54] = 0.718486;
  x[55] = 0.74794;
  x[56] = 0.77599;
  x[57] = 0.802583;
  x[58] = 0.827669;
  x[59] = 0.851202;
  x[60] = 0.873136;
  x[61] = 0.893431;
  x[62] = 0.912049;
  x[63] = 0.928955;
  x[64] = 0.944116;
  x[65] = 0.957505;
  x[66] = 0.969097;
  x[67] = 0.978869;
  x[68] = 0.986803;
  x[69] = 0.992885;
  x[70] = 0.997103;
  x[71] = 0.99945;
  a.rsp(x, r, dr);
  cout << r << endl;
  cout << dr << endl;

  chebyshev b(10.0, 0.5, 4, radius_equal_volume_sphere, (1.5,0.2));
  b.rsp(x, r, dr);
  cout << r << endl;
  cout << dr << endl;

  cylinder c(10.0, 0.5, radius_equal_volume_sphere, (1.5,0.2));
  c.rsp(x, r, dr);
  cout << r << endl;
  cout << dr << endl;

  raindrop d(10.0, 0.5, (1.5,0.2));
  d.rsp(x, r, dr);
  cout << r << endl;
  cout << dr << endl;

  return 0;
}

#endif
