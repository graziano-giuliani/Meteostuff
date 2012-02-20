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

#ifndef __CETEMPS__BESSEL__FUNCTIONS__
#define __CETEMPS__BESSEL__FUNCTIONS__

#include <cmath>
#include <boost/numeric/ublas/vector.hpp>

namespace cetemps {

  using namespace boost::numeric::ublas;

//
// Calculation of spherical Bessel functions of the first/second kind
// j=jr+i*ji of real/complex argument x=xr+i*xi of orders from 1 to nmax
// by using backward recursion. Parametr nnmax determines numerical
// accuracy. u=ur+i*ui - function (1/x)(d/dx)(x*j(x))
//

  template <typename T> 
  void rjb(int iprec, T x, vector <T> &y, vector <T> &u)
  {
    int nmax = y.size();
    // BOOST_UBLAS_CHECK(nmax == u.size( ), external_logic());
    int l = nmax + iprec;

    vector <T> z(l);
    int i, i1;
    T y0, z0, y1, yi, xx, yi1;

    xx = 1.0 / x;
    z(l-1) = 1.0 / ((double) ((2*l) + 1) * xx);
    for (i = 1; i < l; i++)
    {
    	i1 = l-i;
	    z(i1-1) = 1.0 / ((double) ((2*i1) + 1) * xx - z(i1));
    }
    z0 = 1.0 / (xx - z(0));
    y0 = z0 * cos(x) * xx;
    y1 = y0 * z(0);
    u(0) = y0 - y1 * xx;
    y(0) = y1;
    for (i = 1; i < nmax; i++)
    {
	    yi1 = y(i-1);
	    yi = yi1 * z(i);
	    y(i) = yi;
	    u(i) = yi1 - (double) (i+1) * yi * xx;
    }
    return;
  }

  template <typename T> 
  void ryb(T x, vector <T> &y, vector <T> &v)
  {
    int nmax = y.size();
    // BOOST_UBLAS_CHECK(nmax == v.size( ), external_logic());

    T c, s, x1, x2, x3, y1;

    c = cos(x);
    s = sin(x);
    x1 = 1.0 / x;
    x2 = x1 * x1;
    x3 = x2 * x1;
    y1 = -c * x2 - s * x1;
    y(0) = y1;
    y(1) = (x3 * -3.0 + x1) * c - x2 * 3.0 * s;
    for (int i = 2; i < nmax; i++)
	    y(i) = (double) ((i * 2) + 1) * x1 * y(i-1) - y(i-2);
    v(0) = -x1 * (c + y1);
    for (int i = 1; i < nmax; i++)
	    v(i) = y(i - 1) - (double) (i+1) * x1 * y(i);
    return;
  }
}

#endif //__CETEMPS__BESSEL__FUNCTIONS__
