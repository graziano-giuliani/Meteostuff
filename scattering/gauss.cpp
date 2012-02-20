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
#include <iostream>
#include <gauss.h>

namespace cetemps {

  using namespace boost::numeric::ublas;

//
//
//  Purpose: 
//
//    LEGENDRE_COM computes abscissas and weights for
//    Gauss-Legendre quadrature.
//
//  Integration interval:
//
//    [ -1, 1 ]
//
//  Weight function:
//
//    1.
//
//  Integral to approximate:
//
//    Integral ( -1 <= X <= 1 ) F(X) dX.
//
//  Approximate integral:
//
//    sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NORDER, the order of the rule.
//    NORDER must be greater than 0.
//
//    Output, double XTAB[NORDER], the abscissas of the rule.
//
//    Output, double WEIGHT[NORDER], the weights of the rule.
//    The weights are positive, symmetric, and should sum to 2.
//
  void gauss(int norder, vector <double> &xtab, vector <double> &weight)
  {
    BOOST_UBLAS_CHECK(norder <= xtab.size(), external_logic());
    BOOST_UBLAS_CHECK(norder <= weight.size(), external_logic());

    double d1;
    double d2pn;
    double d3pn;
    double d4pn;
    double dp;
    double dpn;
    double e1;
    double fx;
    double h;
    int i;
    int iback;
    int k;
    int m;
    int mp1mi;
    int ncopy;
    int nmove;
    double p;
    double pk;
    double pkm1;
    double pkp1;
    double t;
    double u;
    double v;
    double x0;
    double xtemp;

    if ( norder < 1 )
    {
      std::cout << std::endl
                << "GAUSS - Fatal error!" << std::endl
                << "  Illegal value of NORDER = " << norder << std::endl;
      throw 1;
    }
   
    e1 = ( double ) ( norder * ( norder + 1 ) );
    m = ( norder + 1 ) / 2;
   
    for ( i = 1; i <= ( norder + 1 ) / 2; i++ )
    {
      mp1mi = m + 1 - i;
      t = M_PI * ( double ) ( 4 * i - 1 ) / ( double ) ( 4 * norder + 2 );
      x0 = cos(t) * ( 1.0 - ( 1.0 - 1.0 / 
        ( double ) ( norder ) ) / ( double ) ( 8 * norder * norder ) );
      pkm1 = 1.0;
      pk = x0;

      for ( k = 2; k <= norder; k++ )
      {
        pkp1 = 2.0 * x0 * pk - pkm1 - ( x0 * pk - pkm1 ) / ( double ) ( k );
        pkm1 = pk;
        pk = pkp1;
      }
   
      d1 = ( double ) ( norder ) * ( pkm1 - x0 * pk );
      dpn = d1 / ( 1.0 - x0 * x0 );
      d2pn = ( 2.0 * x0 * dpn - e1 * pk ) / ( 1.0 - x0 * x0 );
      d3pn = ( 4.0 * x0 * d2pn + ( 2.0 - e1 ) * dpn ) / ( 1.0 - x0 * x0 );
      d4pn = ( 6.0 * x0 * d3pn + ( 6.0 - e1 ) * d2pn ) / ( 1.0 - x0 * x0 );

      u = pk / dpn;
      v = d2pn / dpn;
  //
  //  Initial approximation H:
  //
      h = - u * ( 1.0 + 0.5 * u * ( v + u * ( v * v - d3pn 
        / ( 3.0 * dpn ) ) ) );
  //
  //  Refine H using one step of Newton's method:
  //
      p = pk + h * ( dpn + 0.5 * h * ( d2pn + h / 3.0 
        * ( d3pn + 0.25 * h * d4pn ) ) );

      dp = dpn + h * ( d2pn + 0.5 * h * ( d3pn + h * d4pn / 3.0 ) );
      h = h - p / dp;
      xtemp = x0 + h;
      xtab(mp1mi-1) = xtemp;
      fx = d1 - h * e1 * ( pk + 0.5 * h * ( dpn + h / 3.0 
        * ( d2pn + 0.25 * h * ( d3pn + 0.2 * h * d4pn ) ) ) );
      weight(mp1mi-1) = 2.0 * ( 1.0 - xtemp * xtemp ) / ( fx * fx ); 
    }
   
    if ( ( norder % 2 ) == 1 )
    {
      xtab(0) = 0.0;
    }
  //
  //  Shift the data up.
  //
    nmove = ( norder + 1 ) / 2;
    ncopy = norder - nmove;

    for ( i = 1; i <= nmove; i++ )
    {
      iback = norder + 1 - i;
      xtab(iback-1) = xtab(iback-ncopy-1);
      weight(iback-1) = weight(iback-ncopy-1);
    }
  //
  //  Reflect values for the negative abscissas.
  //
    for ( i = 0; i < norder - nmove; i++ )
    {
      xtab(i) = - xtab(norder-1-i);
      weight(i) = weight(norder-1-i);
    }
    return;
  }
}

#ifdef TESTME

#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric::ublas;

int main(int argc, char *argv[])
{
  int order = 10;
  vector <double> x(10);
  vector <double> w(10);

  cetemps::gauss(order, x, w);

  std::cout << x << std::endl;
  std::cout << w << std::endl;
  return(0);
}

#endif
