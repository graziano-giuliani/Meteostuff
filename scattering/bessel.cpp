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

#ifdef TESTME

#include <bessel.h>

#include <iostream>
#include <boost/numeric/ublas/io.hpp>
#include <simplecomplex.h>

using namespace boost::numeric::ublas;
using namespace cetemps;

int main(int argc, char *argv[])
{
  int order = 10;
  double x = 5.50;
  vector <double> y(10);
  vector <double> u(10);

  int iprec = 10;

  rjb(iprec, x, y, u);
  std::cout << y << std::endl;
  std::cout << u << std::endl;
  ryb(x, y, u);
  std::cout << y << std::endl;
  std::cout << u << std::endl;
  std::complex<double> cx(5.0,2.0);
  vector <dcplex> cy(10);
  vector <dcplex> cu(10);
  rjb(iprec, cx, cy, cu);
  std::cout << cy << std::endl;
  std::cout << cu << std::endl;
  return(0);
}

#endif
