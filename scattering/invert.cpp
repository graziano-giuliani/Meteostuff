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

#include <simplecomplex.h>
#include <invert.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace ublas = boost::numeric::ublas;

int main(int argc, char *argv[])
{
  std::size_t n = 3;
  ublas::matrix<dcplex> m(n, n);
  ublas::matrix<dcplex> x(n, n);

  m(0, 0) = 3; m(0, 1) = 2; m(0, 2) = 1;
  m(1, 0) = 1; m(1, 1) = 3; m(1, 2) = 1;
  m(2, 0) = 2; m(2, 1) = 2; m(2, 2) = 1;

  int info = lapack_invert(m, x);
  if (info != 0) return -1;

  std::cout << m << std::endl;
  std::cout << x << std::endl;
  std::cout << prod(m, x) << std::endl;

  return 0;
};

#endif
