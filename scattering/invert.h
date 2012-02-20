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

#ifndef __CETEMPS__MATRIX__INVERT__
#define __CETEMPS__MATRIX__INVERT__

#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/bindings/lapack/lapack.hpp>

namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;

template <typename T, typename F, typename A>
BOOST_UBLAS_INLINE
int lapack_invert(const ublas::matrix<T,F,A> &a, ublas::matrix<T,F,A> &b)
{
  BOOST_UBLAS_CHECK(a.size1() == a.size2(), ublas::external_logic());
  BOOST_UBLAS_CHECK(a.size1() == b.size1(), ublas::external_logic());
  BOOST_UBLAS_CHECK(a.size2() == b.size2(), ublas::external_logic());

  ublas::permutation_matrix<int> piv(a.size1());
  ublas::matrix<T, ublas::column_major, A> ca(a);

  int info = lapack::getrf(ca, piv);
  if (info == 0)
  {
    info = lapack::getri(ca, piv);
    if (info == 0)
      b = ca;
  }
  return info;
}

template <typename T, typename F, typename A>
BOOST_UBLAS_INLINE
int lapack_invert_inplace(ublas::matrix<T,F,A> &a)
{
  BOOST_UBLAS_CHECK(a.size1() == a.size2(), ublas::external_logic());

  ublas::permutation_matrix<int> piv(a.size1());
  ublas::matrix<T, ublas::column_major, A> ca(a);

  int info = lapack::getrf(ca, piv);
  if (info == 0)
  {
    info = lapack::getri(ca, piv);
    if (info == 0)
      a = ca;
  }
  return info;
}

#endif
