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


#define DSFMT_MEXP 19937
#include <dSFMT.h>

#include <iostream>
#include <fstream>
#include <random.h>
#include <cstdlib>

using namespace cetemps;

double cetemps::ranf( )
{
  static bool has_seed = false;

  if (! has_seed)
  {
    int seed;
    try {
      std::ifstream urandom;
      urandom.open("/dev/urandom", std::ios::in);
      urandom.read((char *) &seed, 4);
      urandom.close();
    }
    catch (...)
    {
      dsfmt_gv_init_gen_rand(rand());
    }
    dsfmt_gv_init_gen_rand(seed);
    has_seed = true;
  }

  return dsfmt_gv_genrand_open_open( );
}

#ifdef TESTME

#include <iostream>

int main(int argc, char *argv[])
{
  for (int i = 0; i < 100; i ++)
    std::cout << ranf() << std::endl;
  return 0;
}

#endif
