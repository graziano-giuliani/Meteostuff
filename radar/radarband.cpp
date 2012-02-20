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

#ifdef TESTME

#include <radarband.h>
#include <iostream>

using namespace cetemps;

int main(int argc, char*argv[])
{
  radarband rb;

  std::cout << "wl of 14 cm is " << rb.BandFromWLcm(14) << std::endl;
  std::cout << "wl of 0.01 cm is " << rb.BandFromWLcm(0.01) << std::endl;
  std::cout << "wl of 100000 cm is " << rb.BandFromWLcm(100000) << std::endl;
  std::cout << "wl of 10000 cm is " << rb.BandFromWLcm(10000) << std::endl;
  std::cout << "wl of 0.1 cm is " << rb.BandFromWLcm(0.1) << std::endl;

  std::cout << "fr of 9.6 GHz is " << rb.BandFromFrGHz(9.6) << std::endl;
  std::cout << "fr of 0.001 GHz is " << rb.BandFromFrGHz(0.001) << std::endl;
  std::cout << "fr of 350 GHz is " << rb.BandFromFrGHz(350) << std::endl;
  std::cout << "fr of 0.003 GHz is " << rb.BandFromFrGHz(0.003) << std::endl;
  std::cout << "fr of 300 GHz is " << rb.BandFromFrGHz(300) << std::endl;

  std::cout << "Extremes of X band in GHz : (" << rb.Low_fr_GHz(radar_band_x)
    << ", " << rb.Hi_fr_GHz(radar_band_x) << ")" << std::endl;
  std::cout << "Extremes of C band in cm : (" << rb.Low_wl_cm(radar_band_c)
    << ", " << rb.Hi_wl_cm(radar_band_c) << ")" << std::endl;
  return 0;
}

#endif
