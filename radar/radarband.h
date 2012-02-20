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

#ifndef ___CETEMPS__RADARBANDS___
#define ___CETEMPS__RADARBANDS___

namespace cetemps {

  typedef enum {
    radar_band_hf  = 1,
    radar_band_vhf = 2,
    radar_band_uhf = 3,
    radar_band_l   = 4,
    radar_band_s   = 5,
    radar_band_c   = 6,
    radar_band_x   = 7,
    radat_band_ku  = 8,
    radar_band_k   = 9,
    radar_band_ka  = 10,
    radar_band_mm  = 11
  } t_enum_radar_band;
  
  char bandnames[13][8] = {"LOWER", "HF", "VHF", "UHF", "L",
        "S", "C", "X", "Ku", "K", "Ka", "mm", "HIGHER"};
  float bound_wl_cm[12] = {
        10000, 1000, 100, 30, 15, 8, 4, 2.5, 1.7, 1.2, 0.75, 0.1 };
  float bound_fr_GHz[12] = {
        0.003, 0.03, 0.3, 1, 2, 4, 8, 12, 18, 27, 40, 300 };
  class radarband {
    public:
      char *BandFromWLcm(float wl)
      {
        if (wl > bound_wl_cm[0]) return bandnames[12];
        if (wl < bound_wl_cm[11]) return bandnames[0];
        int ind = 0;
        while (wl < bound_wl_cm[ind+1]) ind++;
        return bandnames[ind+1];
      }
      char *BandFromFrGHz(float fr)
      {
        if (fr < bound_fr_GHz[0]) return bandnames[0];
        if (fr > bound_fr_GHz[11]) return bandnames[12];
        int ind = 0;
        while (fr > bound_fr_GHz[ind+1]) ind++;
        return bandnames[ind+1];
      }
      double Low_wl_cm(t_enum_radar_band b)
      {
        return bound_wl_cm[b];
      }
      double Hi_wl_cm(t_enum_radar_band b)
      {
        return bound_wl_cm[b-1];
      }
      double Low_fr_GHz(t_enum_radar_band b)
      {
        return bound_fr_GHz[b-1];
      }
      double Hi_fr_GHz(t_enum_radar_band b)
      {
        return bound_fr_GHz[b];
      }
  };

};

#endif
