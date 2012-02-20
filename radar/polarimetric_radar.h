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

#ifndef __CETEMPS__POLARIMETRIC__RADAR__
#define __CETEMPS__POLARIMETRIC__RADAR__

#include <hydrometeor.h>
#include <hydropsd.h>
#include <string>

namespace cetemps {

  typedef enum {
    polar_distr_unknown = -1,
    polar_distr_normal = 0,
    polar_distr_equipr = 1,
    polar_distr_random = 2
  } t_enum_distr;

  class polardistr {
    public:
      polardistr( );
      void set_gaussian(int np, double d1, double d2, double m, double s);
      void set_uniform(int np, double d1, double d2);
      void set_random(int np, double d1, double d2);
      int Points( ) { return np+1; }
      double MeanValue( ) { return mean; }
      double StandardDeviation( ) { return stdd; }
      double Norm( ) { return norm; }
      double GetVal(int index);
      double GetWeight(int index);
    private:
      double gaussian_integral(double d1, double d2, double m, double s);
      t_enum_distr type;
      int np;
      double d1;
      double d2;
      double range;
      double deltad;
      double mean;
      double stdd;
      double norm;
  };

  class polarimetric_radar {
    public:
      polarimetric_radar(double freqGHz);
      char *documentation( );
      double freqGHz;
      double wavelen_micr;
      double wavelen_mm;
    private:
      float latitude;          // radar latitude in degrees WGS84
      float longitude;         // radar longitude in degrees WGS84
      float elevation;         // radar elevation in msl
      std::string location;    // common cartographic name of the location
      std::string description; // simple mnemonica of radar
      std::string type;        // type of radar (builder, model, serial, ...)
  };

  typedef enum {
    radar_calculate_unknown_method = -1,
    radar_calculate_tmatrix_method = 0
  } t_enum_calculation_method;

  class radar_observation {
    public:
      radar_observation(polarimetric_radar &p);
      void setup(int ndiam, polardistr ad, polardistr bd);
      void tmatrix_calculate(hydrometeor &h, hydroPSD &p,
                             double theta, double phi);
      double ath;
      double zhh;
      double zvv;
      double zhv;
      double zvh;
      double ldrhv;
      double ldrvh;
      double rhohv;
      double dhv;
      double kdp;
      double athh;
      double atvv;
      double athv;
      double atvh;
      double phihv;
      double rrate;
      double wc;
    private:
      double wl;
      double fr;
      int ndiam;
      double ath_coeff;
      polardistr alpha_distr;
      polardistr beta_distr;
  };

}

#endif
