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

#ifndef __CETEMPS__HYDROMETEOR__
#define __CETEMPS__HYDROMETEOR__

#include <particle.h>
#include <complex>
#include <hydropsd.h>

namespace cetemps {

  // ###############################################
  //
  // ATTENTION !!!
  //
  // All wavelenghts are expressed in millimiters
  // All diameters here are expressed in millimiters
  // All temperatures are expressed in Kelvin
  // All pressures are expressed in Pascal
  // All mixing ratios are expressed in kg/kg
  // -----------------------------------------------
  // ###############################################

  typedef enum {
    pruppacher_beard_1970 = 0,
    keenan_carey_zrnic_may_2001 = 1,
    andsager_beard_laird_1999 = 2,
    beard_chuang_1987 = 3,
    brandes_zhang_vivekanandan_2002 = 4
  } t_enum_axisratio_formula;

  class hydrometeor {
    public:
      virtual t_particle_shape GetShape(double diameter) = 0;
      virtual double ShapeParameter(double diameter, hydroPSD &p) = 0;
      virtual int GetPolynomialOrder(double diameter) = 0;
      virtual double Velocity(double diameter) = 0;
      virtual double PrecipRate(double diameter, double delta, double ndd) = 0;
      virtual double Water(double diameter, double delta, double ndd) = 0;
      virtual int BaseGaussianPoints(double diameter) = 0;
      virtual std::complex <double> RefractionIndex(double wavelenght) = 0;
  }; 

  class hydro_rain : public hydrometeor {
    public:
      hydro_rain(t_enum_axisratio_formula ifrar, double temperature,
                 double density);
      virtual t_particle_shape GetShape(double diameter);
      virtual double ShapeParameter(double diameter, hydroPSD &p);
      virtual int GetPolynomialOrder(double diameter) { return 0; }
      virtual double Velocity(double diameter);
      virtual double PrecipRate(double diameter, double delta, double ndd);
      virtual double Water(double diameter, double delta, double ndd);
      virtual int BaseGaussianPoints(double diameter);
      virtual std::complex <double> RefractionIndex(double wavelenght);
    private:
      t_enum_axisratio_formula ifrar;
      double temperature;
      double density;
  };

  class hydro_hail : public hydrometeor {
    public:
      hydro_hail(double pressure, double mixing_ratio, double temperature,
                 double density);
      virtual t_particle_shape GetShape(double diameter);
      virtual double ShapeParameter(double diameter, hydroPSD &p);
      virtual int GetPolynomialOrder(double diameter) { return 0; }
      virtual double Velocity(double diameter);
      virtual double PrecipRate(double diameter, double delta, double ndd);
      virtual double Water(double diameter, double delta, double ndd);
      virtual int BaseGaussianPoints(double diameter);
      virtual std::complex <double> RefractionIndex(double wavelenght);
    private:
      bool is_dry;
      double density;
      double fincl;
      int iceopt;
      double p;
      double q;
      double t;
  };

  class hydro_graupel : public hydrometeor {
    public:
      hydro_graupel(double pressure, double mixing_ratio, double temperature,
                    double density);
      virtual t_particle_shape GetShape(double diameter);
      virtual double ShapeParameter(double diameter, hydroPSD &p);
      virtual int GetPolynomialOrder(double diameter) { return 0; }
      virtual double Velocity(double diameter);
      virtual double PrecipRate(double diameter, double delta, double ndd);
      virtual double Water(double diameter, double delta, double ndd);
      virtual int BaseGaussianPoints(double diameter);
      virtual std::complex <double> RefractionIndex(double wavelenght);
    private:
      bool is_dry;
      double density;
      double fincl;
      int iceopt;
      double p;
      double q;
      double t;
  };

  class hydro_snow : public hydrometeor {
    public:
      hydro_snow(double pressure, double mixing_ratio, double temperature,
                 double density);
      virtual t_particle_shape GetShape(double diameter);
      virtual double ShapeParameter(double diameter, hydroPSD &p);
      virtual int GetPolynomialOrder(double diameter) { return 0; }
      virtual double Velocity(double diameter);
      virtual double PrecipRate(double diameter, double delta, double ndd);
      virtual double Water(double diameter, double delta, double ndd);
      virtual int BaseGaussianPoints(double diameter);
      virtual std::complex <double> RefractionIndex(double wavelenght);
    private:
      bool is_dry;
      double density; 
      double fincl;
      int iceopt;
      double p;
      double q;
      double t;
  };

  class hydro_ice : public hydrometeor {
    public:
      hydro_ice(double pressure, double mixing_ratio, double temperature,
                double density);
      virtual t_particle_shape GetShape(double diameter);
      virtual double ShapeParameter(double diameter, hydroPSD &p);
      virtual int GetPolynomialOrder(double diameter) { return 0; }
      virtual double Velocity(double diameter);
      virtual double PrecipRate(double diameter, double delta, double ndd);
      virtual double Water(double diameter, double delta, double ndd);
      virtual int BaseGaussianPoints(double diameter);
      virtual std::complex <double> RefractionIndex(double wavelenght);
    private:
      double density; 
      double fincl;
      int iceopt;
      double p;
      double q;
      double t;
  };

  class hydro_cloud : public hydrometeor {
    public:
      hydro_cloud(t_enum_axisratio_formula ifrar, double temperature,
                  double density);
      virtual t_particle_shape GetShape(double diameter)
      { return particle_spheroid; }
      virtual double ShapeParameter(double diameter, hydroPSD &p)
      { return 1.0001; }
      virtual int GetPolynomialOrder(double diameter) { return 0; }
      virtual double Velocity(double diameter) { return 0.0; }
      virtual double PrecipRate(double diameter, double delta, double ndd)
      {
        return 0.0;
      }
      virtual double Water(double diameter, double delta, double ndd);
      virtual int BaseGaussianPoints(double diameter) { return 2; }
      virtual std::complex <double> RefractionIndex(double wavelenght);
    private:
      t_enum_axisratio_formula ifrar;
      double temperature;
      double density;
  };

  typedef enum {
    class_large_drops                = 0,
    class_light_rain                 = 1,
    class_medium_rain                = 2,
    class_heavy_rain                 = 3,
    class_hail_rain                  = 4,
    class_hail                       = 5,
    class_graupel_small_hail         = 6,
    class_dry_snow                   = 7,
    class_wet_snow                   = 8,
    class_ice_crystal                = 9,
    class_drizzle_rain               = 10,
    class_wet_hail                   = 11,
    class_wet_graupel_wet_small_hail = 12,
    class_wet_hail_rain              = 13,
    class_cloud_liquid               = 14
  } t_enum_hess_hydro_class;

  // class hydro_volcanic_hash_ice : public hydrometeor {
  // };
  // class hydro_volcanic_hash_rain : public hydrometeor {
  // };

}

#endif
