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

#include <cmath>
#include <hydrometeor.h>
#include <random.h>
#include <refindex.h>

using namespace cetemps;

hydro_rain::hydro_rain(t_enum_axisratio_formula ifrar, double temperature,
                       double density)
{
  this->ifrar = ifrar;
  this->temperature = temperature;
  this->density = density;
}

t_particle_shape hydro_rain::GetShape(double diam)
{
  return particle_spheroid;
}

double hydro_rain::ShapeParameter(double diam, hydroPSD &p)
{
  switch(ifrar)
  {
    case pruppacher_beard_1970:
      if (diam >= 0.3) return 1.0/(1.03-0.062*diam);
      break;
    case keenan_carey_zrnic_may_2001:
      return (1.0/(0.9939+
              0.00736*diam - 0.018485*diam*diam + 0.001456*diam*diam*diam));
      break;
    case andsager_beard_laird_1999:
      if(diam > 1.0 && diam <= 4.0)
        return (1.0/(1.012-0.01445*diam-0.01028*diam*diam));
      else
        return (1.0/(1.0048 +
                0.00057*diam - 0.02628*diam*diam +
                0.003682*diam*diam*diam - 0.0001677*diam*diam*diam*diam));
      break;
    case beard_chuang_1987:
      if (diam >= 0.3)
        return (1.0/(1.0048+
                0.00057*diam - 0.02628*diam*diam +
                0.003682*diam*diam*diam - 0.0001677*diam*diam*diam*diam));
      break;
    case brandes_zhang_vivekanandan_2002:
      if (diam >= 0.3)
        return (1.0/(0.9951 +
                0.02510*diam - 0.03644*diam*diam +
                0.005030*diam*diam*diam - 0.0002492*diam*diam*diam*diam));
      break;
    default:
      return 1.00001;
      break;
  }
  return 1.00001;
}

double hydro_rain::Velocity(double diam)
{
  return (3.78*pow(diam, 0.67));
}

double hydro_rain::PrecipRate(double diam, double delta, double ndd)
{
  double vel = Velocity(diam);
  return (0.6*M_PI*density*1e-3*pow(diam, 3)*ndd*vel*delta);
}

double hydro_rain::Water(double diam, double delta, double ndd)
{
  return (M_PI/6.0)*density*1e-3*pow(diam, 3)*ndd*delta;
}

int hydro_rain::BaseGaussianPoints(double diam)
{
  if (diam <= 0.5) return 2;
  else if (diam <= 1.0) return 4;
  else if (diam <= 3.0) return 8;
  else return 16;
}

std::complex <double> hydro_rain::RefractionIndex(double wavelenght)
{
  return refindex(wateps(wavelenght, temperature));
}

hydro_hail::hydro_hail(double pressure, double mixing_ratio,
                       double temperature, double density)
{
  p = pressure;
  q = mixing_ratio;
  t = temperature;

  if (temperature > 273.15) // Wet hail
  {
    is_dry = false;
    fincl = 0.05+0.35*ranf();
    this->density = (fincl+(1-fincl))*density+0.2*ranf();
    this->density = density;
    iceopt = 3;
  }
  else
  {
    is_dry = true;
    fincl = 0.0;
    this->density = density;
    iceopt = 0;
  }
}

t_particle_shape hydro_hail::GetShape(double diam)
{
  return particle_spheroid;
}

double hydro_hail::ShapeParameter(double diam, hydroPSD &p)
{
  if (diam < 10) return 1.0+0.1*ranf();
  return 1.0+0.67*ranf();
}

double hydro_hail::Velocity(double diam)
{
  return pow(((4.0*0.9*diam*9.81)/(3*1.123*0.5)), 0.5);
}

double hydro_hail::PrecipRate(double diam, double delta, double ndd)
{
  double vel = Velocity(diam);
  return (0.6*M_PI*density*1e-3*(diam*diam*diam)*ndd*vel*delta*0.5);
}

double hydro_hail::Water(double diam, double delta, double ndd)
{
  return (density*(M_PI/6.0)*1e-3*(diam*diam*diam)*ndd*delta);
}

int hydro_hail::BaseGaussianPoints(double diam)
{
  return 2;
}

std::complex <double> hydro_hail::RefractionIndex(double wavelenght)
{
  return refindex(iceeps(iceopt, density, fincl, wavelenght, p, q, t));
}

t_particle_shape hydro_graupel::GetShape(double diam)
{
  return particle_spheroid;
}

hydro_graupel::hydro_graupel(double pressure, double mixing_ratio,
                             double temperature, double density)
{
  p = pressure;
  q = mixing_ratio;
  t = temperature;

  if (temperature > 273.15) // Wet graupel
  {
    is_dry = false;
    fincl = 0.20+0.25*ranf();
    this->density = (fincl+(1-fincl))*density;
    this->density = density;
    iceopt = 2;
  }
  else
  {
    is_dry = true;
    fincl = 0.05+0.35*ranf();
    this->density = (fincl+(1-fincl))*density;
    this->density = density;
    iceopt = 2;
  }
}

double hydro_graupel::ShapeParameter(double diam, hydroPSD &p)
{
  double dmin = p.minvalue( );
  double dmax = p.maxvalue( );
  return (0.7*((diam-dmin)/(dmax-dmin))+1.0+0.1*ranf());
}

double hydro_graupel::Velocity(double diam)
{
  return pow(((4.0*0.9*diam*9.81)/(3.0*1.123*0.5)), 0.5);
}

double hydro_graupel::PrecipRate(double diam, double delta, double ndd)
{
  double vel = Velocity(diam);
  return (0.6*M_PI*density*1e-3*diam*diam*diam*ndd*vel*delta*0.5);
}

double hydro_graupel::Water(double diam, double delta, double ndd)
{
  return (density*(M_PI/6.0)*1e-3*(diam*diam*diam)*ndd*delta);
}

int hydro_graupel::BaseGaussianPoints(double diam)
{
  return 4;
}

std::complex <double> hydro_graupel::RefractionIndex(double wavelenght)
{
  return refindex(iceeps(iceopt, density, fincl, wavelenght, p, q, t));
}

hydro_snow::hydro_snow(double pressure, double mixing_ratio,
                       double temperature, double density)
{
  p = pressure;
  q = mixing_ratio;
  t = temperature;

  if (temperature > 273.15) // Wet Snow
  {
    is_dry = false;
    fincl = 0.01+0.09*ranf();
    this->density = (1-fincl)*density;
    this->density = density;
    iceopt = 2;
  }
  else
  {
    is_dry = true;
    fincl = 0.45+0.40*ranf();
    this->density = (fincl+(1-fincl))*density;
    this->density = density;
    iceopt = 1;
  }
}

t_particle_shape hydro_snow::GetShape(double diam)
{
  return particle_spheroid;
}

double hydro_snow::ShapeParameter(double diam, hydroPSD &p)
{
  double dmin = p.minvalue( );
  double dmax = p.maxvalue( );
  if (is_dry)
    return (0.2*((diam-dmin)/(dmax-dmin))+1.0+0.1*ranf());
  else
    return (0.8*((diam-dmin)/(dmax-dmin))+1.0+0.1*ranf());
}

double hydro_snow::Velocity(double diam)
{
  return (4.836*pow((diam*1e-3), 0.25));
}

double hydro_snow::PrecipRate(double diam, double delta, double ndd)
{
  double vel = Velocity(diam);
  return (0.6*M_PI*density*1e-3*diam*diam*diam*ndd*vel*delta*0.1);
}

double hydro_snow::Water(double diam, double delta, double ndd)
{
  return (density*(M_PI/6.0)*1e-3*(diam*diam*diam)*ndd*delta);
}

int hydro_snow::BaseGaussianPoints(double diam)
{
  return 8;
}

std::complex <double> hydro_snow::RefractionIndex(double wavelenght)
{
  return refindex(iceeps(iceopt, density, fincl, wavelenght, p, q, t));
}

hydro_ice::hydro_ice(double pressure, double mixing_ratio, double temperature,
                     double density)
{
  p = pressure;
  q = mixing_ratio;
  t = temperature;
  this->density = density;
  fincl = 0.0;
  iceopt = 0;
}

t_particle_shape hydro_ice::GetShape(double diam)
{
  return particle_spheroid;
}

double hydro_ice::ShapeParameter(double diam, hydroPSD &p)
{
  double dmin = p.minvalue( );
  double dmax = p.maxvalue( );
  return 1.0*((diam-dmin)/(dmax-dmin))+1.0+0.1*ranf();
}

double hydro_ice::Velocity(double diam)
{
  return pow(((4.0*0.9*diam*9.81)/(3.0*1.123*0.5)), 0.5);
}

double hydro_ice::PrecipRate(double diam, double delta, double ndd)
{
  double vel = Velocity(diam);
  return (0.6*M_PI*density*1e-3*diam*diam*diam*ndd*vel*delta*0.5);
}

double hydro_ice::Water(double diam, double delta, double ndd)
{
   return density*(M_PI/6.0)*1e-3*(diam*diam*diam)*ndd*delta;
}

int hydro_ice::BaseGaussianPoints(double diam)
{
  return 16;
}

std::complex <double> hydro_ice::RefractionIndex(double wavelenght)
{
  return refindex(iceeps(iceopt, density, fincl, wavelenght, p, q, t));
}

hydro_cloud::hydro_cloud(t_enum_axisratio_formula ifrar, double temperature,
                         double density)
{
  this->ifrar = ifrar;
  this->temperature = temperature;
  this->density = density;
}

double hydro_cloud::Water(double diam, double delta, double ndd)
{
  return (M_PI/6.0)*density*1e-3*(diam*diam*diam)*ndd*delta;
}

std::complex <double> hydro_cloud::RefractionIndex(double wavelenght)
{
  return refindex(wateps(wavelenght, temperature));
}

#ifdef TESTME

#include <hydropsd.h>

using namespace cetemps;

int main(int argc, char *argv[])
{
  double freq = 9.6; //GHz
  double alam = 299792.458/freq*0.001;

  std::cout << "Using wavelenght of " << alam << " mm ("
            << freq << " GHz)" << std::endl;

  // Large drops
  {
    double dzero = 2.50*ranf()+2.0;
    double mu = 1.81*ranf()-0.94;
    double nw = 100*ranf()+5;
    double temp = -5+40*ranf();
    double rmin = 0.3;
    double rmax = 0.5*(dzero*3.2);
    if (rmax > 3.5) rmax = 3.5;
    GeneralizedGammaPSD psd(nw, mu, dzero);
    psd.SetExtremes(rmin*2, rmax*2);

    hydro_rain rain(pruppacher_beard_1970, temp+273.15, 0.997);

    std::complex <double> refindex = rain.RefractionIndex(alam);
    double kappa2 = KappaSquare(refindex);

    std::cout << "Temperature is " << temp << " deg Celsius" << std::endl;
    std::cout << "RI: " << refindex << std::endl;
    std::cout << "K2: " << kappa2 << std::endl;

    double dmin = 2*rmin;
    double daxi = (2*rmax-2*rmin) / 51;
    for (int i = 0; i < 50; i ++)
    {
      double diam = dmin + (i-1) * daxi;
      double ndd = psd.ndd(diam);
      std::cout << "Diam : " << diam << std::endl;
      std::cout << "Ndd : " << ndd << std::endl;
      std::cout << "Shape : " << rain.GetShape(diam) << std::endl;
      std::cout << "Eps: " << rain.ShapeParameter(diam, psd) << std::endl;
      std::cout << "Vel: " << rain.Velocity(diam) << std::endl;
      std::cout << "Rate: " << rain.PrecipRate(diam, daxi, ndd) << std::endl;
      std::cout << "Water: " << rain.Water(diam, daxi, ndd) << std::endl;
      std::cout << "NDGS: " << rain.BaseGaussianPoints(diam) << std::endl;
    }
  }

  // Light rain
  {
    double dzero = 0.850*ranf()+0.550;
    double mu = 4.999999*ranf()-0.999999;
    double nw = 20000*ranf()+1000;
    double temp = 40*ranf();
    double rmin = 0.3;
    double rmax = 0.5*(dzero*3.2);
    if (rmax > 3.5) rmax = 3.5;
    GeneralizedGammaPSD psd(nw, mu, dzero);
    psd.SetExtremes(rmin*2, rmax*2);

    hydro_rain rain(pruppacher_beard_1970, temp+273.15, 0.997);

    std::complex <double> refindex = rain.RefractionIndex(alam);
    double kappa2 = KappaSquare(refindex);

    std::cout << "Temperature is " << temp << " deg Celsius" << std::endl;
    std::cout << "RI: " << refindex << std::endl;
    std::cout << "K2: " << kappa2 << std::endl;

    double dmin = 2*rmin;
    double daxi = (2*rmax-2*rmin) / 51;
    for (int i = 0; i < 50; i ++)
    {
      double diam = dmin + (i-1) * daxi;
      double ndd = psd.ndd(diam);
      std::cout << "Diam : " << diam << std::endl;
      std::cout << "Ndd : " << ndd << std::endl;
      std::cout << "Shape : " << rain.GetShape(diam) << std::endl;
      std::cout << "Eps: " << rain.ShapeParameter(diam, psd) << std::endl;
      std::cout << "Vel: " << rain.Velocity(diam) << std::endl;
      std::cout << "Rate: " << rain.PrecipRate(diam, daxi, ndd) << std::endl;
      std::cout << "Water: " << rain.Water(diam, daxi, ndd) << std::endl;
      std::cout << "NDGS: " << rain.BaseGaussianPoints(diam) << std::endl;
    }
  }

  // Medium rain
  {
    double dzero = 0.500*ranf()+1.400;
    double mu = 4.999999*ranf()-0.999999;
    double nw = 9000*ranf()+1000;
    double temp = 40*ranf();
    double rmin = 0.3;
    double rmax = 0.5*(dzero*3.2);
    if (rmax > 3.5) rmax = 3.5;
    GeneralizedGammaPSD psd(nw, mu, dzero);
    psd.SetExtremes(rmin*2, rmax*2);

    hydro_rain rain(pruppacher_beard_1970, temp+273.15, 0.997);

    std::complex <double> refindex = rain.RefractionIndex(alam);
    double kappa2 = KappaSquare(refindex);

    std::cout << "Temperature is " << temp << " deg Celsius" << std::endl;
    std::cout << "RI: " << refindex << std::endl;
    std::cout << "K2: " << kappa2 << std::endl;

    double dmin = 2*rmin;
    double daxi = (2*rmax-2*rmin) / 51;
    for (int i = 0; i < 50; i ++)
    {
      double diam = dmin + (i-1) * daxi;
      double ndd = psd.ndd(diam);
      std::cout << "Diam : " << diam << std::endl;
      std::cout << "Ndd : " << ndd << std::endl;
      std::cout << "Shape : " << rain.GetShape(diam) << std::endl;
      std::cout << "Eps: " << rain.ShapeParameter(diam, psd) << std::endl;
      std::cout << "Vel: " << rain.Velocity(diam) << std::endl;
      std::cout << "Rate: " << rain.PrecipRate(diam, daxi, ndd) << std::endl;
      std::cout << "Water: " << rain.Water(diam, daxi, ndd) << std::endl;
      std::cout << "NDGS: " << rain.BaseGaussianPoints(diam) << std::endl;
    }
  }

  // Heavy rain
  {
    double dzero = 1.350*ranf()+1.850;
    double mu = 4.999999*ranf()-0.999999;
    double nw = 7000*ranf()+2000;
    double temp = 40*ranf();
    double rmin = 0.3;
    double rmax = 0.5*(dzero*3.2);
    if (rmax > 3.5) rmax = 3.5;
    GeneralizedGammaPSD psd(nw, mu, dzero);
    psd.SetExtremes(rmin*2, rmax*2);

    hydro_rain rain(pruppacher_beard_1970, temp+273.15, 0.997);

    std::complex <double> refindex = rain.RefractionIndex(alam);
    double kappa2 = KappaSquare(refindex);

    std::cout << "Temperature is " << temp << " deg Celsius" << std::endl;
    std::cout << "RI: " << refindex << std::endl;
    std::cout << "K2: " << kappa2 << std::endl;

    double dmin = 2*rmin;
    double daxi = (2*rmax-2*rmin) / 51;
    for (int i = 0; i < 50; i ++)
    {
      double diam = dmin + (i-1) * daxi;
      double ndd = psd.ndd(diam);
      std::cout << "Diam : " << diam << std::endl;
      std::cout << "Ndd : " << ndd << std::endl;
      std::cout << "Shape : " << rain.GetShape(diam) << std::endl;
      std::cout << "Eps: " << rain.ShapeParameter(diam, psd) << std::endl;
      std::cout << "Vel: " << rain.Velocity(diam) << std::endl;
      std::cout << "Rate: " << rain.PrecipRate(diam, daxi, ndd) << std::endl;
      std::cout << "Water: " << rain.Water(diam, daxi, ndd) << std::endl;
      std::cout << "NDGS: " << rain.BaseGaussianPoints(diam) << std::endl;
    }
  }

  // Hail/rain
  {
    double dzero = 1.50*ranf()+1.90;
    double mu = 4.999999*ranf()-0.999999;
    double nw = 7000*ranf()+2000;
    double temp = -10+40*ranf();
    double rmin = 0.5;
    double rmax = 4.0;
    GeneralizedGammaPSD psd(nw, mu, dzero);
    psd.SetExtremes(rmin*2, rmax*2);

    temp = -3.0;
    hydro_hail hail(85000, 0.001, temp+273.15, 0.600);

    std::complex <double> refindex = hail.RefractionIndex(alam);
    double kappa2 = KappaSquare(refindex);

    std::cout << "Temperature is " << temp << " deg Celsius" << std::endl;
    std::cout << "RI: " << refindex << std::endl;
    std::cout << "K2: " << kappa2 << std::endl;

    double dmin = 2*rmin;
    double daxi = (2*rmax-2*rmin) / 51;
    for (int i = 0; i < 50; i ++)
    {
      double diam = dmin + (i-1) * daxi;
      double ndd = psd.ndd(diam);
      std::cout << "Diam : " << diam << std::endl;
      std::cout << "Ndd : " << ndd << std::endl;
      std::cout << "Shape : " << hail.GetShape(diam) << std::endl;
      std::cout << "Eps: " << hail.ShapeParameter(diam, psd) << std::endl;
      std::cout << "Vel: " << hail.Velocity(diam) << std::endl;
      std::cout << "Rate: " << hail.PrecipRate(diam, daxi, ndd) << std::endl;
      std::cout << "Water: " << hail.Water(diam, daxi, ndd) << std::endl;
      std::cout << "NDGS: " << hail.BaseGaussianPoints(diam) << std::endl;
    }
  }

  // Hail
  {
    double m0 = 100+200*ranf();
    double slope = 0.4+0.6*ranf();
    double temp = -20+40*ranf();
    double rmin = 2.5;
    double rmax = 15.0;
    ExponentialPSD psd(m0, slope);
    psd.SetExtremes(rmin*2, rmax*2);

    hydro_hail hail(85000, 0.001, temp+273.15, 0.500);

    std::complex <double> refindex = hail.RefractionIndex(alam);
    double kappa2 = KappaSquare(refindex);

    std::cout << "Temperature is " << temp << " deg Celsius" << std::endl;
    std::cout << "RI: " << refindex << std::endl;
    std::cout << "K2: " << kappa2 << std::endl;

    double dmin = 2*rmin;
    double daxi = (2*rmax-2*rmin) / 51;
    for (int i = 0; i < 50; i ++)
    {
      double diam = dmin + (i-1) * daxi;
      double ndd = psd.ndd(diam);
      std::cout << "Diam : " << diam << std::endl;
      std::cout << "Ndd : " << ndd << std::endl;
      std::cout << "Shape : " << hail.GetShape(diam) << std::endl;
      std::cout << "Eps: " << hail.ShapeParameter(diam, psd) << std::endl;
      std::cout << "Vel: " << hail.Velocity(diam) << std::endl;
      std::cout << "Rate: " << hail.PrecipRate(diam, daxi, ndd) << std::endl;
      std::cout << "Water: " << hail.Water(diam, daxi, ndd) << std::endl;
      std::cout << "NDGS: " << hail.BaseGaussianPoints(diam) << std::endl;
    }
  }

  // Graupel/small hail
  {
    double m0 = 10+260*ranf();
    double slope = 0.5+1.1*ranf();
    double temp = -50+60*ranf();
    double rmin = 0.5;
    double rmax = 2.5;
    ExponentialPSD psd(m0, slope);
    psd.SetExtremes(rmin*2, rmax*2);

    hydro_graupel graup(85000, 0.001, temp+273.15, 0.500);

    std::complex <double> refindex = graup.RefractionIndex(alam);
    double kappa2 = KappaSquare(refindex);

    std::cout << "Temperature is " << temp << " deg Celsius" << std::endl;
    std::cout << "RI: " << refindex << std::endl;
    std::cout << "K2: " << kappa2 << std::endl;

    double dmin = 2*rmin;
    double daxi = (2*rmax-2*rmin) / 51;
    for (int i = 0; i < 50; i ++)
    {
      double diam = dmin + (i-1) * daxi;
      double ndd = psd.ndd(diam);
      std::cout << "Diam : " << diam << std::endl;
      std::cout << "Ndd : " << ndd << std::endl;
      std::cout << "Shape : " << graup.GetShape(diam) << std::endl;
      std::cout << "Eps: " << graup.ShapeParameter(diam, psd) << std::endl;
      std::cout << "Vel: " << graup.Velocity(diam) << std::endl;
      std::cout << "Rate: " << graup.PrecipRate(diam, daxi, ndd) << std::endl;
      std::cout << "Water: " << graup.Water(diam, daxi, ndd) << std::endl;
      std::cout << "NDGS: " << graup.BaseGaussianPoints(diam) << std::endl;
    }
  }

  // Dry snow
  {
    double rate = 0.05+1.0*ranf();
    double m0 = 2.5*1e+03/pow(rate, 0.94);
    double slope = 2.29/pow(rate,0.45);
    double temp = -51+51*ranf();
    double rmin = 0.5;
    double rmax = 7.5;
    ExponentialPSD psd(m0, slope);
    psd.SetExtremes(rmin*2, rmax*2);

    hydro_snow snow(85000, 0.001, temp+273.15, 0.100);

    std::complex <double> refindex = snow.RefractionIndex(alam);
    double kappa2 = KappaSquare(refindex);

    std::cout << "Temperature is " << temp << " deg Celsius" << std::endl;
    std::cout << "RI: " << refindex << std::endl;
    std::cout << "K2: " << kappa2 << std::endl;

    double dmin = 2*rmin;
    double daxi = (2*rmax-2*rmin) / 51;
    for (int i = 0; i < 50; i ++)
    {
      double diam = dmin + (i-1) * daxi;
      double ndd = psd.ndd(diam);
      std::cout << "Diam : " << diam << std::endl;
      std::cout << "Ndd : " << ndd << std::endl;
      std::cout << "Shape : " << snow.GetShape(diam) << std::endl;
      std::cout << "Eps: " << snow.ShapeParameter(diam, psd) << std::endl;
      std::cout << "Vel: " << snow.Velocity(diam) << std::endl;
      std::cout << "Rate: " << snow.PrecipRate(diam, daxi, ndd) << std::endl;
      std::cout << "Water: " << snow.Water(diam, daxi, ndd) << std::endl;
      std::cout << "NDGS: " << snow.BaseGaussianPoints(diam) << std::endl;
    }
  }

  // Wet snow
  {
    double rate = 0.5+1.2*ranf();
    double m0 = 2.5*1e+03/pow(rate, 0.94);
    double slope = 2.29/pow(rate,0.45);
    double temp = -2.5+5*ranf();
    double rmin = 0.5;
    double rmax = 7.5;
    ExponentialPSD psd(m0, slope);
    psd.SetExtremes(rmin*2, rmax*2);

    hydro_snow snow(85000, 0.001, temp+273.15, 0.200);

    std::complex <double> refindex = snow.RefractionIndex(alam);
    double kappa2 = KappaSquare(refindex);

    std::cout << "Temperature is " << temp << " deg Celsius" << std::endl;
    std::cout << "RI: " << refindex << std::endl;
    std::cout << "K2: " << kappa2 << std::endl;

    double dmin = 2*rmin;
    double daxi = (2*rmax-2*rmin) / 51;
    for (int i = 0; i < 50; i ++)
    {
      double diam = dmin + (i-1) * daxi;
      double ndd = psd.ndd(diam);
      std::cout << "Diam : " << diam << std::endl;
      std::cout << "Ndd : " << ndd << std::endl;
      std::cout << "Shape : " << snow.GetShape(diam) << std::endl;
      std::cout << "Eps: " << snow.ShapeParameter(diam, psd) << std::endl;
      std::cout << "Vel: " << snow.Velocity(diam) << std::endl;
      std::cout << "Rate: " << snow.PrecipRate(diam, daxi, ndd) << std::endl;
      std::cout << "Water: " << snow.Water(diam, daxi, ndd) << std::endl;
      std::cout << "NDGS: " << snow.BaseGaussianPoints(diam) << std::endl;
    }
  }

  // Ice crystals
  {
    double m0 = 1+40*ranf();
    double slope = 1.1+1.9*ranf();
    double temp = -72.5+67*ranf();
    double rmin = 0.1;
    double rmax = 2.4;
    ExponentialPSD psd(m0, slope);
    psd.SetExtremes(rmin*2, rmax*2);

    hydro_ice ice(85000, 0.001, temp+273.15, 0.400);

    std::complex <double> refindex = ice.RefractionIndex(alam);
    double kappa2 = KappaSquare(refindex);

    std::cout << "Temperature is " << temp << " deg Celsius" << std::endl;
    std::cout << "RI: " << refindex << std::endl;
    std::cout << "K2: " << kappa2 << std::endl;

    double dmin = 2*rmin;
    double daxi = (2*rmax-2*rmin) / 51;
    for (int i = 0; i < 50; i ++)
    {
      double diam = dmin + (i-1) * daxi;
      double ndd = psd.ndd(diam);
      std::cout << "Diam : " << diam << std::endl;
      std::cout << "Ndd : " << ndd << std::endl;
      std::cout << "Shape : " << ice.GetShape(diam) << std::endl;
      std::cout << "Eps: " << ice.ShapeParameter(diam, psd) << std::endl;
      std::cout << "Vel: " << ice.Velocity(diam) << std::endl;
      std::cout << "Rate: " << ice.PrecipRate(diam, daxi, ndd) << std::endl;
      std::cout << "Water: " << ice.Water(diam, daxi, ndd) << std::endl;
      std::cout << "NDGS: " << ice.BaseGaussianPoints(diam) << std::endl;
    }
  }

  // Drizzle
  {
    double dzero = 0.130*ranf()+0.360;
    double mu = 4.999999*ranf()-0.999999;
    double nw = 7000*ranf()+14000;
    double temp = 40*ranf();
    double rmin = 0.1;
    double rmax = 1.5;
    GeneralizedGammaPSD psd(nw, mu, dzero);
    psd.SetExtremes(rmin*2, rmax*2);

    hydro_rain rain(pruppacher_beard_1970, temp+273.15, 0.997);

    std::complex <double> refindex = rain.RefractionIndex(alam);
    double kappa2 = KappaSquare(refindex);

    std::cout << "Temperature is " << temp << " deg Celsius" << std::endl;
    std::cout << "RI: " << refindex << std::endl;
    std::cout << "K2: " << kappa2 << std::endl;

    double dmin = 2*rmin;
    double daxi = (2*rmax-2*rmin) / 51;
    for (int i = 0; i < 50; i ++)
    {
      double diam = dmin + (i-1) * daxi;
      double ndd = psd.ndd(diam);
      std::cout << "Diam : " << diam << std::endl;
      std::cout << "Ndd : " << ndd << std::endl;
      std::cout << "Shape : " << rain.GetShape(diam) << std::endl;
      std::cout << "Eps: " << rain.ShapeParameter(diam, psd) << std::endl;
      std::cout << "Vel: " << rain.Velocity(diam) << std::endl;
      std::cout << "Rate: " << rain.PrecipRate(diam, daxi, ndd) << std::endl;
      std::cout << "Water: " << rain.Water(diam, daxi, ndd) << std::endl;
      std::cout << "NDGS: " << rain.BaseGaussianPoints(diam) << std::endl;
    }
  }

  // Wet hail
  {
    double m0 = 100+200*ranf();
    double slope = 0.4+0.6*ranf();
    double temp = -5+25*ranf();
    double rmin = 2.5;
    double rmax = 15.0;
    ExponentialPSD psd(m0, slope);
    psd.SetExtremes(rmin*2, rmax*2);

    hydro_hail hail(85000, 0.001, temp+273.15, 0.500);

    std::complex <double> refindex = hail.RefractionIndex(alam);
    double kappa2 = KappaSquare(refindex);

    std::cout << "Temperature is " << temp << " deg Celsius" << std::endl;
    std::cout << "RI: " << refindex << std::endl;
    std::cout << "K2: " << kappa2 << std::endl;

    double dmin = 2*rmin;
    double daxi = (2*rmax-2*rmin) / 51;
    for (int i = 0; i < 50; i ++)
    {
      double diam = dmin + (i-1) * daxi;
      double ndd = psd.ndd(diam);
      std::cout << "Diam : " << diam << std::endl;
      std::cout << "Ndd : " << ndd << std::endl;
      std::cout << "Shape : " << hail.GetShape(diam) << std::endl;
      std::cout << "Eps: " << hail.ShapeParameter(diam, psd) << std::endl;
      std::cout << "Vel: " << hail.Velocity(diam) << std::endl;
      std::cout << "Rate: " << hail.PrecipRate(diam, daxi, ndd) << std::endl;
      std::cout << "Water: " << hail.Water(diam, daxi, ndd) << std::endl;
      std::cout << "NDGS: " << hail.BaseGaussianPoints(diam) << std::endl;
    }
  }

  // Wet graupel/wet small hail
  {
    double m0 = 10+260*ranf();
    double slope = 0.5+1.1*ranf();
    double temp = -5+25*ranf();
    double rmin = 0.5;
    double rmax = 2.5;
    ExponentialPSD psd(m0, slope);
    psd.SetExtremes(rmin*2, rmax*2);

    hydro_graupel graupel(85000, 0.001, temp+273.15, 0.500);

    std::complex <double> refindex = graupel.RefractionIndex(alam);
    double kappa2 = KappaSquare(refindex);

    std::cout << "Temperature is " << temp << " deg Celsius" << std::endl;
    std::cout << "RI: " << refindex << std::endl;
    std::cout << "K2: " << kappa2 << std::endl;

    double dmin = 2*rmin;
    double daxi = (2*rmax-2*rmin) / 51;
    for (int i = 0; i < 50; i ++)
    {
      double diam = dmin + (i-1) * daxi;
      double ndd = psd.ndd(diam);
      std::cout << "Diam : " << diam << std::endl;
      std::cout << "Ndd : " << ndd << std::endl;
      std::cout << "Shape : " << graupel.GetShape(diam) << std::endl;
      std::cout << "Eps: " << graupel.ShapeParameter(diam, psd) << std::endl;
      std::cout << "Vel: " << graupel.Velocity(diam) << std::endl;
      std::cout << "Rate: " << graupel.PrecipRate(diam, daxi, ndd) << std::endl;
      std::cout << "Water: " << graupel.Water(diam, daxi, ndd) << std::endl;
      std::cout << "NDGS: " << graupel.BaseGaussianPoints(diam) << std::endl;
    }
  }

  // Wet hail/rain
  {
    double dzero = 1.50*ranf()+1.90;
    double mu = 4.999999*ranf()-0.999999;
    double nw = 7000*ranf()+2000;
    double temp = -5+35*ranf();
    double rmin = 0.5;
    double rmax = 4.0;
    GeneralizedGammaPSD psd(nw, mu, dzero);
    psd.SetExtremes(rmin*2, rmax*2);

    hydro_rain rain(pruppacher_beard_1970, temp+273.15, 0.900);

    std::complex <double> refindex = rain.RefractionIndex(alam);
    double kappa2 = KappaSquare(refindex);

    std::cout << "Temperature is " << temp << " deg Celsius" << std::endl;
    std::cout << "RI: " << refindex << std::endl;
    std::cout << "K2: " << kappa2 << std::endl;

    double dmin = 2*rmin;
    double daxi = (2*rmax-2*rmin) / 51;
    for (int i = 0; i < 50; i ++)
    {
      double diam = dmin + (i-1) * daxi;
      double ndd = psd.ndd(diam);
      std::cout << "Diam : " << diam << std::endl;
      std::cout << "Ndd : " << ndd << std::endl;
      std::cout << "Shape : " << rain.GetShape(diam) << std::endl;
      std::cout << "Eps: " << rain.ShapeParameter(diam, psd) << std::endl;
      std::cout << "Vel: " << rain.Velocity(diam) << std::endl;
      std::cout << "Rate: " << rain.PrecipRate(diam, daxi, ndd) << std::endl;
      std::cout << "Water: " << rain.Water(diam, daxi, ndd) << std::endl;
      std::cout << "NDGS: " << rain.BaseGaussianPoints(diam) << std::endl;
    }
  }

  // Liquid Cloud
  {
    double dzero = 0.200*ranf()+0.100;
    double mu = 4.999999*ranf()+3.0;
    double nw = 7000*ranf()+20000;
    double temp = 40*ranf();
    double rmin = 0.01;
    double rmax = 0.5;
    GeneralizedGammaPSD psd(nw, mu, dzero);
    psd.SetExtremes(rmin*2, rmax*2);

    hydro_rain rain(pruppacher_beard_1970, temp+273.15, 0.997);

    std::complex <double> refindex = rain.RefractionIndex(alam);
    double kappa2 = KappaSquare(refindex);

    std::cout << "Temperature is " << temp << " deg Celsius" << std::endl;
    std::cout << "RI: " << refindex << std::endl;
    std::cout << "K2: " << kappa2 << std::endl;

    double dmin = 2*rmin;
    double daxi = (2*rmax-2*rmin) / 51;
    for (int i = 0; i < 50; i ++)
    {
      double diam = dmin + (i-1) * daxi;
      double ndd = psd.ndd(diam);
      std::cout << "Diam : " << diam << std::endl;
      std::cout << "Ndd : " << ndd << std::endl;
      std::cout << "Shape : " << rain.GetShape(diam) << std::endl;
      std::cout << "Eps: " << rain.ShapeParameter(diam, psd) << std::endl;
      std::cout << "Vel: " << rain.Velocity(diam) << std::endl;
      std::cout << "Rate: " << rain.PrecipRate(diam, daxi, ndd) << std::endl;
      std::cout << "Water: " << rain.Water(diam, daxi, ndd) << std::endl;
      std::cout << "NDGS: " << rain.BaseGaussianPoints(diam) << std::endl;
    }
  }

  return 0;
}

#endif
