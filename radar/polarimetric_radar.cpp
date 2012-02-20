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

#include <polarimetric_radar.h>
#include <refindex.h>
#include <tmatrix.h>
#include <cmath>
#include <random.h>

using namespace cetemps;

polardistr::polardistr()
{
  np = 0;
  d1 = 0;
  d2 = 0;
  mean = 0;
  stdd = 0;
  norm = 0;
  type = polar_distr_unknown;
}

double polardistr::gaussian_integral(double d1, double d2, double m, double s)
{
  double f1 = 0.5*(1+erf((d1-m)/(s*M_SQRT2)));
  double f2 = 0.5*(1+erf((d2-m)/(s*M_SQRT2)));
  return (f2-f1);
}

void polardistr::set_gaussian(int np, double d1, double d2, double m, double s)
{
  this->np = np;
  this->d1 = d1;
  this->d2 = d2;
  mean = m;
  stdd = s;
  range = (d2-d1);
  deltad = range/np;
  norm = 0.0;
  for (int n = 0; n <= np; n ++)
  {
    double d = (d1+deltad*n);
    norm += gaussian_integral(d-deltad/2.0, d+deltad/2.0, m, s);
  }
  type = polar_distr_normal;
  return;
}

void polardistr::set_uniform(int np, double d1, double d2)
{
  this->np = np;
  this->d1 = d1;
  this->d2 = d2;
  mean = (d1-d2)/2.0;
  stdd = -1;
  range = (d2-d1);
  deltad = range/np;
  norm = (d2-d1);
  type = polar_distr_equipr;
}

void polardistr::set_random(int np, double d1, double d2)
{
  this->np = np;
  this->d1 = d1;
  this->d2 = d2;
  mean = (d1-d2)/2.0;
  stdd = -1;
  range = (d2-d1);
  deltad = range/np;
  norm = (d2-d1);
  type = polar_distr_random;
}

double polardistr::GetVal(int index)
{
  if (type == polar_distr_random) return d1+ranf( )*(d2-d1);
  if (np == 0) return d1;
  return d1+index*deltad;
}

double polardistr::GetWeight(int index)
{
  if (np == 0) return 1.0; // Just one value !

  double d = GetVal(index);
  switch (type)
  {
    case polar_distr_normal:
      return gaussian_integral(d-deltad/2.0, d+deltad/2.0, mean, stdd)/norm;
      break;
    case polar_distr_equipr:
    case polar_distr_random:
      return (deltad/norm);
      break;
    default:
      throw ("Unknown Angular Distribution");
      break;
  }
  return 0.0;
}

polarimetric_radar::polarimetric_radar(double freqGHz)
{
  this->freqGHz = freqGHz;
  wavelen_micr = 299792.458 / freqGHz;
  wavelen_mm = wavelen_micr * 0.001;
}

radar_observation::radar_observation(polarimetric_radar &p)
{
  fr = p.freqGHz;
  wl = p.wavelen_mm;
  ath_coeff = (0.5*wl*wl)/M_PI;
  ndiam = 50;
  alpha_distr.set_uniform(0, 0.0, 0.0);
  beta_distr.set_gaussian(180, 0.0, 180.0, 0.0, 5.0);
}

void radar_observation::setup(int ndiam, polardistr adistr, polardistr bdistr)
{
  this->ndiam = ndiam;
  alpha_distr = adistr;
  beta_distr = bdistr;
  return;
}

void radar_observation::tmatrix_calculate(hydrometeor &h, hydroPSD &p,
                                          double theta, double phi)
{
  zhh = 0;
  zvv = 0;
  zhv = 0;
  zvh = 0;
  ldrhv = 0;
  ldrvh = 0;
  rhohv = 0;
  dhv = 0;
  kdp = 0;
  ath = 0;
  athh = 0;
  athv = 0;
  atvh = 0;
  atvv = 0;
  phihv = 0;
  rrate = 0;
  wc = 0;

  long double z11 = 0, z12 = 0, z13 = 0, z14 = 0;
  long double z21 = 0, z22 = 0, z23 = 0, z24 = 0;
  long double z31 = 0, z32 = 0, z33 = 0, z34 = 0;
  long double z41 = 0, z42 = 0, z43 = 0, z44 = 0;
  std::complex <double> rho(0.0, 0.0);
  std::complex <double> rhodet(0.0, 0.0);
  std::complex <double> cshh2(0.0, 0.0);
  std::complex <double> csvv2(0.0, 0.0);
  std::complex <double> cshv2(0.0, 0.0);
  std::complex <double> csvh2(0.0, 0.0);
  std::complex <double> ckdp(0.0, 0.0);
  std::complex <double> cathh(0.0, 0.0);
  std::complex <double> catvv(0.0, 0.0);
  std::complex <double> cathv(0.0, 0.0);
  std::complex <double> catvh(0.0, 0.0);

  int nd1 = ndiam-1;
  // double *diams = p.regspaced(ndiam);
  double *diams = p.logspaced(ndiam);
  double *diam = new(double[nd1]);
  double *ddiams = new(double[nd1]);
  double pin = M_PI/180.0;
  int *ng = new(int[nd1]);
  scattering_particle *part[nd1];
  for (int id = 0; id < (nd1); id ++)
  {
    ddiams[id] = diams[id+1] - diams[id];
    diam[id]  = diams[id] + ddiams[id] / 2.0;
    switch (h.GetShape(diam[id]))
    {
      case particle_raindrop:
        part[id] = new raindrop(diam[id]/2.0, h.ShapeParameter(diam[id],p),
                             h.RefractionIndex(wl));
        break;
      case particle_cylinder:
        part[id] = new cylinder(diam[id]/2.0, h.ShapeParameter(diam[id],p),
                     radius_equal_volume_sphere, h.RefractionIndex(wl));
        break;
      case particle_spheroid:
        part[id] = new spheroid(diam[id]/2.0, h.ShapeParameter(diam[id],p),
                     radius_equal_volume_sphere, h.RefractionIndex(wl));
        break;
      case particle_chebyshev:
        part[id] = new chebyshev(diam[id]/2.0, h.ShapeParameter(diam[id],p),
                      h.GetPolynomialOrder(diam[id]),
                      radius_equal_volume_sphere, h.RefractionIndex(wl));
        break;
      default:
        throw("Error in shape: Unsupported");
        break;
    }
    ng[id] = h.BaseGaussianPoints(diam[id]);
  }
  delete [] diams;

  radiation rad(wl);
  Tmatrix TM[nd1];

  // Part to be executed in parallel
  #pragma omp parallel
  {
    for (int ipar = 0; ipar < nd1; ipar ++)
    {
      #pragma omp sections
      {
        #pragma omp section
        TM[ipar].Setup(*(part[ipar]), rad, 0.0001, ng[ipar]);
      }
    }
  }
  // Part to be executed in parallel

  for (int id = 0; id < nd1; id ++)
    delete part[id];

  AmplitudeMatrix A;
  PhaseMatrix P;

  for (int id = 0; id < nd1; id ++)
  {
    long double NDdD = ddiams[id]*p.ndd(diam[id]);
    int nalpha = alpha_distr.Points( );
    int nbeta = beta_distr.Points( );

    rrate += h.PrecipRate(diam[id], ddiams[id], p.ndd(diam[id]));
    wc += h.Water(diam[id], ddiams[id], p.ndd(diam[id]));

    for (int ia = 0; ia < nalpha; ia ++)
    {
      long double pdfa = alpha_distr.GetWeight(ia);

      double alpha = alpha_distr.GetVal(ia);

      for (int ib = 0; ib < nbeta; ib ++)
      {
        long double pdfo = beta_distr.GetWeight(ib);

        long double weight = pdfa*pdfo*NDdD;
        if (weight < 1e-20) continue;

        double beta = beta_distr.GetVal(ib);

        // Back propagation

        TM[id].GetAmplitudeAndPhaseMatrices(alpha, beta, theta, phi,
                                         180.0-theta, 180.0+phi, A, P);

        cshh2 += (A.s22*A.s22)*weight;
        csvv2 += (A.s11*A.s11)*weight;
        cshv2 += (A.s21*A.s21)*weight;
        csvh2 += (A.s12*A.s12)*weight;
        rho   += (conj(-A.s11)*A.s22)*weight;
        rhodet += sqrt(cshh2)*sqrt(csvv2);

        z11   += P.z11*weight;
        z12   += P.z12*weight;
        z13   += P.z13*weight;
        z14   += P.z14*weight;
        z21   += P.z21*weight;
        z22   += P.z22*weight;
        z23   += P.z23*weight;
        z24   += P.z24*weight;
        z31   += P.z31*weight;
        z32   += P.z32*weight;
        z33   += P.z33*weight;
        z34   += P.z34*weight;
        z41   += P.z41*weight;
        z42   += P.z42*weight;
        z43   += P.z43*weight;
        z44   += P.z44*weight;

        // Forward propagation

        TM[id].GetAmplitudeAndPhaseMatrices(alpha, beta, theta, phi,
                                         theta, phi, A, P);

        ckdp  += (A.s22-A.s11)*weight;
        cathh += A.s22*weight;
        catvv += A.s11*weight;
        cathv += A.s21*weight;
        catvh += A.s12*weight;
      }
    }

    ath += 4.343*1e-3*(-TM[id].GetQext( )*ath_coeff)*NDdD;
  }

  delete [] diam;
  delete [] ddiams;
  delete [] ng;

  long double zcost = (pow(wl, 4.0)*2*M_PI)/(pow(M_PI, 5)*
                       KappaSquare(h.RefractionIndex(wl)));

  if ((z11-z12-z21+z22) > 0.0)
    zhh = zcost*(z11-z12-z21+z22);
  if ((z11+z12+z21+z22) > 0.0)
    zvv = zcost*(z11+z12+z21+z22);
  if ((z11+z12-z21-z22)> 0.0)
    zhv = zcost*(z11+z12-z21-z22);
  if ((z11-z12+z21-z22) > 0.0)
    zvh = zcost*(z11-z12+z21-z22);
  if (fabs(zhh-zvv) < 1e-20)
    zvv = zhh;
  if ((z11+z12+z21+z22) > 0.0 &&
      (z11+z12-z21-z22) > 0.0 &&
      (z11-z12+z21-z22) > 0.0)
  {
    ldrhv=(z11+z12-z21-z22)/(z11+z12+z21+z22);
    ldrvh=(z11-z12+z21-z22)/(z11-z12-z21+z22);
  }

  if (((z11-z12-z21+z22)*(z11+z12+z21+z22)) > 0.0)
  {
    rhohv=1.0-real(rho/rhodet);
    dhv = acos(rhohv)/pin;
  }

  kdp=1e-3*(180.0/M_PI)*wl*real(ckdp);
  athh=1e-3*8.686*wl*imag(cathh);
  atvv=1e-3*8.686*wl*imag(catvv);
  athv=1e-3*8.686*wl*imag(cathv);
  atvh=1e-3*8.686*wl*imag(catvh);
  
  if (real(cshv2) > 0.0)
    phihv=0.5*atan(imag(cshv2)/real(cshv2));

  return;
}

#ifdef TESTME

#include <radarband.h>

int main(int argc, char *argv[])
{
  radarband rb;
  std::cout << "Using radar frequency of 9.6 GHz ("
    << rb.BandFromFrGHz(9.6) << " band)." << std::endl;
  polarimetric_radar radar(9.6);
  radar_observation obs(radar);

  polardistr alphad, betad;
  alphad.set_random(0, 0, 360);
  betad.set_gaussian(180, 0.0, 180.0, 0.0, 5.0);
  obs.setup(50, alphad, betad);

  hydro_rain rain(beard_chuang_1987, 293.15, 0.997);
  GeneralizedGammaPSD psd_rain(1000, 5, 2.0);
  psd_rain.SetExtremes(2.0, 7);
  obs.tmatrix_calculate(rain, psd_rain, 9.6, 0.0);

  betad.set_uniform(90, 0, 90);
  obs.setup(50, alphad, betad);
 
  ExponentialPSD psd_ice(25, 2.1);
  psd_ice.SetExtremes(0.1*2, 2.4*2);
  hydro_ice ice(45000, 0.001, 230, 0.5);
  obs.tmatrix_calculate(ice, psd_ice, 9.6, 0.0);

  return 0;
}

#endif
