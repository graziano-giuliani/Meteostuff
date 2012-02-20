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
#include <hydropsd.h>

using namespace cetemps;

int factorial(int num)
{
  int result=1;
  for (int i=1; i<=num; ++i)
    result*=i;
  return result;
}

double kummer(const double aval, const double zval)
{
  double eps(1.0e-15);      // Small threshold controling precision
  double z( fabs(zval) );   // We only allow positive values of 'z'
  double a( fabs(aval) );   // We only allow positive values of 'a'
  double den(a);            // Variable to store denominator
  double s(1.0);            // Result will be stored here
  double coef(1.0);         // Initialize coefficient
  while( fabs(coef) > eps)
  {
    coef = coef*z;   // Compute numerator
    den += 1.0;
    coef = coef/den; // Compute coefficient
    s += coef;       // Add new coefficient to result
  }
  return s;
}

double lower_gamma(const double a, const double z)
{
  double zp( fabs(z) );    // We only allow positive values of 'z'
  double ap( fabs(a) );    // We only allow positive values of 'a'
  double s( kummer(ap, zp) );
  return exp(log(zp)*ap) * exp(-zp) * s / ap;
}

double upper_gamma(const double a, const double z)
{
  return ( tgamma(a) - lower_gamma(a, z) );
}

ExponentialPSD::ExponentialPSD(double intercept, double slope)
{
  this->intercept = intercept;
  this->slope = slope;
}

inline double ExponentialPSD::ndd(double d)
{
  return (intercept*exp(-slope*d));
}

inline double ExponentialPSD::sum(double d1, double d2)
{
  return (intercept/slope*(exp(-slope*d1)-exp(-slope*d2)));
}

inline double ExponentialPSD::moment(int order)
{
  return (intercept*factorial(order)/pow(slope, order+1));
}

GeneralizedGammaPSD::GeneralizedGammaPSD(double intercept,
                                         double shape,
                                         double dzero)
{
  this->intercept = intercept;
  this->shape = shape;
  this->dzero = dzero;
  slope = (3.67+shape)/dzero;
  fmu = (6.0/pow(3.67, 4.0))*pow((3.67+shape),(shape+4.0))/tgamma(shape+4.0);
}

inline double GeneralizedGammaPSD::ndd(double d)
{
  return (intercept*fmu*pow(d/dzero, shape)*exp(-slope*d));
}

inline double GeneralizedGammaPSD::sum(double d1, double d2)
{
  return (intercept*fmu/pow(dzero, shape) *
          pow(-1,shape) / pow(-slope, shape+1) *
          (upper_gamma(1+shape,slope*d2)-upper_gamma(1+shape,slope*d1)));
}

inline double GeneralizedGammaPSD::moment(int order)
{
  return (intercept*fmu/pow(slope, order+2)*pow(dzero,order+2-shape)*
          tgamma(order+2));
}

WRFGammaPSD::WRFGammaPSD(double intercept,
    double shape, double slope)
{
  this->intercept = intercept;
  this->shape = shape;
  this->slope = slope;
}

double WRFGammaPSD::ndd(double d)
{
  double dm = d*0.001;
  return 0.001*(intercept*pow(dm, shape)*exp(-slope*dm));
}

double WRFGammaPSD::sum(double d1, double d2)
{
  throw "Not implemented";
}

double WRFGammaPSD::moment(int order)
{
  throw "Not implemented";
}

WRFDoubleGammaPSD::WRFDoubleGammaPSD(double q, double t,
                                     double kap0, double kap1,
                                     double lam0, double lam1, double shape)
{
  double loga, a, b, M2, M3, Mf1, second;
  const double sa[10] = { 5.065339, -0.062659, -3.032362,
                          0.029469, -0.000285,  0.31255,
                          0.000204,  0.003199,  0.0,
                         -0.015952 };
  const double sb[10] = { 0.476221, -0.015896,  0.165977,
                          0.007468, -0.000141,  0.060366,
                          0.000079,  0.000594,  0.0,
                         -0.003577 };
  const double bm_s = 2.0;
  const double cse = bm_s + 1.0;
  this->kap0 = kap0;
  this->kap1 = kap1;
  this->lam0 = lam0;
  this->lam1 = lam1;
  this->shape = shape;
  double tc = t-273.15;

  // Mass power law relations:  mass = am*D**bm, Snow content in kg/m3
  M2 = q/0.069;

  if (bm_s > 2.0-1.E-3 && bm_s < 2.0+1.E-3)
  {
    loga = sa[0] + sa[1]*tc + sa[2]*bm_s + sa[3]*tc*bm_s + sa[4]*tc*tc +
           sa[5]*bm_s*bm_s + sa[6]*tc*tc*bm_s + sa[7]*tc*bm_s*bm_s +
           sa[8]*tc*tc*tc + sa[9]*bm_s*bm_s*bm_s;
    a = pow(10.0, loga);
    b = sb[0] + sb[1]*tc + sb[2]*bm_s + sb[3]*tc*bm_s + sb[4]*tc*tc +
        sb[5]*bm_s*bm_s + sb[6]*tc*tc*bm_s + sb[7]*tc*bm_s*bm_s +
        sb[8]*tc*tc*tc + sb[9]*bm_s*bm_s*bm_s;
    second = pow(M2/a, 1.0/b);
  }
  else
    second = M2;

  loga = sa[0] + sa[1]*tc + sa[2]*cse + sa[3]*tc*cse + sa[4]*tc*tc +
         sa[5]*cse*cse + sa[6]*tc*tc*cse + sa[7]*tc*cse*cse +
         sa[8]*tc*tc*tc + sa[9]*cse*cse*cse;
  a = pow(10.0, loga);
  b = sb[0] + sb[1]*tc + sb[2]*cse + sb[3]*tc*cse + sb[4]*tc*tc +
      sb[5]*cse*cse + sb[6]*tc*tc*cse + sb[7]*tc*cse*cse +
      sb[8]*tc*tc*tc + sb[9]*cse*cse*cse;
  M3 = a * pow(second,b);
  Mf1 = M2/M3;
  mult = pow(M2,4.0)/pow(M3,3.0)*0.001;
  kap1 = kap1*pow(Mf1,shape);
  lam0 = lam0*Mf1;
  lam1 = lam1*Mf1;
}

double WRFDoubleGammaPSD::ndd(double d)
{
  double dm = d*0.001;
  return mult*(kap0*exp(-lam0*dm) + kap1*pow(dm,shape)*exp(-lam1*dm));
}

double WRFDoubleGammaPSD::sum(double d1, double d2)
{
  throw "Not implemented";
}

double WRFDoubleGammaPSD::moment(int order)
{
  throw "Not implemented";
}

void hydroPSD::SetExtremes(double d1, double d2)
{
  if (d1 < d2)
  {
    lower_bound = d1;
    upper_bound = d2;
  }
  else
  {
    lower_bound = d2;
    upper_bound = d1;
  }
  return;
}

double hydroPSD::minvalue( )
{
  return lower_bound;
}

double hydroPSD::maxvalue( )
{
  return upper_bound;
}

#ifdef TESTME

#include <iostream>

int main(int argc, char *argv[])
{
  hydroPSD *p;

  ExponentialPSD epsd(800, 2.7);
  p = (hydroPSD *) &epsd;

  std::cout << "N(2.0) = "       << p->ndd(2.0) << std::endl;
  std::cout << "SUM(2.0,4.0) = " << p->sum(2.0,4.0) << std::endl;
  std::cout << "M1 =           " << p->moment(1) << std::endl;
  std::cout << "M2 =           " << p->moment(2) << std::endl;
  std::cout << "M3 =           " << p->moment(3) << std::endl;

  GeneralizedGammaPSD gpsd(800, 3.0, 2.7);
  p = (hydroPSD *) &gpsd;

  std::cout << "N(2.0) = "       << p->ndd(2.0) << std::endl;
  std::cout << "SUM(2.0,4.0) = " << p->sum(2.0,4.0) << std::endl;
  std::cout << "M1 =           " << p->moment(1) << std::endl;
  std::cout << "M2 =           " << p->moment(2) << std::endl;
  std::cout << "M3 =           " << p->moment(3) << std::endl;
  return 0;
}

#endif
