/***************************************************************************
 *   Copyright (C) 2007 by Graziano Giuliani                               *
 *   graziano.giuliani at poste.it                                         *
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
// $Id: $

#include <meteosalute.h>

#include <cmath>

using namespace std;
using namespace Meteo;

// Calculation of UTCI
// based on calculations as presented at the MC/WG meeting
// of COST Action 730 in Eilat, 2008-09-10
// and decisions after the WG1 meeting in Stuttgart 2009-02-24
// alpha version released for public use after termination of COST Action 730
// Version a 0.002, October 2009
// Copyright (C) 2009  Peter Broede
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

float es(float ta)
{
  // Hardy, R.; ITS-90 Formulations for Vapor Pressure, Frostpoint
  // Temperature, Dewpoint Temperature and Enhancement Factors in the
  // Range -100 to 100 °C; 
  // Proceedings of Third International Symposium on Humidity and Moisture;
  // edited by National Physical Laboratory (NPL), London, 1998, pp. 214-221
  // http://www.thunderscientific.com/tech_info/reflibrary/its90formulas.pdf
  // (retrieved 2008-10-01)
  float es , tk;
  float g[8] = { -2.8365744E3,
                 -6.028076559E3,
                  1.954263612E1,
                 -2.737830188E-2,
                  1.6261698E-5,
                  7.0229056E-10,
                 -1.8680009E-13,
                  2.7150305 };
  tk = ta+273.15; // air temp in K
  es = g[7]*log(tk);
  for ( int i = 0; i < 7; i ++ )
    es = es+g[i]*pow(tk,i-2); 
  es = exp(es)*0.01; // convert Pa to hPa
  return es;
}

// The Universal Thermal Climate Index UTCI in Operational Use
// 
// For the past four years the European Union has funded within the COST
// Action 730 the development of the Universal Thermal Climate Index UTCI.
// This enforced the efforts of ISB Commission 6 on UTCI which already
// started in 2000. 
// 
// The term “universal” must be understood in terms of appropriate for
// all assessments of the outdoor thermal conditions in the major fields of
// human biometeorology such as public weather service, public health
// system, precautionary planning, and climate impact research in the
// health sector. UTCI should become an international standard based on
// recent scientific progress in human response related
// thermo-physiological modelling.
// 
// After accessible models of human thermoregulation had been evaluated,
// the advanced multi-node ‘Fiala’ thermoregulation model was selected,
// extensively validated, and extended for purposes of the project. In the
// next step a state-of-the-art adaptive clothing model was developed and
// integrated. This model considers (i) the behavioural adaptation of
// clothing insulation observed for the general urban population in
// relation to the actual environmental temperature, (ii) the distribution
// of the clothing over different body parts providing local insulation
// values for the different model segments, and (iii) the reduction of
// thermal and evaporative clothing resistances caused by wind and the
// movement of the wearer, who was assumed walking 4 km/h on the level.
// 
// UTCI was then developed following the concept of an equivalent
// temperature. This involved the definition of a reference environment
// with 50% relative humidity (but not exceeding 20 hPa), with still air
// and radiant temperature equalling air temperature, to which all other
// climatic conditions are compared. Equal physiological conditions are
// based on the equivalence of the dynamic physiological response predicted
// by the model for the actual and the reference environment. As this
// dynamic response is multidimensional (body core temperature, sweat rate,
// skin wettedness etc at different exposure times), a single dimensional
// strain index was calculated by principal component analysis. The UTCI
// equivalent temperature for a given combination of wind, radiation,
// humidity and air temperature is then defined as the air temperature of
// the reference environment which produces the same strain index value.
// 
// The associated assessment scale was derived from the simulated
// physiological responses and comprises ten thermal stress categories
// ranging from extreme cold stress to extreme heat stress.
// 
// As calculating the UTCI equivalent temperatures by running the
// thermoregulation model repeatedly could be too time-consuming for
// climate simulations and numerical weather forecasts, several options to
// speed up this calculation were considered. In the first instance
// polynomial regression equations predicting the UTCI equivalent
// temperature values are available as an operational procedure which is
// accessible both as software source code and executable program at the
// project’s website (  HYPERLINK "http://www.utci.org"  www.utci.org ).
// Comparisons to existing thermal stress/strain assessment procedures
// showed good conformity. However, in contrast to these procedures, UTCI
// is based on contemporary science.
// 
// The difficulties of different meteorological data levels (observations,
// numerical simulations, etc.), particularly with respect to the
// calculation of mean radiant temperature, were also assessed. Potential
// applications were identified in the fields of public weather services,
// public health systems, urban planning, tourism & recreation and climate
// impact research. It is recommended to run the UTCI model for the
// fundamental application in Numerical Weather Predictions and climate
// assessments operationally in Regional Specialised Meteorological Centres
// or Regional Climate Centres, respectively.
// 
// Freiburg, Germany, October 15, 2009
// 
// Prof Dr Gerd Jendritzky on behalf of the co- and WG chairs George
// Havenith, UK, Richard de Dear, AU, Philipp Weihs, AT, Ekaterina
// Batchvarova, BG, plus the MC members and experts.
// 
float *Meteosalute::utci(float *t,    // Temperature Celsius
                         float *tmrt, // Mean radiant temperature
                         float *wind, // Wind m/s
                         float *rh)   // Relative humidity %
{
  if (vector_size <= 0) return 0;
  float *utci = new float[vector_size];
  float ta, pa, va, e, dtm;

  for (int i = 0; i < vector_size; i ++)
  {
    if ( t[i] < -50.0 || t[i] > 50.0 ) continue;
    if ( tmrt[i] < t[i]-30.0 || tmrt[i] > t[i]+70.0 ) continue;
    if ( wind[i] < 0.5 || wind[i] > 30.0 ) continue;
    if ( rh[i] <= 0.0 || rh[i] >= 100.0 ) continue;
    ta = t[i];
    e = es(ta);
    pa = (e*rh[i]/100.0)/10.0; // use vapour pressure in kPa
    va = wind[i];
    dtm = tmrt[i] - ta;
    // computed by a 6th order approximating polynomial
    utci[i] = ta+ ( 6.07562052E-01 ) + ( -2.27712343E-02 ) * ta + 
            ( 8.06470249E-04 ) * ta*ta +
            ( -1.54271372E-04 ) * ta*ta*ta +
            ( -3.24651735E-06 ) * ta*ta*ta*ta +
            ( 7.32602852E-08 ) * ta*ta*ta*ta*ta +
            ( 1.35959073E-09 ) * ta*ta*ta*ta*ta*ta +
            ( -2.25836520E+00 ) * va +
            ( 8.80326035E-02 ) * ta*va +
            ( 2.16844454E-03 ) * ta*ta*va +
            ( -1.53347087E-05 ) * ta*ta*ta*va +
            ( -5.72983704E-07 ) * ta*ta*ta*ta*va +
            ( -2.55090145E-09 ) * ta*ta*ta*ta*ta*va +
            ( -7.51269505E-01 ) * va*va +
            ( -4.08350271E-03 ) * ta*va*va +
            ( -5.21670675E-05 ) * ta*ta*va*va +
            ( 1.94544667E-06 ) * ta*ta*ta*va*va +
            ( 1.14099531E-08 ) * ta*ta*ta*ta*va*va +
            ( 1.58137256E-01 ) * va*va*va +
            ( -6.57263143E-05 ) * ta*va*va*va +
            ( 2.22697524E-07 ) * ta*ta*va*va*va +
            ( -4.16117031E-08 ) * ta*ta*ta*va*va*va +
            ( -1.27762753E-02 ) * va*va*va*va +
            ( 9.66891875E-06 ) * ta*va*va*va*va +
            ( 2.52785852E-09 ) * ta*ta*va*va*va*va +
            ( 4.56306672E-04 ) * va*va*va*va*va +
            ( -1.74202546E-07 ) * ta*va*va*va*va*va +
            ( -5.91491269E-06 ) * va*va*va*va*va*va +
            ( 3.98374029E-01 ) * dtm +
            ( 1.83945314E-04 ) * ta*dtm +
            ( -1.73754510E-04 ) * ta*ta*dtm +
            ( -7.60781159E-07 ) * ta*ta*ta*dtm +
            ( 3.77830287E-08 ) * ta*ta*ta*ta*dtm +
            ( 5.43079673E-10 ) * ta*ta*ta*ta*ta*dtm +
            ( -2.00518269E-02 ) * va*dtm +
            ( 8.92859837E-04 ) * ta*va*dtm +
            ( 3.45433048E-06 ) * ta*ta*va*dtm +
            ( -3.77925774E-07 ) * ta*ta*ta*va*dtm +
            ( -1.69699377E-09 ) * ta*ta*ta*ta*va*dtm +
            ( 1.69992415E-04 ) * va*va*dtm +
            ( -4.99204314E-05 ) * ta*va*va*dtm +
            ( 2.47417178E-07 ) * ta*ta*va*va*dtm +
            ( 1.07596466E-08 ) * ta*ta*ta*va*va*dtm +
            ( 8.49242932E-05 ) * va*va*va*dtm +
            ( 1.35191328E-06 ) * ta*va*va*va*dtm +
            ( -6.21531254E-09 ) * ta*ta*va*va*va*dtm +
            ( -4.99410301E-06 ) * va*va*va*va*dtm +
            ( -1.89489258E-08 ) * ta*va*va*va*va*dtm +
            ( 8.15300114E-08 ) * va*va*va*va*va*dtm +
            ( 7.55043090E-04 ) * dtm*dtm +
            ( -5.65095215E-05 ) * ta*dtm*dtm +
            ( -4.52166564E-07 ) * ta*ta*dtm*dtm +
            ( 2.46688878E-08 ) * ta*ta*ta*dtm*dtm +
            ( 2.42674348E-10 ) * ta*ta*ta*ta*dtm*dtm +
            ( 1.54547250E-04 ) * va*dtm*dtm +
            ( 5.24110970E-06 ) * ta*va*dtm*dtm +
            ( -8.75874982E-08 ) * ta*ta*va*dtm*dtm +
            ( -1.50743064E-09 ) * ta*ta*ta*va*dtm*dtm +
            ( -1.56236307E-05 ) * va*va*dtm*dtm +
            ( -1.33895614E-07 ) * ta*va*va*dtm*dtm +
            ( 2.49709824E-09 ) * ta*ta*va*va*dtm*dtm +
            ( 6.51711721E-07 ) * va*va*va*dtm*dtm +
            ( 1.94960053E-09 ) * ta*va*va*va*dtm*dtm +
            ( -1.00361113E-08 ) * va*va*va*va*dtm*dtm +
            ( -1.21206673E-05 ) * dtm*dtm*dtm +
            ( -2.18203660E-07 ) * ta*dtm*dtm*dtm +
            ( 7.51269482E-09 ) * ta*ta*dtm*dtm*dtm +
            ( 9.79063848E-11 ) * ta*ta*ta*dtm*dtm*dtm +
            ( 1.25006734E-06 ) * va*dtm*dtm*dtm +
            ( -1.81584736E-09 ) * ta*va*dtm*dtm*dtm +
            ( -3.52197671E-10 ) * ta*ta*va*dtm*dtm*dtm +
            ( -3.36514630E-08 ) * va*va*dtm*dtm*dtm +
            ( 1.35908359E-10 ) * ta*va*va*dtm*dtm*dtm +
            ( 4.17032620E-10 ) * va*va*va*dtm*dtm*dtm +
            ( -1.30369025E-09 ) * dtm*dtm*dtm*dtm +
            ( 4.13908461E-10 ) * ta*dtm*dtm*dtm*dtm +
            ( 9.22652254E-12 ) * ta*ta*dtm*dtm*dtm*dtm +
            ( -5.08220384E-09 ) * va*dtm*dtm*dtm*dtm +
            ( -2.24730961E-11 ) * ta*va*dtm*dtm*dtm*dtm +
            ( 1.17139133E-10 ) * va*va*dtm*dtm*dtm*dtm +
            ( 6.62154879E-10 ) * dtm*dtm*dtm*dtm*dtm +
            ( 4.03863260E-13 ) * ta*dtm*dtm*dtm*dtm*dtm +
            ( 1.95087203E-12 ) * va*dtm*dtm*dtm*dtm*dtm +
            ( -4.73602469E-12 ) * dtm*dtm*dtm*dtm*dtm*dtm +
            ( 5.12733497E+00 ) * pa +
            ( -3.12788561E-01 ) * ta*pa +
            ( -1.96701861E-02 ) * ta*ta*pa +
            ( 9.99690870E-04 ) * ta*ta*ta*pa +
            ( 9.51738512E-06 ) * ta*ta*ta*ta*pa +
            ( -4.66426341E-07 ) * ta*ta*ta*ta*ta*pa +
            ( 5.48050612E-01 ) * va*pa +
            ( -3.30552823E-03 ) * ta*va*pa +
            ( -1.64119440E-03 ) * ta*ta*va*pa +
            ( -5.16670694E-06 ) * ta*ta*ta*va*pa +
            ( 9.52692432E-07 ) * ta*ta*ta*ta*va*pa +
            ( -4.29223622E-02 ) * va*va*pa +
            ( 5.00845667E-03 ) * ta*va*va*pa +
            ( 1.00601257E-06 ) * ta*ta*va*va*pa +
            ( -1.81748644E-06 ) * ta*ta*ta*va*va*pa +
            ( -1.25813502E-03 ) * va*va*va*pa +
            ( -1.79330391E-04 ) * ta*va*va*va*pa +
            ( 2.34994441E-06 ) * ta*ta*va*va*va*pa +
            ( 1.29735808E-04 ) * va*va*va*va*pa +
            ( 1.29064870E-06 ) * ta*va*va*va*va*pa +
            ( -2.28558686E-06 ) * va*va*va*va*va*pa +
            ( -3.69476348E-02 ) * dtm*pa +
            ( 1.62325322E-03 ) * ta*dtm*pa +
            ( -3.14279680E-05 ) * ta*ta*dtm*pa +
            ( 2.59835559E-06 ) * ta*ta*ta*dtm*pa +
            ( -4.77136523E-08 ) * ta*ta*ta*ta*dtm*pa +
            ( 8.64203390E-03 ) * va*dtm*pa +
            ( -6.87405181E-04 ) * ta*va*dtm*pa +
            ( -9.13863872E-06 ) * ta*ta*va*dtm*pa +
            ( 5.15916806E-07 ) * ta*ta*ta*va*dtm*pa +
            ( -3.59217476E-05 ) * va*va*dtm*pa +
            ( 3.28696511E-05 ) * ta*va*va*dtm*pa +
            ( -7.10542454E-07 ) * ta*ta*va*va*dtm*pa +
            ( -1.24382300E-05 ) * va*va*va*dtm*pa +
            ( -7.38584400E-09 ) * ta*va*va*va*dtm*pa +
            ( 2.20609296E-07 ) * va*va*va*va*dtm*pa +
            ( -7.32469180E-04 ) * dtm*dtm*pa +
            ( -1.87381964E-05 ) * ta*dtm*dtm*pa +
            ( 4.80925239E-06 ) * ta*ta*dtm*dtm*pa +
            ( -8.75492040E-08 ) * ta*ta*ta*dtm*dtm*pa +
            ( 2.77862930E-05 ) * va*dtm*dtm*pa +
            ( -5.06004592E-06 ) * ta*va*dtm*dtm*pa +
            ( 1.14325367E-07 ) * ta*ta*va*dtm*dtm*pa +
            ( 2.53016723E-06 ) * va*va*dtm*dtm*pa +
            ( -1.72857035E-08 ) * ta*va*va*dtm*dtm*pa +
            ( -3.95079398E-08 ) * va*va*va*dtm*dtm*pa +
            ( -3.59413173E-07 ) * dtm*dtm*dtm*pa +
            ( 7.04388046E-07 ) * ta*dtm*dtm*dtm*pa +
            ( -1.89309167E-08 ) * ta*ta*dtm*dtm*dtm*pa +
            ( -4.79768731E-07 ) * va*dtm*dtm*dtm*pa +
            ( 7.96079978E-09 ) * ta*va*dtm*dtm*dtm*pa +
            ( 1.62897058E-09 ) * va*va*dtm*dtm*dtm*pa +
            ( 3.94367674E-08 ) * dtm*dtm*dtm*dtm*pa +
            ( -1.18566247E-09 ) * ta*dtm*dtm*dtm*dtm*pa +
            ( 3.34678041E-10 ) * va*dtm*dtm*dtm*dtm*pa +
            ( -1.15606447E-10 ) * dtm*dtm*dtm*dtm*dtm*pa +
            ( -2.80626406E+00 ) * pa*pa +
            ( 5.48712484E-01 ) * ta*pa*pa +
            ( -3.99428410E-03 ) * ta*ta*pa*pa +
            ( -9.54009191E-04 ) * ta*ta*ta*pa*pa +
            ( 1.93090978E-05 ) * ta*ta*ta*ta*pa*pa +
            ( -3.08806365E-01 ) * va*pa*pa +
            ( 1.16952364E-02 ) * ta*va*pa*pa +
            ( 4.95271903E-04 ) * ta*ta*va*pa*pa +
            ( -1.90710882E-05 ) * ta*ta*ta*va*pa*pa +
            ( 2.10787756E-03 ) * va*va*pa*pa +
            ( -6.98445738E-04 ) * ta*va*va*pa*pa +
            ( 2.30109073E-05 ) * ta*ta*va*va*pa*pa +
            ( 4.17856590E-04 ) * va*va*va*pa*pa +
            ( -1.27043871E-05 ) * ta*va*va*va*pa*pa +
            ( -3.04620472E-06 ) * va*va*va*va*pa*pa +
            ( 5.14507424E-02 ) * dtm*pa*pa +
            ( -4.32510997E-03 ) * ta*dtm*pa*pa +
            ( 8.99281156E-05 ) * ta*ta*dtm*pa*pa +
            ( -7.14663943E-07 ) * ta*ta*ta*dtm*pa*pa +
            ( -2.66016305E-04 ) * va*dtm*pa*pa +
            ( 2.63789586E-04 ) * ta*va*dtm*pa*pa +
            ( -7.01199003E-06 ) * ta*ta*va*dtm*pa*pa +
            ( -1.06823306E-04 ) * va*va*dtm*pa*pa +
            ( 3.61341136E-06 ) * ta*va*va*dtm*pa*pa +
            ( 2.29748967E-07 ) * va*va*va*dtm*pa*pa +
            ( 3.04788893E-04 ) * dtm*dtm*pa*pa +
            ( -6.42070836E-05 ) * ta*dtm*dtm*pa*pa +
            ( 1.16257971E-06 ) * ta*ta*dtm*dtm*pa*pa +
            ( 7.68023384E-06 ) * va*dtm*dtm*pa*pa +
            ( -5.47446896E-07 ) * ta*va*dtm*dtm*pa*pa +
            ( -3.59937910E-08 ) * va*va*dtm*dtm*pa*pa +
            ( -4.36497725E-06 ) * dtm*dtm*dtm*pa*pa +
            ( 1.68737969E-07 ) * ta*dtm*dtm*dtm*pa*pa +
            ( 2.67489271E-08 ) * va*dtm*dtm*dtm*pa*pa +
            ( 3.23926897E-09 ) * dtm*dtm*dtm*dtm*pa*pa +
            ( -3.53874123E-02 ) * pa*pa*pa +
            ( -2.21201190E-01 ) * ta*pa*pa*pa +
            ( 1.55126038E-02 ) * ta*ta*pa*pa*pa +
            ( -2.63917279E-04 ) * ta*ta*ta*pa*pa*pa +
            ( 4.53433455E-02 ) * va*pa*pa*pa +
            ( -4.32943862E-03 ) * ta*va*pa*pa*pa +
            ( 1.45389826E-04 ) * ta*ta*va*pa*pa*pa +
            ( 2.17508610E-04 ) * va*va*pa*pa*pa +
            ( -6.66724702E-05 ) * ta*va*va*pa*pa*pa +
            ( 3.33217140E-05 ) * va*va*va*pa*pa*pa +
            ( -2.26921615E-03 ) * dtm*pa*pa*pa +
            ( 3.80261982E-04 ) * ta*dtm*pa*pa*pa +
            ( -5.45314314E-09 ) * ta*ta*dtm*pa*pa*pa +
            ( -7.96355448E-04 ) * va*dtm*pa*pa*pa +
            ( 2.53458034E-05 ) * ta*va*dtm*pa*pa*pa +
            ( -6.31223658E-06 ) * va*va*dtm*pa*pa*pa +
            ( 3.02122035E-04 ) * dtm*dtm*pa*pa*pa +
            ( -4.77403547E-06 ) * ta*dtm*dtm*pa*pa*pa +
            ( 1.73825715E-06 ) * va*dtm*dtm*pa*pa*pa +
            ( -4.09087898E-07 ) * dtm*dtm*dtm*pa*pa*pa +
            ( 6.14155345E-01 ) * pa*pa*pa*pa +
            ( -6.16755931E-02 ) * ta*pa*pa*pa*pa +
            ( 1.33374846E-03 ) * ta*ta*pa*pa*pa*pa +
            ( 3.55375387E-03 ) * va*pa*pa*pa*pa +
            ( -5.13027851E-04 ) * ta*va*pa*pa*pa*pa +
            ( 1.02449757E-04 ) * va*va*pa*pa*pa*pa +
            ( -1.48526421E-03 ) * dtm*pa*pa*pa*pa +
            ( -4.11469183E-05 ) * ta*dtm*pa*pa*pa*pa +
            ( -6.80434415E-06 ) * va*dtm*pa*pa*pa*pa +
            ( -9.77675906E-06 ) * dtm*dtm*pa*pa*pa*pa +
            ( 8.82773108E-02 ) * pa*pa*pa*pa*pa +
            ( -3.01859306E-03 ) * ta*pa*pa*pa*pa*pa +
            ( 1.04452989E-03 ) * va*pa*pa*pa*pa*pa +
            ( 2.47090539E-04 ) * dtm*pa*pa*pa*pa*pa +
            ( 1.48348065E-03 ) * pa*pa*pa*pa*pa*pa;
  }
  return utci;
}

#ifdef TEST_ME

#include <iostream>

int main(int argc, char *argv[])
{
  const int V_SIZE = 1;
  Meteosalute m(V_SIZE);
  float t[V_SIZE] = { 30.0};
  float tmr[V_SIZE] = { 40.0};
  float wind[V_SIZE] = { 2.0};
  float rh[V_SIZE] = { 88.0};
  float *newvar = m.utci(t,tmr,wind,rh);

  float *res = newvar;
  for (int i = 0; i < V_SIZE; i ++)
  {
    std::cout << res[i] << std::endl;
  }
  return 0;
}

#endif
