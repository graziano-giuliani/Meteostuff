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

#include <cmath>

#include <psychrometric.h>
#include <meteosalute.h>

using namespace std;
using namespace Meteo;

static inline float metabolism(float t)
{
  return (-3.0909 * t + 203.64); // Base metabolism in function of T
}

float Meteosalute::pmv_hoppe_iso(float t, float rh, float wind, float mtrad,
                                 float iclo)
{
  const float eta = 0.01; // Mechanical efficiency
  const float age = 35.0; // Age
  const float mbody = 75.0; // Weigth in kg
  const float ht = 1.75; // Heigth in m
  const float tcl = 30.005;
  const int MAX_LOOP = 200;
  const int MAX_LOOP_HALF = MAX_LOOP / 2;
  const float tcl_eps = 0.05;
  const float eps = 0.97;
  const float sigm = 5.67e-8;
  float vpa, fcl, metm, metf, metb, h, aef, p1, tcl1, tcl2;
  float hc, diff, abhc, abtcl, difhc, tsk, esw, rsum, csum;
  float erel, eres, ed, load, ts;

  const float adu = 0.203 * powf(mbody, 0.425) * powf(ht, 0.725);
  const float metbf = 3.19 * powf(mbody, (3.0/4.0)) * (1.0 + 0.004* (30.0- age)
      +0.018 * ((ht * 100.0/ powf(mbody, (1.0/3.0))) - 42.1)); // Women
  const float metbm = 3.45 * powf(mbody, (3.0/4.0)) * (1.0 + 0.004* (30.0- age)
      +0.010 * ((ht * 100.0/ powf(mbody, (1.0/3.0))) - 43.4)); // Men

  vpa = func_dbrh2e(t, rh) / 10.0;
  fcl = 1.0 + iclo * 0.15;
  metb = metabolism(t);
  metf = metbf + metb;
  metm = metbm + metb;
  metb = (metf+metm)/2.0;
  h = metb * (1.0 - eta);
  aef = 0.71 * fcl * adu;

  p1 = 35.7 - 0.032 * (metb / (adu * 1.16))* (1 - eta);

  tcl1 = tcl;
  for (int i = 0; i < MAX_LOOP; i ++)
  {
    if (i < MAX_LOOP_HALF)
    {
      hc = 12.06 * sqrtf(wind);
      abhc = 0.0;
    }
    else
    {
      hc = 2.38 * powf(fabsf(tcl1 - t), 4.0);
      abhc = 0.6 * fabsf(powf((tcl1 - t), -0.75));
    }
    tcl2 = p1 - 0.155 * iclo * (3.94 * 0.00000001* fcl *(powf((tcl1 + 273.2),
                                                              4.0)- powf((mtrad
        + 273.2), 4.0))+fcl * hc* (tcl1 - t));
    diff = fabsf(tcl1 - tcl2);
    if (diff < tcl_eps)
      break;
    abtcl = -0.155 * iclo * (4.0 * 3.94* 0.00000001* fcl *powf((tcl1+ 273.2),
                                                               3.0) + fcl * hc
        - fcl *(tcl1 - t)* abhc)- 1.0;
    tcl1 = tcl1 - (tcl2 - tcl1) / abtcl;
    difhc = (12.06 * sqrtf(wind)) - (2.38 * (powf(fabsf(t - tcl1), 0.25)));
    if (difhc > 0.0&& i == MAX_LOOP_HALF)
      break;
  }
  tsk = 35.7 - (0.028 * h / adu);
  esw = 0.42 * adu * (h / adu - 58.08);
  esw = esw < 0.0 ? 0.0 : esw;
  rsum = aef * eps * sigm * (powf((tcl1 + 273.2), 4.0) - powf((mtrad + 273.2),
                                                              4.0));
  csum = adu * fcl * hc * (tcl1 - t);
  erel = 0.0023 * metb * (44.0 - 0.75 * vpa);
  eres = 0.0014 * metb * (34.0 - t);
  ed = 0.406 * adu * (1.92 * tsk - 25.3- 0.75 * vpa);
  load = (h - ed - erel - eres - esw - rsum - csum) / adu;
  ts = (0.303 * expf(-0.036 * (metb / adu)) + 0.028);
  return (ts * load);
}

float *Meteosalute::pmv_hoppe_iso(float *t, float *rh, float *wind,
                                  float *mtrad, float *iclo)
{
  if (vector_size <= 0)
    return 0;
  float *pmv = new float[vector_size];

  for (int i = 0; i < vector_size; i ++)
  {
    pmv[i] = pmv_hoppe_iso(t[i], rh[i], wind[i], mtrad[i], iclo[i]);
  }
  return pmv;
}

#ifdef TEST_ME

#include <iostream>

int main(int argc, char *argv[])
{
  const int V_SIZE = 1;
  Meteosalute m(V_SIZE);
  float t[V_SIZE] =
  { 23.0};
  float rh[V_SIZE] =
  { 85.0};
  float wind[V_SIZE] =
  { 5.0};
  float mtrad[V_SIZE] =
  { 300.0};
  float iclo[V_SIZE] =
  { 0.3};
  float *pmv = m.pmv_hoppe_iso(t, rh, wind, mtrad, iclo);

  float *res = pmv;
  for (int i = 0; i < V_SIZE; i ++)
  {
    std::cout << res[i] << std::endl;
  }
  return 0;
}

#endif
