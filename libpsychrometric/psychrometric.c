/************************************************************************/
/*                            Copyright 1991                            */
/*                        University of Kentucky                        */
/*                Department of Agricultural Engineering                */
/*                              S.H. Zhang                              */
/*                              R.S. Gates                              */
/*                                                                      */
/*          Users may copy freely if this notice is remained            */
/*                  Suggestions:  gates @ engr.uky.edu                  */
/************************************************************************/
#include <math.h>
#include <stdio.h>

#include "psychrometric.h"

/*----------------------------------------------------------------------*/
/* Calculates the dew_point temperature [C].                            */
/* Based on formula in the ASHRAE Handbook of Fundamentals.             */
/* Applise for temperatures between -60 and +70 [C].                    */
/* Inputs are dry-bulb temperature [C] and water vapor pressure [kPa].  */
/*----------------------------------------------------------------------*/
double dewpt(double t, double kpa) {
  double p, dpt;
  p=log(1000.0*kpa);
  if (t>=0.0)
    dpt=(-60.45)+7.0322*p+.37*p*p;
  else
    dpt=(-35.957)-1.8726*p+1.1689*p*p;
  return (dpt);
}

/*----------------------------------------------------------------------*/
/* Calculates partial pressure of water vapor at saturation [kPa] .     */
/* Based on the formulae provided by ASHRAE Handbook of Fundamentals.   */
/* Applies for temperatures between -40 and +120 [C].                   */
/* input: tdb---dry_bulb temperature  [C]                               */
/* output: function=patial saturation pressure  [kPa]                   */
/*----------------------------------------------------------------------*/
double psat(double tdb) {
  double t, y;

  t=tdb+273.15;
  if (t>273.15) {
    y=exp(-5800.2206/t+1.3914993-.048640239*t+.41764768e-4*pow(t, 2.0)
        -.14452093e-7*pow(t, 3.0)+6.5459673*log(t))/1000.0;
  } else {
    y=exp(-5674.5359/t+6.3925247-.9677843e-2*t+.62215701e-6*pow(t, 2.0)
        +.20747825e-8*pow(t, 3.0)-.9484024e-12*pow(t, 4.0)+4.1635019
        *log(t))/1000.0;
  }
  return (y);
}

/*----------------------------------------------------------------------*/
/* calculates the enthalpy of air/water vapor mixture [kJ/kg]           */
/* Based on furmula in the ASHRAE Handbook of Fundamentals.             */
/* Inputs are dry-bulb temperature [C] and humidity ratio (unitless).   */
/*----------------------------------------------------------------------*/
double enthal(double tdb, double w) {
  double y;
  y=1.006*tdb+w*(2501.0+1.805*tdb);
  return (y);
}

/*----------------------------------------------------------------------*/
/* spvol() gets the specific volume of air/water vapor mixture [m^3/kg] */
/* input: tc--- dry_bulb temperature [C]                                */
/*        w---- humidity ratio                                          */
/*        pa---- air pressure  [kPa]                                    */
/* output: function = specific volume                                   */
/*----------------------------------------------------------------------*/
double spvol(double tc, double w, double pa) {
  double t, y;
  t=tc+273.16;
  y=(287.055*t*(1+1.6078*w))/(pa*1000.0);
  return (y);
}

/*----------------------------------------------------------------------*/
/* calculates the humidity ratio of air/water vapor mixture.            */
/* Based on formula in the ASHRAE Handbook of Fundamentals.             */
/* Inputs is water vapor pressure [kPa]. Assume aipr=101.325 [kPa]      */
/* output : function = humidity ratio                                   */
/*----------------------------------------------------------------------*/
double humrat(double p, double pa) {
  double y;
  y=.62198*p/(pa-p);
  return (y);
}

/*----------------------------------------------------------------------*/
/*  wetblb()                                                            */
/* It solves for the wet-bulb temperature iteratively,using enthalpy.   */
/* Based on formula in the ASHRAE Gandbook of Fundamentals.             */
/* Inputs:dry-blub temperature [C]. humidity ratio (unitless)           */
/* enthalpy [kJ/kg]. and air pressure [kPa].                            */
/* Wet-bulb temperature is calculated to the nearest .01 C.             */
/*----------------------------------------------------------------------*/
double wetblb(double tstar, double w, double hh, double pa) {
  double t, dum1, dum3, hwstar, p, wstar, hstar, twb;

  t=tstar+1.0;
  dum3=10.0;
  do {
    dum3=(-dum3/10.0);
    do {
      t=t+dum3;
      hwstar=4.186*t;
      p=psat(t);
      wstar=humrat(p, pa);
      hstar=enthal(t, wstar);
      dum1=hstar-hh-(wstar-w)*hwstar;
    } while (dum1*dum3<=0.0);
  } while (fabs(dum3)>=.000001);

  twb=t;
  return twb;
}

/*----------------------------------------------------------------------*/
/* Calculates relative humidity (decimal) of air/water vapor mixture.   */
/* based on formula in the ASHRAE Handbook of Fundamentals .            */
/* Inputs are degree of saturation---mu (decimal),                      */
/*                      and water vapor pressure [kPa].                 */
/*----------------------------------------------------------------------*/
double relhum(double mu, double pws, double pa) {
  double rh;

  rh=mu/(1.0-(1.0-mu)*pws/pa);
  return (rh);
}

/*----------------------------------------------------------------------*/
/* Calculates the degree of sat. (decimal) of air/water vapor mixture . */
/* Based on formula in the ASHRAE Handbook of Fundamentals.             */
/* Inputs are relative humidity (decimal), water vapor pressure [kPa].  */
/*----------------------------------------------------------------------*/
double degsat(double rh, double pws, double pa) {
  double mu;
  mu=rh*(1.0-pws/pa)/(1.0-rh*pws/pa);
  return (mu);
}

/*----------------------------------------------------------------------*/
/* Inputs are : dry-bulb temperature [C] and relative humidity (decimal)*/
/* All 11 properties are calculated and stored in sp sturcture.         */
/*----------------------------------------------------------------------*/
void db_rh(double tdb, double rh, t_psyprops *sp) {
  sp->tdb=tdb;
  sp->rh=rh;
  sp->ppsat=psat(tdb);
  sp->pp=rh*sp->ppsat;
  sp->w=humrat(sp->pp, sp->airp);
  sp->tdp=dewpt(sp->tdb, sp->pp);
  sp->h=enthal(sp->tdb, sp->w);
  sp->spvol=spvol(sp->tdb, sp->w, sp->airp);
  sp->dens=(1.0/sp->spvol)*(1.0+sp->w);
  sp->twb=wetblb(sp->tdb, sp->w, sp->h, sp->airp);
  sp->mu=degsat(sp->rh, sp->ppsat, sp->airp);
  return;
}

/*----------------------------------------------------------------------*/
/* Inputs are dry-bulb temperature [C] and humidity ratio.              */
/* All 11 properties are calculated and stored in sp structure.         */
/*----------------------------------------------------------------------*/
void db_w(double tdb1, double w1, t_psyprops *sp) {
  double ws1;

  sp->tdb=tdb1;
  sp->ppsat=psat(tdb1);
  ws1=.62198*sp->ppsat/(sp->airp - sp->ppsat);

  if (ws1<w1) {
    sp->w=ws1;
    sp->h=enthal(tdb1, sp->w);
    sp->spvol=spvol(tdb1, sp->w, sp->airp);
    sp->dens=(1.0/sp->spvol)*(1.0+sp->w);
    sp->mu=1.0;
    sp->rh=1.0;
    sp->twb=tdb1;
    sp->tdp=tdb1;
    sp->pp=sp->ppsat;
  } else {
    sp->w=w1;
    sp->mu=w1/ws1;
    sp->h=enthal(tdb1, w1);
    sp->spvol=spvol(tdb1, w1, sp->airp);
    sp->dens=(1.0+w1)/sp->spvol;
    sp->rh=relhum(sp->mu, sp->ppsat, sp->airp);
    sp->twb=wetblb(tdb1, w1, sp->h, sp->airp);
    sp->pp=sp->rh*sp->ppsat;
    sp->tdp=dewpt(tdb1, sp->pp);
  }
  return;
}

/*---------------------------------------------------------------------  */
/* Inputs are :dry-bulb temperature[C] and wet-bulb temperature[C].      */
/* All 11 properties are calculated and stored in sp structure.          */
/*---------------------------------------------------------------------- */
void db_wb(double tdb1, double twb1, t_psyprops *sp) {
  double dummy1, dummy2;

  if (twb1>tdb1)
    sp->twb=tdb1;
  else
    sp->twb=twb1;
  dummy2=psat(twb1);
  dummy1=humrat(dummy2, sp->airp);
  sp->w=((2501.0-2.381*twb1)*dummy1-tdb1+twb1)/(2501.0+1.805*tdb1-4.186*twb1);
  /* see ASHRAE Fundamentals , ch.5   */
  if (sp->w<0.0)
    sp->w=.00000001; /* check for intersection */
  sp->tdb=tdb1;
  sp->h=enthal(tdb1, sp->w);
  sp->spvol=spvol(tdb1, sp->w, sp->airp);
  sp->dens=(1.0/sp->spvol)*(1.0+sp->w);
  sp->ppsat=psat(tdb1);
  sp->rh=(sp->airp*sp->w/sp->ppsat)/(sp->w+.62198);
  sp->pp=sp->rh*sp->ppsat;
  sp->mu=degsat(sp->rh, sp->ppsat, sp->airp);
  sp->tdp=dewpt(tdb1, sp->pp);
  return;
}

/*-------------------------------------------------------------------------*/
/* Inputs are dry-bulb and dew point temperature [C].                      */
/* All 11 properties are calculated and stored in sp sturcture.            */
/*-------------------------------------------------------------------------*/
void db_dp(double tdb1, double tdp1, t_psyprops *sp) {
  sp->tdp=tdp1;
  sp->ppsat=psat(tdb1);
  sp->tdb=tdb1;
  if (tdp1<=0.0)
    sp->pp=(exp((-7.0322+sqrt(49.452-1.48*(-60.45-tdp1)))/2.3378))/1000.0;
  else
    sp->pp=(exp((1.8726+sqrt(3.05066+4.6757*(35.957+tdp1)))/2.3378))/1000.0;

  if (tdb1<tdp1) {
    sp->pp=sp->ppsat;
    sp->tdp=tdb1;
  }

  sp->rh=sp->pp/sp->ppsat;
  if (sp->rh>1.0)
    sp->rh=1.0;
  sp->w=humrat(sp->pp, sp->airp);
  sp->h=enthal(tdb1, sp->w);
  sp->mu=degsat(sp->rh, sp->ppsat, sp->airp);
  sp->spvol=spvol(tdb1, sp->w, sp->airp);
  sp->dens=(1.0/sp->spvol)*(1.0+sp->w);
  sp->twb=wetblb(tdb1, sp->w, sp->h, sp->airp);
  return;
}

/*----------------------------------------------------------------------*/
/* Inputs are dry-bulb temperature[C] and enthalpy [kg/kJ].             */
/* All 11 properties are calculated and stored in sp structure .        */
/*----------------------------------------------------------------------*/
void db_h(double tdb1, double h1, t_psyprops *sp) {
  sp->tdb=tdb1;
  sp->h=h1;
  sp->ppsat=psat(tdb1);
  sp->w=(h1-1.006*tdb1)/(2501.0+1.805*tdb1);
  if (sp->w>0.0) {
    sp->spvol=spvol(tdb1, sp->w, sp->airp);
    sp->dens=(1.0/sp->spvol)*(1.0+sp->w);
    sp->rh=(sp->w*sp->airp/sp->ppsat)/(sp->w+.62198);
    if (sp->rh>1.0) {
      fprintf(stderr, "The ENTHAL value is too high.\n");
      fprintf(stderr, "It does not intersect the dry-bulb temp. line.\n");
    } else {
      sp->pp=sp->rh*sp->ppsat;
      sp->twb=wetblb(tdb1, sp->w, h1, sp->airp);
      sp->tdp=dewpt(tdb1, sp->pp);
      sp->mu=degsat(sp->rh, sp->ppsat, sp->airp);
    }
  } else {
    fprintf(stderr, "The ENTHAL value is too low.\n");
    fprintf(stderr, "It does not intersect the dry-bulb temperature line.\n");
  }
  return;
}

/*----------------------------------------------------------------------*/
/* Inputs are: relative humidity (decimal) and humidity ratio.          */
/* All 11 properties are calculated and stored in sp structure.         */
/*----------------------------------------------------------------------*/
void rh_w(double rh1, double w1, t_psyprops *sp) {
  double dummy1, dummy2;

  sp->rh=rh1;
  sp->w=w1;
  sp->ppsat=(sp->airp*w1/rh1)/(w1+.62198);
  sp->pp=sp->ppsat*rh1;
  sp->tdb=273.16;
  dummy2=sp->tdb;

  if (sp->ppsat<0.6112) {
    do {
      dummy1=dummy2;
      dummy2=6238.64/(24.2779-.344438*log(dummy1)-log(sp->ppsat));
    } while (fabs(fabs(dummy1)-fabs(dummy2))>=.01);
  } else {
    do {
      dummy1=dummy2;
      dummy2=7511.52/(89.63121-log(sp->ppsat)+.02399897*dummy1
          -1.1654551e-5*pow(dummy1, 2.0)-1.2810336e-8
          *pow(dummy1, 3.0)+2.0998405e-11*pow(dummy1, 4.0)-12.150799
          *log(dummy1));
    } while (fabs(fabs(dummy1)-fabs(dummy2))>=.01);
  }
  sp->tdb=dummy2-273.16;
  if (sp->tdb<-40.0)
    sp->tdb=(-40.0);
  sp->h=enthal(sp->tdb, sp->w);
  sp->tdp=dewpt(sp->tdb, sp->pp);
  sp->spvol=spvol(sp->tdb, sp->w, sp->airp);
  sp->dens=(1.0/sp->spvol)*(1.0+sp->w);
  sp->twb=wetblb(sp->tdb, sp->w, sp->h, sp->airp);
  sp->mu=degsat(sp->rh, sp->ppsat, sp->airp);
  return;
}

/*-----------------------------------------------------------------------*/
/* Inputs are: relative humidity (decimal) and wet-bulb temperature [C]. */
/* All 11 properties are calculated and stored in sp structure.          */
/*-----------------------------------------------------------------------*/
void rh_wb(double rh1, double twb1, t_psyprops *sp) {
  double dummy1, dummy2, dummy3, wstar;

  sp->rh=rh1;
  sp->twb=twb1;
  dummy2=twb1-1.0;
  dummy3=(-10.0);
  sp->ppsat=psat(twb1);
  wstar=humrat(sp->ppsat, sp->airp);

  do {
    dummy3=(-dummy3/10.0);
    do {
      dummy2=dummy2+dummy3;
      sp->w=((2501.0-2.381*twb1)*wstar-(dummy2-twb1))/(2501.0+1.805
          *dummy2-4.186*twb1);
      sp->pp=sp->airp*sp->w/(.62198+sp->w);
      dummy1=sp->pp/psat(dummy2);
    } while (dummy3*(dummy1-rh1)>=0.0);
  } while (fabs(dummy3)>=1.0e-6);

  sp->tdb=dummy2;
  sp->ppsat=psat(sp->tdb);
  sp->pp=rh1*sp->ppsat;
  sp->tdp=dewpt(sp->tdb, sp->pp);
  sp->w=humrat(sp->pp, sp->airp);
  sp->h=enthal(sp->tdb, sp->w);
  sp->spvol=spvol(sp->tdb, sp->w, sp->airp);
  sp->dens=(1.0/sp->spvol)*(1.0+sp->w);
  sp->mu=degsat(rh1, sp->ppsat, sp->airp);
  return;
}

/*----------------------------------------------------------------------*/
/* Inputs are: relative humidity (decimal) and dew-point temperature [C]*/
/* All 11 properties are calculated and stored in sp structure .        */
/*----------------------------------------------------------------------*/
void rh_dp(double rh1, double tdp1, t_psyprops *sp) {
  double dummy1;

  sp->rh=rh1;
  sp->tdp=tdp1;

  if (tdp1<=0.0)
    sp->pp=(exp((-7.0322+sqrt(49.452-1.48*(-60.45-tdp1)))/.74))/1000.0;
  else
    sp->pp=(exp((1.8726+sqrt(3.5066+4.6757*(35.957+tdp1)))/2.3378))/1000.0;
  sp->tdb=tdp1-1.0;
  dummy1=(-10.);

  do {
    dummy1=(-dummy1/10.0);
    do {
      sp->tdb=sp->tdb+dummy1;
    } while ( ((sp->pp/psat(sp->tdb)-rh1)*dummy1)>=0.0);
  } while (fabs(dummy1)>=1.0e-6);

  sp->ppsat=psat(sp->tdb);
  sp->w=humrat(sp->pp, sp->airp);
  sp->h=enthal(sp->tdb, sp->w);
  sp->spvol=spvol(sp->tdb, sp->w, sp->airp);
  sp->twb=wetblb(sp->tdb, sp->w, sp->h, sp->airp);
  sp->mu=degsat(rh1, sp->ppsat, sp->airp);
  sp->dens=(1.0/sp->spvol)*(1.0+sp->w);
  return;
}

/*-----------------------------------------------------------------------*/
/* Inputs are : relative humidity (decimal) and enthalpy [kg/kJ].        */
/* All 11 properties are calculated and stored  in sp structure.         */
/*-----------------------------------------------------------------------*/
void rh_h(double rh1, double h1, t_psyprops *sp) {
  double dummy1, dummy2;

  sp->rh=rh1;
  sp->h=h1;
  sp->tdb=h1/1.1+1.0;
  dummy2=10.0;

  do {
    dummy2=(-dummy2/10.0);
    do {
      sp->tdb=sp->tdb+dummy2;
      sp->ppsat=psat(sp->tdb);
      sp->pp=rh1*sp->ppsat;
      sp->w=humrat(sp->pp, sp->airp);
      dummy1=enthal(sp->tdb, sp->w);
    } while (((dummy1-h1)*dummy2)<=0.0);
  } while (fabs(dummy2)>=.01);

  sp->ppsat=psat(sp->tdb);
  sp->pp=rh1*sp->ppsat;
  sp->tdp=dewpt(sp->tdb, sp->pp);
  sp->w=humrat(sp->pp, sp->airp);
  sp->twb=wetblb(sp->tdb, sp->w, h1, sp->airp);
  sp->mu=degsat(rh1, sp->ppsat, sp->airp);
  sp->spvol=spvol(sp->tdb, sp->w, sp->airp);
  sp->dens=(1.0/sp->spvol)*(1.0+sp->w);
  return;
}

/*-----------------------------------------------------------------------*/
/* Inputs are : humidity ratio and wet-bulb temperature [C] .            */
/* All 11 properties are calculated and stored in sp structure.          */
/*-----------------------------------------------------------------------*/
void w_wb(double w1, double twb1, t_psyprops *sp) {
  double dummy1, dummy2, wstar;

  sp->ppsat=psat(twb1); /*Note! It is different from original Pascal */

  wstar=humrat(sp->ppsat, sp->airp);

  sp->w=w1;
  sp->twb=twb1;

  sp->tdb=twb1-1.0;
  dummy2=(-10.0);

  do {
    dummy2=(-dummy2/10.0);
    do {
      sp->tdb=sp->tdb+dummy2;
      sp->h=enthal(sp->tdb, w1);
      dummy1=wetblb(sp->tdb, w1, sp->h, sp->airp);
    } while ((dummy2*(dummy1-twb1)<= 0.0) || (wstar < w1));
  } while (fabs(dummy2)>=0.01); /* differs from original Pascal*/

  sp->ppsat=psat(sp->tdb);
  sp->spvol=spvol(sp->tdb, w1, sp->airp);
  sp->dens=(1.0/sp->spvol)*(1.0+w1);
  sp->pp=(sp->airp*w1)/(w1+.62198);
  sp->rh=sp->pp/sp->ppsat;
  sp->mu=degsat(sp->rh, sp->ppsat, sp->airp);
  sp->tdp=dewpt(sp->tdb, sp->pp);
  return;
}

/*-----------------------------------------------------------------------*/
/* Inputs are : humidity ratio and enthalpy [kg/kJ].                     */
/* All 11 properties are calculated and stored in sp structure.          */
/*-----------------------------------------------------------------------*/
void w_h(double w1, double h1, t_psyprops *sp) {
  sp->w=w1;
  sp->h=h1;
  sp->pp=(sp->airp*w1)/(w1+.62198);
  sp->tdb=(h1-2501.0*w1)/(1.006+1.805*w1);
  sp->tdp=dewpt(sp->tdb, sp->pp);
  sp->ppsat=psat(sp->tdb);
  sp->spvol=spvol(sp->tdb, w1, sp->airp);
  sp->dens=(1.0/sp->spvol)*(1.0+w1);
  sp->twb=wetblb(sp->tdb, w1, h1, sp->airp);
  sp->rh=sp->pp/sp->ppsat;
  sp->mu=degsat(sp->rh, sp->ppsat, sp->airp);
  return;
}

/*------------------------------------------------------------------------*/
/* Inputs are: wet-bulb and dew-point temperatures [C].                   */
/* All 11 properties are calculated and stored  in sp structure.          */
/*------------------------------------------------------------------------*/
void wb_dp(double twb1, double tdp1, t_psyprops *sp) {
  double dummy1, dummy2;

  sp->tdp=tdp1;
  if (tdp1<=0.0)
    sp->pp=(exp((-7.0322+sqrt(49.452-1.48*(-60.45-tdp1)))/.74))/1000.0;
  else
    sp->pp=(exp((1.8726+sqrt(3.5066+4.7657*(35.957+tdp1)))/2.3378))/1000.0;

  sp->w=humrat(sp->pp, sp->airp);
  sp->tdb=twb1-1.0;
  dummy2=(-10.);

  do {
    dummy2=(-dummy2/10.0);
    do {
      sp->tdb=sp->tdb+dummy2;
      sp->h=enthal(sp->tdb, sp->w);
      dummy1=wetblb(sp->tdb, sp->w, sp->h, sp->airp);
    } while (dummy2*(dummy1-twb1)<=0.0);
  } while (fabs(dummy2)>=0.01);

  sp->ppsat=psat(sp->tdb);
  sp->rh=sp->pp/sp->ppsat;
  sp->mu=degsat(sp->tdb, sp->ppsat, sp->airp);
  sp->spvol=spvol(sp->tdb, sp->w, sp->airp);
  sp->dens=(1.0/sp->spvol)*(1.0+sp->w);
  sp->h=enthal(sp->tdb, sp->w);
  sp->twb=wetblb(sp->tdb, sp->w, sp->h, sp->airp);
  return;
}

/*----------------------------------------------------------------------*/
/* Inputs are : wet-blub temperature [C] and enthalpy [kg/kJ].          */
/* All 11 properties are calculated and stored in sp structure.         */
/*----------------------------------------------------------------------*/
void wb_h(double twb1, double h1, t_psyprops *sp) {
  double dummy1, dummy2;

  sp->h=h1;
  sp->twb=twb1;

  sp->tdb=twb1-1.0;
  dummy2=(-10.0);

  do {
    dummy2=(-dummy2/10.0);
    do {
      sp->tdb=sp->tdb+dummy2;
      sp->w=(h1-1.006*sp->tdb)/(2501.0+1.805*sp->tdb);
      dummy1=wetblb(sp->tdb, sp->w, h1, sp->airp);
    } while (dummy2*(dummy1-twb1)<=0.0);
  } while (fabs(dummy2)>=0.01);

  sp->ppsat=psat(sp->tdb);
  sp->w=(h1-1.006*sp->tdb)/(2501.0+1.805*sp->tdb);
  sp->pp=(sp->airp * sp->w)/(sp->w+.62198);
  sp->rh=sp->pp/sp->ppsat;
  sp->mu=degsat(sp->tdb, sp->ppsat, sp->airp);
  sp->tdp=dewpt(sp->tdb, sp->pp);
  sp->spvol=spvol(sp->tdb, sp->w, sp->airp);
  sp->dens=(1.0/sp->spvol)*(1.0+sp->w);
  return;
}

/*----------------------------------------------------------------------*/
/*Inputs are : dewpoint temperature [C], and enthalpy  [kg/kJ].         */
/* All 11 properties are calculated and stored in sp structure.         */
/*----------------------------------------------------------------------*/
void dp_h(double tdp1, double h1, t_psyprops *sp) {
  sp->tdp=tdp1;
  sp->h=h1;

  if (tdp1<=0.0)
    sp->pp=(exp((-7.0322+sqrt(49.452-1.48*(-60.45-tdp1)))/.74))/1000.0;
  else
    sp->pp=(exp((1.8726+sqrt(3.5066+4.6757*(35.957+tdp1)))/2.3378))/1000.0;

  sp->w=.62198*sp->pp/(sp->airp - sp->pp);
  sp->tdb=(h1-2501.0*sp->w)/(1.006+1.805*sp->w);
  if (sp->tdb<-40.0)
    sp->tdb=(-40.0);
  sp->ppsat=psat(sp->tdb);
  sp->rh=sp->pp/sp->ppsat;
  sp->spvol=spvol(sp->tdb, sp->w, sp->airp);
  sp->dens=(1.0/sp->spvol)*(1.0+sp->w);
  sp->twb=wetblb(sp->tdb, sp->w, h1, sp->airp);
  sp->mu=degsat(sp->rh, sp->ppsat, sp->airp);
  return;
}

/*-----------------------------------------------------------------------*/

