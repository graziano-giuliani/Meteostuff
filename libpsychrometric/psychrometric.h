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

#ifndef PSYCHROMETRIC_H_
#define PSYCHROMETRIC_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef struct __psyprops {
  double tdb;   /* dry-bulb temperature [C].                   */
  double tdp;   /* dewpoint temperature [C].                   */
  double h;   /* enthalpy of air/water vapor mixture [kg/kJ].*/
  double mu;   /* degree of saturation (decimal) of moist air.*/
  double ppsat;   /* partial pressure of water vapor at          */
      /* saturation [kPa].                           */
  double pp;   /* water vapor pressure [kPa].                 */
  double rh;   /* relative humidity of moist air (decimal).   */
  double spvol;   /* specific volume of air/water vapor          */
      /* mixture [m^3/kg].                           */
  double twb;   /* wet-bulb temperature [C].                   */
  double w;   /* humidity ratio of air/water vapor mixture.  */
  double dens;   /* air density [kg/m^3].                       */
  double airp;   /* air pressure   ( 101.325 [kPa])             */
} t_psyprops;

/*----------------------------------------------------------------------*/
/* Calculates the dew_point temperature [C].                            */
/* Based on formula in the ASHRAE Handbook of Fundamentals.             */
/* Applise for temperatures between -60 and +70 [C].                    */
/* Inputs are dry-bulb temperature [C] and water vapor pressure [kPa].  */
/*----------------------------------------------------------------------*/
double dewpt(double t, double kpa);

/*----------------------------------------------------------------------*/
/* Calculates partial pressure of water vapor at saturation [kPa] .     */
/* Based on the formulae provided by ASHRAE Handbook of Fundamentals.   */
/* Applies for temperatures between -40 and +120 [C].                   */
/* input: tdb---dry_bulb temperature  [C]                               */
/* output: function=patial saturation pressure  [kPa]                   */
/*----------------------------------------------------------------------*/
double psat(double tdb);

/*----------------------------------------------------------------------*/
/* calculates the enthalpy of air/water vapor mixture [kJ/kg]           */
/* Based on furmula in the ASHRAE Handbook of Fundamentals.             */
/* Inputs are dry-bulb temperature [C] and humidity ratio (unitless).   */
/*----------------------------------------------------------------------*/
double enthal(double tdb, double w);

/*----------------------------------------------------------------------*/
/* spvol() gets the specific volume of air/water vapor mixture [m^3/kg] */
/* input: tc--- dry_bulb temperature [C]                                */
/*        w---- humidity ratio                                          */
/*        pa---- air pressure  [kPa]                                    */
/* output: function = specific volume                                   */
/*----------------------------------------------------------------------*/
double spvol(double tc, double w, double pa);

/*----------------------------------------------------------------------*/
/* calculates the humidity ratio of air/water vapor mixture.            */
/* Based on formula in the ASHRAE Handbook of Fundamentals.             */
/* Inputs is water vapor pressure [kPa]. Assume aipr=101.325 [kPa]      */
/* output : function = humidity ratio                                   */
/*----------------------------------------------------------------------*/
double humrat(double p, double pa);

/*----------------------------------------------------------------------*/
/*  wetblb()                                                            */
/* It solves for the wet-bulb temperature iteratively,using enthalpy.   */
/* Based on formula in the ASHRAE Gandbook of Fundamentals.             */
/* Inputs:dry-blub temperature [C]. humidity ratio (unitless)           */
/* enthalpy [kJ/kg]. and air pressure [kPa].                            */
/* Wet-bulb temperature is calculated to the nearest .01 C.             */
/*----------------------------------------------------------------------*/
double wetblb(double tstar, double w, double hh, double pa);

/*----------------------------------------------------------------------*/
/* Calculates relative humidity (decimal) of air/water vapor mixture.   */
/* based on formula in the ASHRAE Handbook of Fundamentals .            */
/* Inputs are degree of saturation---mu (decimal),                      */
/*                      and water vapor pressure [kPa].                 */
/*----------------------------------------------------------------------*/
double relhum(double mu, double pws, double pa);

/*----------------------------------------------------------------------*/
/* Calculates the degree of sat. (decimal) of air/water vapor mixture . */
/* Based on formula in the ASHRAE Handbook of Fundamentals.             */
/* Inputs are relative humidity (decimal), water vapor pressure [kPa].  */
/*----------------------------------------------------------------------*/
double degsat(double rh, double pws, double pa);

/*----------------------------------------------------------------------*/
/* Inputs are : dry-bulb temperature [C] and relative humidity (decimal)*/
/* All 11 properties are calculated and stored in sp sturcture.         */
/*----------------------------------------------------------------------*/
void db_rh(double tdb, double rh, t_psyprops *sp);

/*----------------------------------------------------------------------*/
/* Inputs are dry-bulb temperature [C] and humidity ratio.              */
/* All 11 properties are calculated and stored in sp structure.         */
/*----------------------------------------------------------------------*/
void db_w(double tdb1, double w1, t_psyprops *sp);

/*----------------------------------------------------------------------*/
/* Inputs are :dry-bulb temperature[C] and wet-bulb temperature[C].     */
/* All 11 properties are calculated and stored in sp structure.         */
/*----------------------------------------------------------------------*/
void db_wb(double tdb1, double twb1, t_psyprops *sp);

/*----------------------------------------------------------------------*/
/* Inputs are dry-bulb and dew point temperature [C].                   */
/* All 11 properties are calculated and stored in sp sturcture.         */
/*----------------------------------------------------------------------*/
void db_dp(double tdb1, double tdp1, t_psyprops *sp);

/*----------------------------------------------------------------------*/
/* Inputs are dry-bulb temperature[C] and enthalpy [kg/kJ].             */
/* All 11 properties are calculated and stored in sp structure .        */
/*----------------------------------------------------------------------*/
void db_h(double tdb1, double h1, t_psyprops *sp);

/*----------------------------------------------------------------------*/
/* Inputs are: relative humidity (decimal) and humidity ratio.          */
/* All 11 properties are calculated and stored in sp structure.         */
/*----------------------------------------------------------------------*/
void rh_w(double rh1, double w1, t_psyprops *sp);

/*----------------------------------------------------------------------*/
/* Inputs are: relative humidity (decimal) and wet-bulb temperature [C].*/
/* All 11 properties are calculated and stored in sp structure.         */
/*----------------------------------------------------------------------*/
void rh_wb(double rh1, double twb1, t_psyprops *sp);

/*----------------------------------------------------------------------*/
/* Inputs are: relative humidity (decimal) and dew-point temperature [C]*/
/* All 11 properties are calculated and stored in sp structure .        */
/*----------------------------------------------------------------------*/
void rh_dp(double rh1, double tdp1, t_psyprops *sp);

/*----------------------------------------------------------------------*/
/* Inputs are : relative humidity (decimal) and enthalpy [kg/kJ].       */
/* All 11 properties are calculated and stored  in sp structure.        */
/*----------------------------------------------------------------------*/
void rh_h(double rh1, double h1, t_psyprops *sp);

/*----------------------------------------------------------------------*/
/* Inputs are : humidity ratio and wet-bulb temperature [C] .           */
/* All 11 properties are calculated and stored in sp structure.         */
/*----------------------------------------------------------------------*/
void w_wb(double w1, double twb1, t_psyprops *sp);

/*----------------------------------------------------------------------*/
/* Inputs are : humidity ratio and enthalpy [kg/kJ].                    */
/* All 11 properties are calculated and stored in sp structure.         */
/*----------------------------------------------------------------------*/
void w_h(double w1, double h1, t_psyprops *sp);

/*----------------------------------------------------------------------*/
/* Inputs are: wet-bulb and dew-point temperatures [C].                 */
/* All 11 properties are calculated and stored  in sp structure.        */
/*----------------------------------------------------------------------*/
void wb_dp(double twb1, double tdp1, t_psyprops *sp);

/*----------------------------------------------------------------------*/
/* Inputs are : wet-blub temperature [C] and enthalpy [kg/kJ].          */
/* All 11 properties are calculated and stored in sp structure.         */
/*----------------------------------------------------------------------*/
void wb_h(double twb1, double h1, t_psyprops *sp);

/*----------------------------------------------------------------------*/
/*Inputs are : dewpoint temperature [C], and enthalpy  [kg/kJ].         */
/* All 11 properties are calculated and stored in sp structure.         */
/*----------------------------------------------------------------------*/
void dp_h(double tdp1, double h1, t_psyprops *sp);

/*----------------------------------------------------------------------*/
/*Conversion functions from one parameter to another.                   */
/*----------------------------------------------------------------------*/
double func_dprh2db(double dp, double rh);
double func_wh2db(double w, double h);
double func_wbrh2db(double wb, double rh, double p);
double func_wbw2db(double wb, double w, double p);
double func_wbh2db(double wb, double h, double p);
double func_hrh2db(double h, double rh, double p);
double func_e2dp(double e);
double func_dbdp2wb(double db, double dp, double p);
double func_dbw2wb(double db, double w, double p);
double func_dbrh2wb(double db, double rh, double p);
double func_ees2rh(double e, double es);
double func_e2w(double e, double p);
double func_dbwb2w(double db, double wb, double p);
double func_dbh2w(double db, double h);
double func_dbrh2w(double db, double rh, double p);
double func_dbw2h(double db, double w);
double func_dbrh2e(double db, double rh);
double func_w2e(double w, double p);
double func_dp2e(double dp);
double func_rhws2ds(double rh, double ws);
double func_t2es(double t);
double func_dbw2v(double db, double w, double p);
double func_wv2da(double w, double v);
double func_es2ws(double es, double p);
double func_w2q(double w);
double func_wv2X(double w, double v);

#ifdef __cplusplus  
}
#endif

#endif /*PSYCHROMETRIC_H_*/
