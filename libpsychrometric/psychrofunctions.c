#include <math.h>
#include "psychrometric.h"

#define C14     6.54
#define C15     14.526
#define C16     0.7389
#define C17     0.09486
#define C18     0.4569

double func_dprh2db(double dp, double rh) {
  double e, es;
  e = func_t2es(dp);
  es = e/rh;
  return (func_e2dp(es));
}

double func_wh2db(double w, double h) {
  return ((h-2501*w)/(1.006+1.805*w));
}

double func_wbrh2db(double wb, double rh, double p) {
  double w, w1, db, db1, db2, twb;
  int i;
  db = wb;
  w = func_dbrh2w(db, rh, p);
  w1 = w;
  while (w1<=w) {
    w1 = w;
    w = func_dbrh2w(db, rh, p);
    db = func_wbw2db(wb, w, p);
  }
  i=0;
  db1 = wb;
  db2 = db;
  db = db1+(db2-db1)/2;
  while (((db2-db1)>0.001)&&(i<1000)) {
    twb = func_dbrh2wb(db, rh, p);
    if (twb>wb) { /* overestimate the dry-bulb temperature */
      db2 = db; /* the wb should be between wb1,twb */
      db = db1+(db2-db1)/2;
    } else { /* underestimate the wb */
      db1 = db;
      db = db1+(db2-db1)/2;
    }
    i ++;
  }
  return db;
}

double func_wbw2db(double wb, double w, double p) {
  double wstar, es, t1, t2;
  es = func_t2es(wb);
  wstar = func_e2w(es, p);
  t1 = (2501-2.381*wb)*wstar-2501*w+4.186*wb*w+wb;
  t2 = 1.805*w+1;
  return (t1/t2);
}

double func_wbh2db(double wb, double h, double p) {
  double w, db, db2, diff, h2, h1;
  db2 = wb;
  diff = 100;
  h1=h2=h;
  while (diff>0.001) {
    w = func_dbh2w(db2, h);
    db = db2;
    db2 = func_wbw2db(wb, w, p);
    h1 = h2;
    h2 = func_dbw2h(db2, w);
    if ((db2>db)&&(h2>h1)) /* overestimate */
      db2 = db2-(db2-db)/2;
    else
      /* under estimate */
      db2 = db2+(db2-db)/2;
    diff = fabs(h2-h);
  }
  return db;
}

double func_hrh2db(double h, double rh, double p) {
  double w, e, db, db1, db2, rh2;
  int gotLowBound, gotHiBound, i;
  gotLowBound=gotHiBound=0;
  db = h/2;
  db1 = db2 = db;
  rh2=rh;
  while (!gotLowBound||!gotHiBound) {
    e = func_dbrh2e(db, rh);
    w = func_e2w(e, p);
    db = func_wh2db(w, h);
    rh2 = e/func_t2es(db);
    if (rh2<rh) { /* over estimate db */
      if (!gotHiBound) {
        db2 = db;
        gotHiBound = 1;
      }
    } else /* under estimate */
    if (!gotLowBound) {
      db1 = db;
      gotLowBound = 1;
    }
  }
  db = db1+(db2-db1)/2;
  i=0;
  while (((db2-db1)>0.001)&&(i<1000)) {
    w = func_dbh2w(db, h);
    e = func_w2e(w, p);
    rh2 = e/func_t2es(db);
    if (rh2<rh) { /* overestimate the dry-bulb temperature */
      db2 = db; /* the wb should be between wb1,twb */
      db = db1+(db2-db1)/2;
    } else { /* underestimate the wb */
      db1 = db;
      db = db1+(db2-db1)/2;
    }
    i++;
  }
  return db;
}

double func_e2dp(double e) {
  double af, t1, t2, s1, s2, s3;
  af = log(e);
  s1 = pow(af, (double)2);
  s2 = pow(af, (double)3);
  s3 = pow(e, 0.1984);
  t1 = C14+C15*af+C16*s1+C17*s2+C18*s3;
  t2 = 6.09+12.608*af+0.4959*s1;
  return (t1>0 ? t1 : t2);
}

double func_dbdp2wb(double db, double dp, double p) {
  double e, w, wb, wb1, wb2, twb, tw;
  int i;
  i=0;
  e = func_t2es(dp);
  w = func_e2w(e, p);
  wb1 = dp;
  wb2 = db;
  twb = wb1+(wb2-wb1)/2;
  while (((wb2-wb1)>0.001)&&(i<1000)) {
    tw = func_dbwb2w(db, twb, p);
    if (tw>w) { /* overestimate the wetblb temperature */
      wb2 = twb; /* the wb should be between wb1,twb */
      twb = wb1+(wb2-wb1)/2;
    } else { /* underestimate the wb */
      wb1 = twb;
      twb = wb1+(wb2-wb1)/2;
    }
    i ++;
  }
  wb = twb;
  return (wb);
}

double func_dbw2wb(double db, double w, double p) {
  double e, dp;
  e = func_w2e(w, p);
  dp = func_e2dp(e);
  return (func_dbdp2wb(db, dp, p));
}

double func_dbrh2wb(double db, double rh, double p) {
  double e, w;
  e = func_t2es(db)*rh;
  w = func_e2w(e, p);
  return (func_dbw2wb(db, w, p));
}

double func_ees2rh(double e, double es) {
  return (e/es);
}

double func_e2w(double e, double p) {
  return (0.62198*e/(p-e));
}

double func_dbwb2w(double db, double wb, double p) {
  double t1, t2, wstar, estar;
  estar = func_t2es(wb);
  wstar = func_e2w(estar, p);
  t1 = (2501-2.381*wb)*wstar-(db-wb);
  t2 = 2501+1.805*db-4.186*wb;
  return (t1/t2);
}

double func_dbh2w(double db, double h) {
  return ((h-1.006*db)/(2501+1.805*db));
}

double func_dbrh2w(double db, double rh, double p) {
  return (func_e2w(func_t2es(db)*rh, p));
}

double func_dbw2h(double db, double w) {
  return (1.006*db+w*(2501.0+1.805*db));
}

double func_dbrh2e(double db, double rh) {
  return (func_t2es(db)*rh);
}

double func_w2e(double w, double p) {
  return (p*w/(0.62198+w));
}

double func_dp2e(double dp) {
  return (func_t2es(dp));
}

double func_rhws2ds(double rh, double ws) {
  return (rh/(1+(1-rh)*ws/0.62198));
}

double func_t2es(double t) {
  t += 273.15;
  if (t>273.15)
    return (exp(-5800.2206/t+1.3914993-.048640239*t+(.41764768e-4)*pow(t,
        2.0)-(.14452093e-7)*pow(t, 3.0)+6.5459673*log(t))/1000.0);
  else
    return (exp(-5674.5359/t+6.3925247-(.9677843e-2)*t+(.62215701e-6)*pow(
        t, 2.0)+(.20747825e-8)*pow(t, 3.0)-(.9484024e-12)*pow(t, 4.0)
        +4.1635019*log(t))/1000.0);
}

double func_dbw2v(double db, double w, double p) {
  db += 273.16;
  return ((287.055*db*(1+1.6078*w))/(p*1000.0));
}

double func_wv2da(double w, double v) {
  return ((1+w)/v);
}

double func_es2ws(double es, double p) {
  return (func_e2w(es, p));
}

double func_w2q(double w) {
  return (w/(1+w));
}

double func_wv2X(double w, double v) {
  return (w/v);
}
