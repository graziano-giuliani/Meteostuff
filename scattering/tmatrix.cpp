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
#include <gauss.h>
#include <bessel.h>
#include <tmatrix.h>
#include <invert.h>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include "boost/format.hpp"

#ifndef TESTME

namespace cetemps {

  void Tmatrix::Setup(scattering_particle &p, radiation &rad,
                      double ddelt, int base_gaussian_divisions)
  {
    t_calculus_space t;

    nmax = p.nmax(rad.wavelenght);
    ngmax = base_gaussian_divisions*nmax;
    nmax2 = 2*nmax;
    ngmax2 = 2*ngmax;
    resize_problem(t);

    ncheck = 0;
    if (p.shape == particle_spheroid ||
        p.shape == particle_cylinder) ncheck = 1;
    if (p.shape == particle_chebyshev)
    {
      chebyshev *c = (chebyshev *) &p;
       if ((c->get_polynomial_degree( )%2) == 0) ncheck = 1;
    }

    double qext1 = 0.0;
    double qsca1 = 0.0;
    double dsca = 0.0;
    double dext = 0.0;
    bool done_convergence = false;

    do
    {
      Constant(p, t);
      Vary(p, rad, t);
      Tmatr0(t);

      qext = 0.0;
      qsca = 0.0;
      for (int i = 0; i < nmax; i ++)
      {
        double fact = (double) (2*(i+1)+1);
        qext += fact * real(t.ct_t(i,i) + t.ct_t(nmax+i,nmax+i));
        qsca += fact * (norm(t.ct_t(i,i)) + norm(t.ct_t(nmax+i,nmax+i)));
      }
      dsca = fabs((qsca1-qsca)/qsca);
      dext = fabs((qext1-qext)/qext);
      done_convergence = (dsca <= ddelt && dext <= ddelt);
      if (! done_convergence)
      {
        nmax ++;
        ngmax = base_gaussian_divisions*nmax;
        nmax2 = 2*nmax;
        ngmax2 = 2*ngmax;
        resize_problem(t);
      }
      qext1 = qext;
      qsca1 = qsca;
    } while (! done_convergence);

    ngmax ++;
    ngmax2 = 2*ngmax;
    resize_problem(t);

    done_convergence = false;
    do
    {
      Constant(p, t);
      Vary(p, rad, t);
      Tmatr0(t);
      qext = 0.0;
      qsca = 0.0;
      for (int i = 0; i < nmax; i ++)
      {
        double fact = (double) (2*(i+1)+1);
        qext += fact * real(t.ct_t(i,i) + t.ct_t(nmax+i,nmax+i));
        qsca += fact * (norm(t.ct_t(i,i)) + norm(t.ct_t(nmax+i,nmax+i)));
      }
      dsca = fabs((qsca1-qsca)/qsca);
      dext = fabs((qext1-qext)/qext);
      done_convergence = (dsca <= ddelt && dext <= ddelt);
      if (! done_convergence)
      {
        ngmax ++;
        ngmax2 = 2*ngmax;
        resize_problem(t);
      }
      qext1 = qext;
      qsca1 = qsca;
    } while (! done_convergence);

    qsca = 0.0;
    qext = 0.0;

    tmat_t11.resize(extents[nmax+1][nmax][nmax]);
    tmat_t12.resize(extents[nmax+1][nmax][nmax]);
    tmat_t21.resize(extents[nmax+1][nmax][nmax]);
    tmat_t22.resize(extents[nmax+1][nmax][nmax]);

    for (int n2 = 0; n2 < nmax; n2 ++)
    {
      int nn2 = nmax+n2;
      for (int n1 = 0; n1 < nmax; n1 ++)
      {
        int nn1 = nmax+n1;
        dcplex zz1 = t.ct_t(n1,n2);
        dcplex zz2 = t.ct_t(n1,nn2);
        dcplex zz3 = t.ct_t(nn1,n2);
        dcplex zz4 = t.ct_t(nn1,nn2);
        tmat_t11[0][n1][n2] = zz1;
        tmat_t12[0][n1][n2] = zz2;
        tmat_t21[0][n1][n2] = zz3;
        tmat_t22[0][n1][n2] = zz4;
        qsca += norm(zz1) + norm(zz2) + norm(zz3) + norm(zz4);
      }
    }
    for (int i = 0; i < 2*nmax; i ++)
      qext += real(t.ct_t(i,i));

    for (int m = 0; m < nmax; m ++)
    {
      Tmatr(m, t);
      int nm = nmax-m;
      qsca1 = 0.0;
      for (int n2 = 0; n2 < nm; n2 ++)
      {
        int nn2 = n2+m;
        int n22 = n2+nm;
        for (int n1 = 0; n1 < nm; n1 ++)
        {
          int nn1 = n1+m;
          int n11 = n1+nm;
          dcplex zz1 = t.ct_t(n1,n2);
          dcplex zz2 = t.ct_t(n1,n22);
          dcplex zz3 = t.ct_t(n11,n2);
          dcplex zz4 = t.ct_t(n11,n22);
          tmat_t11[m+1][nn1][nn2] = zz1;
          tmat_t12[m+1][nn1][nn2] = zz2;
          tmat_t21[m+1][nn1][nn2] = zz3;
          tmat_t22[m+1][nn1][nn2] = zz4;
          qsca1 += 2.0 * (norm(zz1)+norm(zz2)+norm(zz3)+norm(zz4));
        }
      }
      qext1 = 0.0;
      for (int n = 0; n < 2*nm; n ++)
        qext1 += 2.0 * real(t.ct_t(n,n));
      qext += qext1;
      qsca += qsca1;
    }
  }

  int Tmatrix::resize_problem(t_calculus_space &t)
  {
    t.x.resize(2*ngmax);
    t.w.resize(2*ngmax);
    t.an.resize(nmax);
    t.s.resize(2*ngmax);
    t.ss.resize(2*ngmax);
    t.r.resize(2*ngmax);
    t.dr.resize(2*ngmax);
    t.ddr.resize(2*ngmax);
    t.cdr.resize(2*ngmax);
    t.ann.resize(nmax, nmax);
    t.ct_t.resize(2*nmax, 2*nmax);
    t.t99_11.resize(nmax, nmax);
    t.t99_12.resize(nmax, nmax);
    t.t99_21.resize(nmax, nmax);
    t.t99_22.resize(nmax, nmax);
    t.t99_g11.resize(nmax, nmax);
    t.t99_g12.resize(nmax, nmax);
    t.t99_g21.resize(nmax, nmax);
    t.t99_g22.resize(nmax, nmax);
    t.ctt_q.resize(2*nmax, 2*nmax);
    t.ctt_rgq.resize(2*nmax, 2*nmax);
    t.bessely.resize(2*ngmax, nmax);
    t.besseldy.resize(2*ngmax, nmax);
    t.besselj.resize(2*ngmax, nmax);
    t.besseldj.resize(2*ngmax, nmax);
    t.cbesselj.resize(2*ngmax, nmax);
    t.cbesseldj.resize(2*ngmax, nmax);
    return 0;
  }

  // Calculation of Wigner D functios by means of the upward recurrence 
  // relation (see D.A. Varshalovich, A.N. Moskalev, and V.K. Khersonsky,
  // "Quantum Theory of Angular Momentum", World Scientific, Singapore, 1998)

  void Tmatrix::WignerD(double x, int m,
                        Vdouble &dv1, Vdouble &dv2)
  {
    int umax = dv1.size( );
    // BOOST_UBLAS_CHECK(umax == dv2.size( ), external_logic());
    zero_vector<double> zero_vec(umax);

    dv1 = zero_vec;
    dv2 = zero_vec;

    double qs = sqrt(1.0 - x * x);
    double qs1 = 1.0 / qs;

    if (m == 0)
    {
      double d1 = 1.0;
      double d2 = x;
      for (int n = 1; n <= umax; n++)
      {
        double qn = (double) n;
        double qn1 = (double) (n + 1);
        double qn2 = (double) (2*n + 1);
        double d3 = (qn2 * x * d2 - qn * d1) / qn1;
        double der = qs1 * (qn1 * qn / qn2) * (-d1 + d3);
        dv1(n-1) = d2;
        dv2(n-1) = der;
        d1 = d2;
        d2 = d3;
      }
      return;
    }
    double qmm = (double) (m * m);
    double a = 1.0;
    for (int i = 1; i <= m; i++)
      a = a * sqrt((double) (i*2 - 1) / (double) (i*2)) * qs;
    double d1 = 0.0;
    double d2 = a;
    for (int n = m; n <= umax; n++)
    {
      double qn = (double) n;
      double qn2 = (double) (2*n + 1);
      double qn1 = (double) (n + 1);
      double qnm = sqrt(qn * qn - qmm);
      double qnm1 = sqrt(qn1 * qn1 - qmm);
      double d3 = (qn2 * x * d2 - qnm * d1) / qnm1;
      double der = qs1 * (-qn1 * qnm * d1 + qn * qnm1 * d3) / qn2;
      dv1(n-1) = d2;
      dv2(n-1) = der;
      d1 = d2;
      d2 = d3;
    }
    return;
  }

  void Tmatrix::WignerDAm(double x, int m,
                          Vdouble &dv1, Vdouble &dv2)
  {
    int umax = dv1.size( );
    // BOOST_UBLAS_CHECK(umax == dv2.size( ), external_logic());
    zero_vector<double> zero_vec(umax);

    dv1 = zero_vec;
    dv2 = zero_vec;

    double qs = sqrt(1.0 - x * x);
    double qs1 = 1.0 / qs;
    double dsi = qs1;

    double dx = fabs(x);
    if (fabs(1.0-dx) <= 1e-10)
    {
      if (m != 1) return;
      for (int i = 0; i < umax; i ++)
      {
        double dn = (double) ((i+1)*(i+2));
        dn = 0.5*sqrt(dn);
        if (x < 0.0) dn = dn*pow(-1,i+2);
        dv1(i) = dn;
        if (x < 0.0) dn=-dn;
        dv2(i) = dn;
      }
      return;
    }

    if (m == 0)
    {
      double d1 = 1.0;
      double d2 = x;
      for (int n = 1; n <= umax; n++)
      {
        double qn = (double) n;
        double qn1 = (double) (n + 1);
        double qn2 = (double) (2*n + 1);
        double d3 = (qn2 * x * d2 - qn * d1) / qn1;
        double der = qs1 * (qn1 * qn / qn2) * (-d1 + d3);
        dv1(n-1) = d2*dsi;
        dv2(n-1) = der;
        d1 = d2;
        d2 = d3;
      }
      return;
    }
    double qmm = (double) (m * m);
    double a = 1.0;
    for (int i = 1; i <= m; i++)
      a = a * sqrt((double) (i*2 - 1) / (double) (i*2)) * qs;
    double d1 = 0.0;
    double d2 = a;
    for (int n = m; n <= umax; n++)
    {
      double qn = (double) n;
      double qn2 = (double) (2*n + 1);
      double qn1 = (double) (n + 1);
      double qnm = sqrt(qn * qn - qmm);
      double qnm1 = sqrt(qn1 * qn1 - qmm);
      double d3 = (qn2 * x * d2 - qnm * d1) / qnm1;
      double der = qs1 * (-qn1 * qnm * d1 + qn * qnm1 * d3) / qn2;
      dv1(n-1) = d2*dsi;
      dv2(n-1) = der;
      d1 = d2;
      d2 = d3;
    }
    return;
  }

  void Tmatrix::Constant(scattering_particle &p, t_calculus_space &t)
  {
    Vdouble dd(ngmax2);
    for (int n = 0; n < nmax; n++)
    {
      int nn = (n+1)*(n+2);
      t.an(n) = (double) nn;
      double d = sqrt((double) ((2*n) + 3) / (double) nn);
      dd(n) = d;
      for (int n1 = 0; n1 <= n; n1++)
      {
        double ddd = d * dd(n1) * 0.5;
        t.ann(n, n1) = ddd;
        t.ann(n1, n) = ddd;
      }
    }
    if (p.shape != particle_cylinder)
    {
      gauss(ngmax2, t.x, t.w);
    }
    else // a cylinder
    {
      Vdouble x1(ngmax2);
      Vdouble w1(ngmax2);
      Vdouble x2(ngmax2);
      Vdouble w2(ngmax2);

      cylinder *cyl = (cylinder*)&p;
      int ng1 = (int) ((double) (ngmax) / 2.);
      int ng2 = ngmax - ng1;
      gauss(ng1, x1, w1);
      gauss(ng2, x2, w2);
      double xx = -cos(atan(cyl->get_hv_dimension_ratio()));
      for (int i = 0; i < ng1; i++)
      {
        t.w(i) = (xx + 1.0) * 0.5 * w1(i);
        t.x(i) = (xx + 1.0) * 0.5 * x1(i) + (xx - 1.0) * 0.5;
      }
      for (int i = 0; i < ng2; i++)
      {
        t.w(i + ng1) = xx * -0.5 * w2(i);
        t.x(i + ng1) = xx * -0.5 * x2(i) + xx * 0.5;
      }
      for (int i = 0; i < ngmax; i++)
      {
        t.w(ngmax2 - i - 1) = t.w(i);
        t.x(ngmax2 - i - 1) = -t.x(i);
      }
    }
    for (int i = 0; i < ngmax; i++)
    {
      double y = t.x(i);
      y = 1.0 / (1.0 - y * y);
      t.ss(i) = y;
      t.ss(ngmax2 - i - 1) = y;
      y = sqrt(y);
      t.s(i) = y;
      t.s(ngmax2 - i - 1) = y;
    }
    return;
  }

  void Tmatrix::Vary(scattering_particle &p,
      radiation &rad, t_calculus_space &t)
  {
    int umax = t.r.size();
    p.rsp(t.x, t.r, t.dr);

    xp = (M_PI * 2.0) / rad.wavelenght;
    ppi = xp*xp;
    pi =  ppi*p.refractive_index;

    double ta = 0.0;

    Vdouble z(ngmax2);
    Vdcplex cz(ngmax2);
    dcplex pr = pow(p.refractive_index,-1);
    for (int i = 0; i < umax; i ++)
    {
      double vv = sqrt(t.r(i))*xp;
      ta = fmax(0.0, vv);
      t.ddr(i) = 1.0 / vv;
      t.cdr(i) = pr * t.ddr(i);
      z(i) = vv;
      cz(i) = vv * p.refractive_index;
    }

    double tb = ta * abs(p.refractive_index);
    tb = fmax(tb, ngmax2);
    int iprec1 = (int) (1.2 * sqrt(fmax(ta,ngmax2)) + 3.0);
    int iprec2 = (int) (tb+4.0*pow(tb,1.0/3.0)+1.2*sqrt(tb));
    iprec2 = iprec2-ngmax2+5;

    Vdouble a(nmax);
    Vdouble b(nmax);
    Vdcplex ca(nmax);
    Vdcplex cb(nmax);
    for (int i = 0; i < ngmax2; i ++)
    {
      matrix_row<matrix<double> > rowbj (t.besselj, i);
      matrix_row<matrix<double> > rowbdj (t.besseldj, i);
      matrix_row<matrix<double> > rowby (t.bessely, i);
      matrix_row<matrix<double> > rowbdy (t.besseldy, i);
      matrix_row<matrix<dcplex> > colcbj (t.cbesselj, i);
      matrix_row<matrix<dcplex> > colcbdj (t.cbesseldj, i);

      rjb(iprec1, z(i), a, b);
      rowbj = a;
      rowbdj = b;
      ryb(z(i), a, b);
      rowby = a;
      rowbdy = b;
      rjb(iprec2, cz(i), ca, cb);
      colcbj = ca;
      colcbdj = cb;
    }
    return;
  }

  void Tmatrix::Tmatr0(t_calculus_space &t)
  {
    int ngss = ncheck == 1 ? ngmax : ngmax2;
    double factor = ngmax2/ngss;

    const dcplex I(0.0,1.0);

    Vdouble sig(nmax2);
    sig(0) = -1.0;
    for (int i = 1; i < nmax2; i ++)
      sig(i) = -sig(i-1);

    Mdouble d1(ngmax2, nmax);
    Mdouble d2(ngmax2, nmax);
    for (int i = 0; i < ngmax; i ++)
    {
      Vdouble dv1(nmax);
      Vdouble dv2(nmax);
      int i1 = ngmax+i;
      int i2 = ngmax-i-1;
      WignerD(t.x(i1), 0, dv1, dv2);
      for (int j = 0; j < nmax; j ++)
      {
        double si = sig(j);
        double dd1 = dv1(j);
        double dd2 = dv2(j);
        d1(i1,j) = dd1;
        d2(i1,j) = dd2;
        d1(i2,j) = dd1*si;
        d2(i2,j) = -dd2*si;
      }
    }

    Vdouble rr(ngss);
    for (int i = 0; i < ngss; i ++)
      rr(i) = t.w(i)*t.r(i);

    for (int i = 0; i < nmax; i ++)
    {
      double an1;
      an1 = t.an(i);
      for (int j = 0; j < nmax; j ++)
      {
        if (ncheck == 1 && sig(i+j) > 0)
        {
          t.t99_12(i,j) = 0.0;
          t.t99_21(i,j) = 0.0;
          t.t99_g12(i,j) = 0.0;
          t.t99_g21(i,j) = 0.0;
          continue;
        }
        double an2;
        an2 = t.an(j);
        dcplex a12 = 0.0;
        dcplex a21 = 0.0;
        dcplex g12 = 0.0;
        dcplex g21 = 0.0;
        for (int k = 0; k < ngss; k ++)
        {
          double d1n1 = d1(k, i);
          double d2n1 = d2(k, i);
          double d1n2 = d1(k, j);
          double d2n2 = d2(k, j);
          double xa12 = d1n1*d2n2;
          double xa21 = d2n1*d1n2;
          double xa22 = d2n1*d2n2;

          double qj1  = t.besselj(k, i);
          double qy1  = t.bessely(k, i);
          dcplex qj2  = t.cbesselj(k, j);
          dcplex qdj2 = t.cbesseldj(k, j);
          double qdj1 = t.besseldj(k, i);
          double qdy1 = t.besseldy(k, i);

          dcplex c1   = qj2 * qj1;
          dcplex b1   = c1 + ((qj2*I)*qy1);
          dcplex c2   = qj2 * qdj1;
          dcplex b2   = c2 + ((qj2*I)*qdy1);
          double ddri = t.ddr(k);
          dcplex c3   = c1 * ddri;
          dcplex b3   = b1 * ddri;
          dcplex c4   = qdj2 * qj1;
          dcplex b4   = c4 + ((qdj2*I)*qy1);
          dcplex cdri = t.cdr(k);
          dcplex c5   = c1 * cdri;
          dcplex b5   = b1 * cdri;
          double uri  = t.dr(k);
          double rri  = rr(k);
          double f1   = rri * xa22;
          double f2   = rri * uri * an1 * xa12;
          a12  = a12 + f1 * b2 + f2 * b3;
          g12  = g12 + f1 * c2 + f2 * c3;
          f2   = rri * uri * an2 * xa21;
          a21  = a21 + f1 * b4 + f2 * b5;
          g21  = g21 + f1 * c4 + f2 * c5;
        }
        double an12 = t.ann(i, j) * factor;
        t.t99_12(i,j) = a12 * an12;
        t.t99_21(i,j) = a21 * an12;
        t.t99_g12(i,j) = g12 * an12;
        t.t99_g21(i,j) = g21 * an12;
      }
    }

    for (int i = 0; i < nmax; i ++)
    {
      int kk1 = i+nmax;
      for (int j = 0; j < nmax; j ++)
      {
        int kk2 = j+nmax;
        dcplex ta12 = -I*t.t99_12(i,j);
        dcplex tg12 = -I*t.t99_g12(i,j);
        dcplex ta21 = I*t.t99_21(i,j);
        dcplex tg21 = I*t.t99_g21(i,j);
        t.ctt_q(i, j) = pi * ta21 + ppi * ta12;
        t.ctt_rgq(i, j) = pi * tg21 + ppi * tg12;
        t.ctt_q(i, kk2) = 0.0;
        t.ctt_rgq(i, kk2) = 0.0;
        t.ctt_q(kk1, j) = 0.0;
        t.ctt_rgq(kk1, j) = 0.0;
        t.ctt_q(kk1,kk2) = pi * ta12 + ppi * ta21;
        t.ctt_rgq(kk1,kk2) = pi * tg12 + ppi * tg21;
      }
    }
    
    TT(nmax, t);

    return;
  }

  void Tmatrix::Tmatr(int m, t_calculus_space &t)
  {
    int ngss = ncheck == 1 ? ngmax : ngmax2;
    double factor = ngmax2/ngss;

    const dcplex I(0.0,1.0);

    Vdouble sig(nmax2);
    sig(0) = -1.0;
    for (int i = 1; i < nmax2; i ++)
      sig(i) = -sig(i-1);

    Mdouble d1(ngmax2, nmax);
    Mdouble d2(ngmax2, nmax);
    for (int i = 0; i < ngmax; i ++)
    {
      Vdouble dv1(nmax);
      Vdouble dv2(nmax);
      int i1, i2;
      i1 = ngmax+i;
      i2 = ngmax-i-1;
      WignerD(t.x(i1), m+1, dv1, dv2);
      for (int j = 0; j < nmax; j ++)
      {
        double si = sig(j);
        double dd1 = dv1(j);
        double dd2 = dv2(j);
        d1(i1,j) = dd1;
        d2(i1,j) = dd2;
        d1(i2,j) = dd1*si;
        d2(i2,j) = -dd2*si;
      }
    }

    Vdouble rr(ngss);
    Vdouble ds(ngss);
    Vdouble dss(ngss);
    for (int i = 0; i < ngss; i ++)
    {
      double wr = t.w(i)*t.r(i);
      ds(i) = t.s(i)*(m+1)*wr;
      dss(i) = t.ss(i)*(m+1)*(m+1);
      rr(i) = wr;
    }

    for (int i = m; i < nmax; i ++)
    {
      double an1;
      an1 = t.an(i);
      for (int j = m; j < nmax; j ++)
      {
        double an2;
        an2 = t.an(j);
        dcplex a11 = 0.0;
        dcplex a12 = 0.0;
        dcplex a21 = 0.0;
        dcplex a22 = 0.0;
        dcplex g11 = 0.0;
        dcplex g12 = 0.0;
        dcplex g21 = 0.0;
        dcplex g22 = 0.0;
        for (int k = 0; k < ngss; k ++)
        {
          double d1n1 = d1(k, i);
          double d2n1 = d2(k, i);
          double d1n2 = d1(k, j);
          double d2n2 = d2(k, j);
          double xa11 = d1n1*d1n2;
          double xa12 = d1n1*d2n2;
          double xa21 = d2n1*d1n2;
          double xa22 = d2n1*d2n2;
          double aa1  = xa12+xa21;
          double aa2  = xa11*dss(k)+xa22;

          double qj1  = t.besselj(k, i);
          double qy1  = t.bessely(k, i);
          dcplex qj2  = t.cbesselj(k, j);
          dcplex qdj2 = t.cbesseldj(k, j);
          double qdj1 = t.besseldj(k, i);
          double qdy1 = t.besseldy(k, i);

          dcplex c1   = qj2 * qj1;
          dcplex b1   = c1 + ((qj2*I)*qy1);
          dcplex c2   = qj2 * qdj1;
          dcplex b2   = c2 + ((qj2*I)*qdy1);
          double ddri = t.ddr(k);
          dcplex c3   = c1 * ddri;
          dcplex b3   = b1 * ddri;
          dcplex c4   = qdj2 * qj1;
          dcplex b4   = c4 + ((qdj2*I)*qy1);
          dcplex cdri = t.cdr(k);
          dcplex c5   = c1 * cdri;
          dcplex b5   = b1 * cdri;
          dcplex c6   = qdj2 * qdj1;
          dcplex b6   = c6 + ((qdj2*I)*qdy1);
          dcplex c7   = c4 * ddri;
          dcplex b7   = b4 * ddri;
          dcplex c8   = c2 * cdri;
          dcplex b8   = b2 * cdri;
          double uri  = t.dr(k);
          double rri  = rr(k);
          double dsi  = ds(k);
          if (ncheck == 1 && sig(i+j) < 0)
          {
            double f1   = rri * aa2;
            double f2   = rri * uri * an1 * xa12;
            a12  = a12 + f1 * b2 + f2 * b3;
            g12  = g12 + f1 * c2 + f2 * c3;
            f2   = rri * uri * an2 * xa21;
            a21  = a21 + f1 * b4 + f2 * b5;
            g21  = g21 + f1 * c4 + f2 * c5;
            continue;
          }
          double e1 = dsi * aa1;
          a11 = a11 + e1 * b1;
          g11 = g11 + e1 * c1;
          double e2 = dsi * uri * xa11;
          double e3 = e2 * an2;
          e2 = e2 * an1;
          a22 = a22 + e1 * b6 + e2 * b7 + e3 * b8;
          g22 = g22 + e1 * c6 + e2 * c7 + e3 * c8;
          if (ncheck == 1)
            continue;
          double f1   = rri * aa2;
          double f2   = rri * uri * an1 * xa12;
          a12  = a12 + f1 * b2 + f2 * b3;
          g12  = g12 + f1 * c2 + f2 * c3;
          f2   = rri * uri * an2 * xa21;
          a21  = a21 + f1 * b4 + f2 * b5;
          g21  = g21 + f1 * c4 + f2 * c5;
        }
        double an12 = t.ann(i, j) * factor;
        t.t99_11(i,j) = a11 * an12;
        t.t99_12(i,j) = a12 * an12;
        t.t99_21(i,j) = a21 * an12;
        t.t99_22(i,j) = a22 * an12;
        t.t99_g11(i,j) = g11 * an12;
        t.t99_g12(i,j) = g12 * an12;
        t.t99_g21(i,j) = g21 * an12;
        t.t99_g22(i,j) = g22 * an12;
      }
    }

    int nm = nmax-m;
    for (int i = m; i < nmax; i ++)
    {
      int k1 = i-m;
      int kk1 = k1+nm;
      for (int j = m; j < nmax; j ++)
      {
        int k2 = j-m;
        int kk2 = k2+nm;
        dcplex ta11 = -t.t99_11(i,j);
        dcplex ta12 = -I*t.t99_12(i,j);
        dcplex ta21 = I*t.t99_21(i,j);
        dcplex ta22 = -t.t99_22(i,j);
        dcplex tg11 = -t.t99_g11(i,j);
        dcplex tg12 = -I*t.t99_g12(i,j);
        dcplex tg21 = I*t.t99_g21(i,j);
        dcplex tg22 = -t.t99_g22(i,j);
        t.ctt_q(k1, k2) = pi * ta21 + ppi * ta12;
        t.ctt_q(k1, kk2) = pi * ta11 + ppi * ta22;
        t.ctt_q(kk1, k2) = pi * ta22 + ppi * ta11;
        t.ctt_q(kk1,kk2) = pi * ta12 + ppi * ta21;
        t.ctt_rgq(k1, k2) = pi * tg21 + ppi * tg12;
        t.ctt_rgq(k1, kk2) = pi * tg11 + ppi * tg22;
        t.ctt_rgq(kk1, k2) = pi * tg22 + ppi * tg11;
        t.ctt_rgq(kk1,kk2) = pi * tg12 + ppi * tg21;
      }
    }

    TT(nm, t);
    return;
  }

  void Tmatrix::TT(int m, t_calculus_space &t)
  {
    int m2 = 2*m;

    Mdcplex tq(m2, m2);
    tq = subrange(t.ctt_q, 0,m2, 0,m2);

    if ((lapack_invert_inplace(tq)) != 0)
    {
      throw ("NASTY ERROR");
    }

    for (int i = 0; i < m2; i ++)
    {
      for (int j = 0; j < m2; j ++)
      {
        t.ct_t(i, j) = 0.0;
        for (int k = 0; k < m2; k ++)
        {
          t.ct_t(i, j) = t.ct_t(i, j) - (tq(k, j) * t.ctt_rgq(i, k));
        }
      }
    }

    return;
  }

  int Tmatrix::GetAmplitudeAndPhaseMatrices(double particle_alpha,
                                            double particle_beta,
                                            double incident_theta,
                                            double incident_phi,
                                            double scattered_theta,
                                            double scattered_phi,
                                            AmplitudeMatrix &A,
                                            PhaseMatrix &P)
  {
    if (particle_alpha  < 0.0 || particle_alpha  > 360.0 ||
        particle_beta   < 0.0 || particle_beta   > 180.0 ||
        incident_theta  < 0.0 || incident_theta  > 180.0 ||
        incident_phi    < 0.0 || incident_phi    > 360.0 ||
        scattered_theta < 0.0 || scattered_theta > 180.0 ||
        scattered_phi   < 0.0 || scattered_phi   > 180.0)
      return -1;

    Ampl(particle_alpha, particle_beta, incident_theta, incident_phi,
         scattered_theta, scattered_phi, A);

    const ddcplex I(0.0,1.0);
    P.z11=0.5*real((A.s11*conj(A.s11)+A.s12*conj(A.s12)
             +A.s21*conj(A.s21)+A.s22*conj(A.s22)));
    P.z12=0.5*real((A.s11*conj(A.s11)-A.s12*conj(A.s12)
             +A.s21*conj(A.s21)-A.s22*conj(A.s22)));
    P.z13=real(-A.s11*conj(A.s12)-A.s22*conj(A.s21));
    P.z14=real((I)*(A.s11*conj(A.s12)-A.s22*conj(A.s21)));
    P.z21=0.5*real((A.s11*conj(A.s11)+A.s12*conj(A.s12)
             -A.s21*conj(A.s21)-A.s22*conj(A.s22)));
    P.z22=0.5*real((A.s11*conj(A.s11)-A.s12*conj(A.s12)
             -A.s21*conj(A.s21)+A.s22*conj(A.s22)));
    P.z23=real(-A.s11*conj(A.s12)+A.s22*conj(A.s21));
    P.z24=real((I)*(A.s11*conj(A.s12)+A.s22*conj(A.s21)));
    P.z31=real(-A.s11*conj(A.s21)-A.s22*conj(A.s12));
    P.z32=real(-A.s11*conj(A.s21)+A.s22*conj(A.s12));
    P.z33=real(A.s11*conj(A.s22)+A.s12*conj(A.s21));
    P.z34=real((-I)*(A.s11*conj(A.s22)+A.s21*conj(A.s12)));
    P.z41=real((I)*(A.s21*conj(A.s11)+A.s22*conj(A.s12)));
    P.z42=real((I)*(A.s21*conj(A.s11)-A.s22*conj(A.s12)));
    P.z43=real((-I)*(A.s22*conj(A.s11)-A.s12*conj(A.s21)));
    P.z44=real(A.s22*conj(A.s11)-A.s12*conj(A.s21));

    return 0;
  }

  void Tmatrix::Ampl(double particle_alpha, double particle_beta,
                     double incident_theta, double incident_phi,
                     double scattered_theta, double scattered_phi,
                     AmplitudeMatrix &A)
  {
    const double degtorad=0.01745329251994329576;
    double littlenumber = 1e-8;

    double alph   = particle_alpha*degtorad;
    double bet    = particle_beta*degtorad;
    double thetl  = incident_theta*degtorad;
    double phil   = incident_phi*degtorad;
    double thetl1 = scattered_theta*degtorad;
    double phil1  = scattered_phi*degtorad;

    if (thetl  < M_PI_2) thetl  = thetl + littlenumber;
    if (thetl  > M_PI_2) thetl  = thetl - littlenumber;
    if (thetl1 < M_PI_2) thetl1 = thetl1 + littlenumber;
    if (thetl1 > M_PI_2) thetl1 = thetl1 - littlenumber;
    if (phil   < M_PI)   phil   = phil + littlenumber;
    if (phil   > M_PI)   phil   = phil - littlenumber;
    if (phil1  < M_PI)   phil1  = phil1 + littlenumber;
    if (phil1  > M_PI)   phil1  = phil1 - littlenumber;
    if (bet <= M_PI_2 && M_PI_2 - bet <= littlenumber) bet = bet - littlenumber;
    if (bet > M_PI_2  && bet - M_PI_2 <= littlenumber) bet = bet + littlenumber;

    double cb    = cos(bet);
    double sb    = sin(bet);
    double ct    = cos(thetl);
    double st    = sin(thetl);
    double cp    = cos(phil - alph);
    double sp    = sin(phil - alph);
    double ctp   = ct*cb + st*sb*cp;
    double thetp = acos(ctp);
    double cpp   = cb*st*cp - sb*ct;
    double spp   = st*sp;
    double phip  = atan(spp / cpp);
    if (phip > 0.0 && sp < 0.0) phip = phip + M_PI;
    if (phip < 0.0 && sp > 0.0) phip = phip + M_PI;
    if (phip < 0.0) phip = phip + 2.0 * M_PI;

    double ct1    = cos(thetl1);
    double st1    = sin(thetl1);
    double cp1    = cos(phil1 - alph);
    double sp1    = sin(phil1 - alph);
    double ctp1   = ct1*cb + st1*sb*cp1;
    double thetp1 = acos(ctp1);
    double cpp1   = cb*st1*cp1 - sb*ct1;
    double spp1   = st1*sp1;
    double phip1  = atan(spp1 / cpp1);
    if (phip1 > 0.0 && sp1 < 0.0) phip1 = phip1 + M_PI;
    if (phip1 < 0.0 && sp1 > 0.0) phip1 = phip1 + M_PI;
    if (phip1 < 0.0) phip1 = phip1 + 2.0 * M_PI;

    // Compute matrix beta, eq. (21)

    double ca = cos(alph);
    double sa = sin(alph);

    Mdouble b(3,3);
    b(0,0) = ca*cb;
    b(0,1) = sa*cb;
    b(0,2) = -sb;
    b(1,0) = -sa;
    b(1,1) = ca;
    b(1,2) = 0.0;
    b(2,0) = ca*sb;
    b(2,1) = sa*sb;
    b(2,2) = cb;

    // Compute matrices al and al1, eq. (14)

    cp  = cos(phil);
    sp  = sin(phil);
    cp1 = cos(phil1);
    sp1 = sin(phil1);

    Mdouble al(3,2);
    Mdouble al1(3,2);

    al(0,0)  = ct*cp;
    al(0,1)  = -sp;
    al(1,0)  = ct*sp;
    al(1,1)  = cp;
    al(2,0)  = -st;
    al(2,1)  = 0.0;
    al1(0,0) = ct1*cp1;
    al1(0,1) = -sp1;
    al1(1,0) = ct1*sp1;
    al1(1,1) = cp1;
    al1(2,0) = -st1;
    al1(2,1) = 0.0;

    // Compute matrices ap^(-1) and ap1^(-1), eq. (15)
    
    ct  = ctp;
    st  = sin(thetp);
    cp  = cos(phip);
    sp  = sin(phip);
    ct1 = ctp1;
    st1 = sin(thetp1);
    cp1 = cos(phip1);
    sp1 = sin(phip1);

    Mdouble ap(2,3);
    ap(0,0) = ct*cp;
    ap(0,1) = ct*sp;
    ap(0,2) = -st;
    ap(1,0) = -sp;
    ap(1,1) = cp;
    ap(1,2) = 0.0;

    Mdouble ap1(2,3);
    ap1(0,0) = ct1*cp1;
    ap1(0,1) = ct1*sp1;
    ap1(0,2) = -st1;
    ap1(1,0) = -sp1;
    ap1(1,1) = cp1;
    ap1(1,2) = 0.0;

    // Compute matrices r and r^(-1), eq. (13)

    Mddouble c(3,2);
    Mddouble r(2,2);
    Mddouble r1(2,2);

    c  = prod(b, al);
    r  = prod(ap, c);
    c  = prod(b, al1);
    r1 = prod(ap1, c);
    
    long double d = 1.0 / (r1(0,0)*r1(1,1)-r1(0,1)*r1(1,0));

    long double x = r1(0,0);
    r1(0,0) = r1(1,1)*d;
    r1(0,1) = -r1(0,1)*d;
    r1(1,0) = -r1(1,0)*d;
    r1(1,1) = x*d;

    const dcplex I(0.0, 1.0);
    Mdcplex cal(nmax, nmax);

    for (int nn=0; nn < nmax; nn ++)
    {
      for (int n=0; n < nmax; n ++)
      {
        dcplex cn(pow(I, (nn-n-1)));
        double dnn = (double) ((2*(n+1)+1)*(2*(nn+1)+1));
        dnn = dnn / (double) ((n+1)*(nn+1)*(n+2)*(nn+2));
        dnn = sqrt(dnn);
        cal(n,nn) = cn * dnn;
      }
    }

    A.s11 = 0.0;
    A.s12 = 0.0;
    A.s21 = 0.0;
    A.s22 = 0.0;

    double ph    = phip1-phip;
    for (int m = 0; m < nmax; m ++)
    {
      Vdouble dv1(nmax);
      Vdouble dv2(nmax);
      Vdouble dv01(nmax);
      Vdouble dv02(nmax);
      WignerDAm(ctp1, m, dv1, dv2);
      WignerDAm(ctp, m, dv01, dv02);
      double fc = 2.0 * cos(m*ph);
      double fs = 2.0 * sin(m*ph);
      for (int nn = 0; nn < nmax; nn ++)
      {
        double dv1nn = m*dv01(nn);
        double dv2nn = dv02(nn);
        for (int n = 0; n < nmax; n ++)
        {
          double dv1n = m*dv1(n);
          double dv2n = dv2(n);
          ddcplex ct11 = tmat_t11[m][n][nn];
          ddcplex ct22 = tmat_t22[m][n][nn];
          if (m == 0)
          {
            ddcplex cn = cal(n,nn) * dv2n * dv2nn;
            A.s11 = A.s11 + cn*ct22;
            A.s22 = A.s22 + cn*ct11;
          }
          else
          {
            ddcplex ct12 = tmat_t12[m][n][nn];
            ddcplex ct21 = tmat_t21[m][n][nn];
            ddcplex cn1 = cal(n,nn)*fc;
            ddcplex cn2 = cal(n,nn)*fs;
            long double d11 = dv1n*dv1nn;
            long double d12 = dv1n*dv2nn;
            long double d21 = dv2n*dv1nn;
            long double d22 = dv2n*dv2nn;
            A.s11 += (ct11*d11+ct21*d21+ct12*d12+ct22*d22)*cn1;
            A.s12 += (ct11*d12+ct21*d22+ct12*d11+ct22*d21)*cn2;
            A.s21 -= (ct11*d21+ct21*d11+ct12*d22+ct22*d12)*cn2;
            A.s22 += (ct11*d22+ct21*d12+ct12*d21+ct22*d11)*cn1;
          }
        }
      }
    }
    A.s11 /= xp;
    A.s12 /= xp;
    A.s21 /= xp;
    A.s22 /= xp;
    ddcplex cvv = A.s11*r(0,0)+A.s12*r(1,0);
    ddcplex cvh = A.s11*r(0,1)+A.s12*r(1,1);
    ddcplex chv = A.s21*r(0,0)+A.s22*r(1,0);
    ddcplex chh = A.s21*r(0,1)+A.s22*r(1,1);
    A.s11 = r1(0,0)*cvv+r1(0,1)*chv;
    A.s12 = r1(0,0)*cvh+r1(0,1)*chh;
    A.s21 = r1(1,0)*cvv+r1(1,1)*chv;
    A.s22 = r1(1,0)*cvh+r1(1,1)*chh;
    return;
  }

  std::ostream& operator<< (std::ostream& os, PhaseMatrix &P)
  {
    using boost::format;

    os << "Phase Matrix" << std::endl
       << "#################################################" << std::endl;
    format fm("%+14.6f %+14.6f %+14.6f %+14.6f\n");
    fm % P.z11 % P.z12 % P.z13 % P.z14;
    os << fm;
    fm % P.z21 % P.z22 % P.z23 % P.z24;
    os << fm;
    fm % P.z31 % P.z32 % P.z33 % P.z34;
    os << fm;
    fm % P.z41 % P.z42 % P.z43 % P.z44;
    os << fm;
    os << "#################################################" << std::endl;
    return os;
  }

  std::ostream& operator<< (std::ostream& os, AmplitudeMatrix &A)
  {
    using boost::format;

    os << "Amplitude Matrix" << std::endl
       << "#################################################" << std::endl;
    format fm("%+15.8f %+15.8f\n%+15.8f %+15.8f\n");
    fm % A.s11 % A.s12 % A.s21 % A.s22;
    os << fm;
    os << "#################################################" << std::endl;
    return os;
  }
}

#endif

#ifdef TESTME

#include <iostream>
#include <complex>
#include <boost/numeric/ublas/io.hpp>

using namespace cetemps;
using namespace boost::numeric::ublas;

int main(int argc, char *argv[])
{
  double wl = 6.283185;

  radiation rad(wl);
  std::cout << "Using radiation of wavelenght " << wl << std::endl;

  std::complex <double> rindex(1.5, 0.02);
  spheroid p(10.0000, 0.50, radius_equal_surface_area_sphere, rindex);
  std::cout << (scattering_particle &) p << std::endl;

  Tmatrix T;
  
  T.Setup((scattering_particle &) p, rad, 0.0001, 2);

  AmplitudeMatrix A;
  PhaseMatrix P;
  T.GetAmplitudeAndPhaseMatrices(145.00, 52.00, 56.00,
                                 114.00, 65.00, 128.00, A, P);

  std::cout << A << std::endl;
  std::cout << P << std::endl;

  return 0;
}

#endif
