/***************************************************************************
 *   Copyright (C) 2008-2009 Graziano Giuliani                             *
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

#ifndef __CETEMPS__SCATTERING__TMATRIX__
#define __CETEMPS__SCATTERING__TMATRIX__

#include <particle.h>
#include <radiation.h>
#include <simplecomplex.h>
#include "boost/multi_array.hpp"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <iostream>

namespace cetemps {

  using namespace boost::numeric::ublas;
  using namespace boost;

  class AmplitudeMatrix {
    public:
      AmplitudeMatrix( ) { s11 = s12 = s21 = s22 = 0.0; }
      ddcplex s11; ddcplex s12; ddcplex s21; ddcplex s22;
      friend std::ostream& operator<< (std::ostream& os, AmplitudeMatrix &A);
  };

  class PhaseMatrix {
    public:
      PhaseMatrix( ) { z11 = z12 = z13 = z14 = 0.0;
                       z21 = z22 = z23 = z24 = 0.0;
                       z31 = z32 = z33 = z34 = 0.0;
                       z41 = z42 = z43 = z44 = 0.0; }
      double z11; double z12; double z13; double z14;
      double z21; double z22; double z23; double z24;
      double z31; double z32; double z33; double z34;
      double z41; double z42; double z43; double z44;
      friend std::ostream& operator<< (std::ostream& os, PhaseMatrix &P);
  };

  typedef vector <double> Vdouble;
  typedef vector <dcplex> Vdcplex;
  typedef matrix <double> Mdouble;
  typedef matrix <long double> Mddouble;
  typedef matrix <dcplex> Mdcplex;

  typedef struct calculus_space {
    Vdouble x;
    Vdouble w;
    Vdouble an;
    Vdouble s;
    Vdouble ss;
    Vdouble r;
    Vdouble dr;
    Vdouble ddr;
    Vdcplex cdr;
    Mdouble ann;
    Mdcplex ct_t;
    Mdcplex t99_11;
    Mdcplex t99_12;
    Mdcplex t99_21;
    Mdcplex t99_22;
    Mdcplex t99_g11;
    Mdcplex t99_g12;
    Mdcplex t99_g21;
    Mdcplex t99_g22;
    Mdcplex ctt_q;
    Mdcplex ctt_rgq;
    Mdouble bessely;
    Mdouble besseldy;
    Mdouble besselj;
    Mdouble besseldj;
    Mdcplex cbesselj;
    Mdcplex cbesseldj;
  } t_calculus_space;

  class Tmatrix {
    public:
      void Setup(scattering_particle &p, radiation &rad,
                 double ddelt, int base_gaussian_divisions);
      double GetQscat( ) { return qsca; }
      double GetQext( ) { return qext; }

      int GetAmplitudeAndPhaseMatrices(double particle_alpha,
                                       double particle_beta,
                                       double incident_theta,
                                       double incident_phi,
                                       double scattered_theta,
                                       double scattered_phi,
                                       AmplitudeMatrix &A,
                                       PhaseMatrix &P);

    private:

      // The scattering T-matrix values
 
      multi_array <scplex, 3> tmat_t11;
      multi_array <scplex, 3> tmat_t12;
      multi_array <scplex, 3> tmat_t21;
      multi_array <scplex, 3> tmat_t22;

      int resize_problem(t_calculus_space &t);
      void WignerD(double x, int m, Vdouble &dv1, Vdouble &dv2);
      void WignerDAm(double x, int m, Vdouble &dv1, Vdouble &dv2);
      void Constant(scattering_particle &p, t_calculus_space &t);
      void Vary(scattering_particle &p, radiation &rad, t_calculus_space &t);
      void Tmatr0(t_calculus_space &t);
      void Tmatr(int m, t_calculus_space &t);
      void TT(int m, t_calculus_space &t);
      void Ampl(double particle_alpha, double particle_beta,
                double incident_theta, double incident_phi,
                double scattered_theta, double scattered_phi,
                AmplitudeMatrix &A);
      int nmax;
      int nmax2;
      int ngmax;
      int ngmax2;
      int ncheck;
      double ppi;
      dcplex pi;
      double qext;
      double qsca;
      double xp;
  };

}

#endif
