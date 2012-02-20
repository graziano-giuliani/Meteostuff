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

#ifndef __CETEMPS__HYDROPSD__
#define __CETEMPS__HYDROPSD__

namespace cetemps {

  typedef enum {
    exponential_psd = 0,
    generalized_gamma_psd = 1,
    wrf_gamma_psd = 2
  } t_enum_PSD;

  class hydroPSD {
    public:
      virtual double sum(double d1, double d2) = 0;
      virtual double ndd(double d) = 0;
      virtual double moment(int order) = 0;
      void SetExtremes(double d1, double d2);
      double minvalue( );
      double maxvalue( );
      double *regspaced(int nvals)
      {
        double *interval = new double[nvals];

        double lmin = lower_bound;
        double lmax = upper_bound;
        double delta = (lmax - lmin) / (nvals-1);
        for (int i = 0; i < nvals; i ++)
          interval[i] = lmin+i*delta;
        return interval;
      }
      double *log10spaced(int nvals)
      {
        double *interval = new double[nvals];

        double lmin = log10(lower_bound);
        double lmax = log10(upper_bound);
        double delta = (lmax - lmin) / (nvals-1);
        for (int i = 0; i < nvals; i ++)
          interval[i] = pow(10, lmin+i*delta);
        return interval;
      }
      double *logspaced(int nvals)
      {
        double *interval = new double[nvals];

        double lmin = log(lower_bound);
        double lmax = log(upper_bound);
        double delta = (lmax - lmin) / (nvals-1);
        for (int i = 0; i < nvals; i ++)
          interval[i] = exp(lmin+i*delta);
        return interval;
      }
    private:
      double lower_bound;
      double upper_bound;
  };

  class ExponentialPSD : public hydroPSD {
    public:
      ExponentialPSD(double intercept, double slope);
      virtual double ndd(double d);
      virtual double sum(double d1, double d2);
      virtual double moment(int order);
    private:
      double intercept;
      double slope;
  };

  class GeneralizedGammaPSD : public hydroPSD {
    public:
      GeneralizedGammaPSD(double intercept, double shape, double mean_diameter);
      virtual double ndd(double d);
      virtual double sum(double d1, double d2);
      virtual double moment(int order);
    private:
      double intercept;
      double shape;
      double dzero;
      double slope;
      double fmu;
  };

  class WRFGammaPSD : public hydroPSD {
    public:
      WRFGammaPSD(double intercept, double shape, double slope);
      virtual double ndd(double d);
      virtual double sum(double d1, double d2);
      virtual double moment(int order);
    private:
      double intercept;
      double shape;
      double slope;
  };

  class WRFDoubleGammaPSD : public hydroPSD {
    public:
      WRFDoubleGammaPSD(double q, double t,
                        double kap0, double kap1,
                        double lam0, double lam1, double shape);
      virtual double ndd(double d);
      virtual double sum(double d1, double d2);
      virtual double moment(int order);
    private:
      double kap0 , kap1;
      double lam0 , lam1;
      double shape;
      double mult;
  };

}

#endif
