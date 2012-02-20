/***************************************************************************
 *   Copyright (C) 2008-2010 by Graziano Giuliani                          *
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

#ifndef __GRIBLIB_H__
#define __GRIBLIB_H__

#include <stddef.h>

#define UNDEFINED 9.999e20

#define UINT2(a,b)     ((int) ((a << 8) + (b)))
#define UINT4(a,b,c,d) ((int) ((a << 24) + (b << 16) + (c << 8) + (d)))
#define UINT3(a,b,c)   ((int) ((a << 16) + (b << 8) + (c)))
#define INT2(a,b)      ((1-(int) ((unsigned) (a & 0x80) >> 6)) * \
                         (int) (((a & 0x7f) << 8) + b))
#define INT3(a,b,c)    ((1-(int) ((unsigned) (a & 0x80) >> 6)) * \
                         (int) (((a & 127) << 16)+(b<<8)+c))

void gbits(unsigned char *, long *, long, size_t, long, long);
void gbit(unsigned char *, long *, long, size_t);
long uint2(unsigned char *p);
long uint4(unsigned char *p);
unsigned int uint_n(unsigned char *p, int n);
long int4(unsigned char *p);
long int2(unsigned char *p);
int int_n(unsigned char *p, int n);

double int_power(double x, int y);
double ibm2flt(unsigned char *ibm);
double scaled2flt(int scale_factor, int scale_value);
int missing_points(unsigned char *bitmap, int n);
float ieee2flt(unsigned char *ieee);
float ieee2flt_nan(unsigned char *ieee);
int sub_missing_values(unsigned char *p, float *missing1, float *missing2);

typedef struct {
  long ndata;
  long packing;
  long precision;
  long bitmap_flag;
  unsigned char *mask;
  double reference;
  double bin_scale;
  double dec_scale;
  int nbits;
  long ngroups;
  int ctable_5_4;
  int ctable_5_5;
  int ctable_5_6;
  int extra_octets;
  long ref_group_width;
  long ref_group_length;
  long nbit_group_width;
  int group_length_factor;
  int nbits_group_len;
  long len_last;
  int npnts;
  int n_sub_missing;
  float missing1;
  float missing2;
  unsigned int mv, mvl;
  double *rlv;
} t_g2_pinfo;

typedef struct {
  long ndata;
  int nbits;
  double reference;
  double scale;
  unsigned char *bitmap;
} t_g1_pinfo;

int unpk_grib1(float *data, unsigned char *bds, unsigned char *bitmap,
                int n_bits, long n, double ref, double scale);

int unpk_grib2(float *data, t_g2_pinfo &pcki, unsigned char *sec7);

#endif
