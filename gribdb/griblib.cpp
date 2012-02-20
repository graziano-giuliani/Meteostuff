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

#include <iostream>

#include <griblib.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>

#include <jasper/jasper.h>
#include <png.h>

void rd_bitstream(unsigned char *p, unsigned int *u, int n_bits, int n);
void rd_var_len_bitstream(unsigned char *p, unsigned int *u, int n);
int dec_png(unsigned char *image, long *w, long *h, unsigned char *data);
int unpk_complex(float *data, t_g2_pinfo &pkci, unsigned char *sec7);
int unpk_run_length(float *data, t_g2_pinfo &pkci, unsigned char *sec7);

static int bitsum[256] = {
    8, 7, 7, 6, 7, 6, 6, 5, 7, 6, 6, 5, 6, 5, 5, 4,
    7, 6, 6, 5, 6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3,
    7, 6, 6, 5, 6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3,
    6, 5, 5, 4, 5, 4, 4, 3, 5, 4, 4, 3, 4, 3, 3, 2,
    7, 6, 6, 5, 6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3,
    6, 5, 5, 4, 5, 4, 4, 3, 5, 4, 4, 3, 4, 3, 3, 2,
    6, 5, 5, 4, 5, 4, 4, 3, 5, 4, 4, 3, 4, 3, 3, 2,
    5, 4, 4, 3, 4, 3, 3, 2, 4, 3, 3, 2, 3, 2, 2, 1,
    7, 6, 6, 5, 6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3,
    6, 5, 5, 4, 5, 4, 4, 3, 5, 4, 4, 3, 4, 3, 3, 2,
    6, 5, 5, 4, 5, 4, 4, 3, 5, 4, 4, 3, 4, 3, 3, 2,
    5, 4, 4, 3, 4, 3, 3, 2, 4, 3, 3, 2, 3, 2, 2, 1,
    6, 5, 5, 4, 5, 4, 4, 3, 5, 4, 4, 3, 4, 3, 3, 2,
    5, 4, 4, 3, 4, 3, 3, 2, 4, 3, 3, 2, 3, 2, 2, 1,
    5, 4, 4, 3, 4, 3, 3, 2, 4, 3, 3, 2, 3, 2, 2, 1,
    4, 3, 3, 2, 3, 2, 2, 1, 3, 2, 2, 1, 2, 1, 1, 0};

static void user_read_data(png_structp , png_bytep , png_uint_32 );

double ibm2flt(unsigned char *ibm)
{
  int positive, power;
  unsigned int abspower;
  long int mant;
  double value, exp;

  mant = (ibm[1] << 16) + (ibm[2] << 8) + ibm[3];
  if (mant == 0) return 0.0;

  positive = (ibm[0] & 0x80) == 0;
  power = (int) (ibm[0] & 0x7f) - 64;
  abspower = power > 0 ? power : -power;

  /* calc exp */
  exp = 16.0;
  value = 1.0;
  while (abspower)
  {
    if (abspower & 1) value *= exp;
    exp = exp * exp;
    abspower >>= 1;
  }

  if (power < 0) value = 1.0 / value;
  value = value * mant / 16777216.0;
  if (positive == 0) value = -value;
  return value;
}

void gbit(unsigned char *in, long *iout, long iskip, size_t nbyte)
{
  gbits(in, iout, iskip, nbyte, 0, 1);
  return;
}

//
// Get bits - unpack bits:  Extract arbitrary size values from a packed
//    bit string, right justifying each value in the unpacked iout array.
//
// *in   = pointer to character array input
// *iout = pointer to unpacked array output
// iskip = initial number of bits to skip
// nbyte = number of bits to take
// nskip = additional number of bits to skip on each iteration
// n     = number of iterations
//
void gbits(unsigned char *in, long *iout, long iskip, size_t nbyte,
           long nskip, long n)
{
  long i,tbit,bitcnt,ibit,itmp;
  long nbit,index;
  static long ones[]={1,3,7,15,31,63,127,255};

  // nbit is the start position of the field in bits
  nbit = iskip;
  for (i=0; i<n; i++)
  {
    bitcnt = nbyte;
    index=nbit/8;
    ibit=nbit%8;
    nbit = nbit + nbyte + nskip;

    // first byte
    tbit= ( bitcnt < (8-ibit) ) ? bitcnt : 8-ibit;  // find min
    itmp = (int)*(in+index) & ones[7-ibit];
    if (tbit != 8-ibit) itmp >>= (8-ibit-tbit);
    index++;
    bitcnt = bitcnt - tbit;
    // now transfer whole bytes
    while (bitcnt >= 8)
    {
      itmp = itmp<<8 | (int)*(in+index);
      bitcnt = bitcnt - 8;
      index++;
    }
    // get data from last byte
    if (bitcnt > 0)
      itmp = ( itmp << bitcnt ) |
             ( ((int)*(in+index) >> (8-bitcnt)) & ones[bitcnt-1] );
    *(iout+i) = itmp;
  }
  return;
}

long uint2(unsigned char *p)
{
  return (p[0] << 8) + p[1];
}

long uint4(unsigned char *p)
{
  return ((p[0] << 24) + (p[1] << 16) + (p[2] << 8) + p[3]);
}

unsigned  int uint_n(unsigned char *p, int n)
{
  unsigned int i;
  i = 0;
  while (n-- > 0) {
    i = (i << 8) + *p++;
  }
  return i;
}

long int2(unsigned char *p)
{
  long i;
  if (p[0] & 0x80)
    i = -(((p[0] & 0x7f) << 8) + p[1]);
  else
    i = (p[0] << 8) + p[1];
  return i;
}

long int4(unsigned char *p)
{
  long i;
  if (p[0] & 0x80)
    i = -(((p[0] & 0x7f) << 24) + (p[1] << 16) + (p[2] << 8) + p[3]);
  else
    i = (p[0] << 24) + (p[1] << 16) + (p[2] << 8) + p[3];
  return i;
}

int int_n(unsigned char *p, int n)
{
  int i, sign;

  if (n == 0) return 0;
  sign = *p;
  i = *p++ & 127;
  while (n-- > 1) {
    i = i * 256 + (int) *p++;
  }
  if (sign & 0x80) i = -i;
  return i;
}

void rd_bitstream(unsigned char *p, unsigned int *u, int n_bits, int n)
{
  static unsigned int ones[]={0, 1,3,7,15,31,63,127,255};

  unsigned int tbits;
  int i, t_bits;

  if (n_bits > 25) throw "WHOAAA! Error in rd_bitstream";

  if (n_bits == 0) {
    for (i = 0; i < n; i++) {
      u[i] = 0;
    }
    return;
  }

  tbits = t_bits = 0;

  for (i = 0; i < n; i++) {
    tbits = tbits & ones[t_bits];
    while (t_bits < n_bits) {
      t_bits += 8;
      tbits = (tbits << 8) | *p++;
    }
    t_bits -= n_bits;
    u[i] = (int) (tbits >> t_bits);
  }
  return;
}

void rd_var_len_bitstream(unsigned char *p, unsigned int *u, int n)
{
  static unsigned int ones[]={0, 1,3,7,15,31,63,127,255};
  unsigned int tbits, n_bits;
  int i, t_bits;

  tbits = t_bits = 0;
  for (i = 0; i < n; i++)
  {
    n_bits = u[i];
    if (n_bits == 0) {
      u[i] = 0;
    }
    else if (n_bits > 25) {
      throw "WHOAAA! Error in rd_var_len_bitstream";
    }
    else {
      tbits = tbits & ones[t_bits];
      while (t_bits < (int) n_bits) {
        tbits = (tbits << 8) | *p++;
        t_bits += 8;
      }
      t_bits -= n_bits;
      u[i] = (int) (tbits >> t_bits);
    }
  }
  return;
}

double int_power(double x, int y)
{
  double value;

  if (y < 0)
  {
    y = -y;
    x = 1.0 / x;
  }
  value = 1.0;

  while (y)
  {
    if (y & 1)
      value *= x;
    x = x * x;
    y >>= 1;
  }
  return value;
}

double scaled2flt(int scale_factor, int scale_value)
{
   if (scale_factor == 0) return (double) scale_value;
   if (scale_factor < 0) return scale_value * int_power(10.0, -scale_factor);
   return scale_value / int_power(10.0, scale_factor);
}

int missing_points(unsigned char *bitmap, int n)
{
  int count;
  unsigned int tmp;
  if (bitmap == NULL) return 0;

  count = 0;
  while (n >= 8)
  {
    tmp = *bitmap++;
    n -= 8;
    count += bitsum[tmp];
  }
  tmp = *bitmap | ((1 << (8 - n)) - 1);
  count += bitsum[tmp];
  return count;
}


int unpk_grib1(float *data, unsigned char *bds, unsigned char *bitmap,
                int n_bits, long n, double ref, double scale)
{
  unsigned char *bits;

  static const unsigned int mask[] = {0,1,3,7,15,31,63,127,255};
  static const unsigned int map_masks[8] = {128, 64, 32, 16, 8, 4, 2, 1};
  static const double shift[9] = {  1.0,   2.0,   4.0,   8.0,  16.0,
                                   32.0,  64.0, 128.0, 256.0};
  int i, mask_idx, t_bits, c_bits, j_bits;
  unsigned int j, map_mask, tbits, jmask, bbits;
  double jj;

  if ((bds[3] & 64) != 0) // Do not handle complex packing
    return -1;

  if (bds[3] & 128) // Harmonic
  {
    bits = bds + 15;
    /* fill in global mean */
    *data++ = (ibm2flt(bds+11));
    n -= 1;
  }
  else
    bits = bds + 11;

  tbits = bbits = 0;

  /* assume integer has 32+ bits */
  if (n_bits <= 25)
  {
    jmask = (1 << n_bits) - 1;
    t_bits = 0;

    if (bitmap)
    {
      for (i = 0; i < n; i++)
      {
        /* check bitmap */
        mask_idx = i & 7;
        if (mask_idx == 0) bbits = *bitmap++;
        if ((bbits & map_masks[mask_idx]) == 0)
        {
          *data++ = UNDEFINED;
          continue;
        }
        while (t_bits < n_bits)
        {
          tbits = (tbits * 256) + *bits++;
          t_bits += 8;
        }
        t_bits -= n_bits;
        j = (tbits >> t_bits) & jmask;
        *data++ = ref + scale*j;
      }
    }
    else
    {
      for (i = 0; i < n; i++)
      {
        if (n_bits - t_bits > 8)
        {
          tbits = (tbits << 16) | (bits[0] << 8) | (bits[1]);
          bits += 2;
          t_bits += 16;
        }
        while (t_bits < n_bits)
        {
          tbits = (tbits * 256) + *bits++;
          t_bits += 8;
        }
        t_bits -= n_bits;
        data[i] = (tbits >> t_bits) & jmask;
      }
      /* at least this vectorizes :) */
      for (i = 0; i < n; i++)
        data[i] = ref + scale*data[i];
    }
  }
  else
  {
    /* older unoptimized code, not often used */
    c_bits = 8;
    map_mask = 128;
    while (n-- > 0)
    {
      if (bitmap)
      {
        j = (*bitmap & map_mask);
        if ((map_mask >>= 1) == 0)
        {
          map_mask = 128;
          bitmap++;
        }
        if (j == 0)
        {
          *data++ = UNDEFINED;
          continue;
        }
      }

      jj = 0.0;
      j_bits = n_bits;
      while (c_bits <= j_bits)
      {
        if (c_bits == 8)
        {
          jj = jj * 256.0  + (double) (*bits++);
          j_bits -= 8;
        }
        else
        {
          jj = (jj * shift[c_bits]) + (double) (*bits & mask[c_bits]);
          bits++;
          j_bits -= c_bits;
          c_bits = 8;
        }
      }
      if (j_bits)
      {
        c_bits -= j_bits;
        jj = (jj * shift[j_bits]) + (double) ((*bits >> c_bits) & mask[j_bits]);
      }
      *data++ = ref + scale*jj;
    }
  }
  return 0;
}

float ieee2flt(unsigned char *ieee)
{
  double fmant;
  int exp;

  if ((ieee[0] & 127) == 0 && ieee[1] == 0 && ieee[2] == 0 && ieee[3] == 0)
    return (float) 0.0;

  exp = ((ieee[0] & 127) << 1) + (ieee[1] >> 7);
  fmant = (double) ((int) ieee[3] + (int) (ieee[2] << 8) +
          (int) ((ieee[1] | 128) << 16));
  if (ieee[0] & 128) fmant = -fmant;
  return (float) (ldexp(fmant, (int) (exp - 128 - 22)));
}

float ieee2flt_nan(unsigned char *ieee)
{
  double fmant;
  int exp;

  if ((ieee[0] & 127) == 0 && ieee[1] == 0 && ieee[2] == 0 && ieee[3] == 0)
    return (float) 0.0;

  exp = ((ieee[0] & 127) << 1) + (ieee[1] >> 7);
  if (exp == 255) return (float) UNDEFINED;
  fmant = (double) ((int) ieee[3] + (int) (ieee[2] << 8) +
          (int) ((ieee[1] | 128) << 16));
  if (ieee[0] & 128) fmant = -fmant;
  return (float) (ldexp(fmant, (int) (exp - 128 - 22)));
}

int sub_missing_values(unsigned char *p, float *missing1, float *missing2)
{
  int ct55 = p[22];
  int ct51 = p[20];

  if (ct51 == 0)               // ieee
  {
    if (p[23] == 255 && p[24] == 255 && p[25] == 255 && p[26] == 255)
      *missing1 = UNDEFINED;
    else *missing1 = ieee2flt(p+23);
    if (ct55 == 2)
    {
      if (p[27] == 255 && p[28] == 255 && p[29] == 255 && p[30] == 255)
        *missing1 = UNDEFINED;
      else *missing2 = ieee2flt(p+27);
    }
  }
  else if (ct51 == 1)          // integer
  {
    if (p[23] == 255 && p[24] == 255 && p[25] == 255 && p[26] == 255)
      *missing1 = UNDEFINED;
    else *missing1 = (float) int4(p+23);
    if (ct55 == 2)
    {
      if (p[27] == 255 && p[28] == 255 && p[29] == 255 && p[30] == 255)
        *missing1 = UNDEFINED;
      else *missing2 = (float) int4(p+27);
    }
  }
  return ct55;
}

int unpk_complex(float *data, t_g2_pinfo &pkci, unsigned char *sec7)
{
  double ref_val = pkci.reference*pkci.dec_scale;
  double factor = pkci.bin_scale*pkci.dec_scale;

  if (pkci.packing == 3 && (pkci.ctable_5_6 != 1 && pkci.ctable_5_6 != 2))
    throw "Unsupported spatial differencing order";

  unsigned char *mask_pointer = pkci.mask;

  if (pkci.ngroups == 0)
  {
    if (pkci.bitmap_flag == 255)
    {
      for (int i = 0; i < pkci.ndata; i++) data[i] = ref_val;
      return 0;
    }
    if (pkci.bitmap_flag == 0 || pkci.bitmap_flag == 254)
    {
      unsigned char mask = 0;
      for (int i = 0; i < pkci.ndata; i++)
      {
        if ((i & 7) == 0) mask = *mask_pointer++;
        data[i] = (mask & 128) ?  ref_val : UNDEFINED;
        mask <<= 1;
      }
      return 0;
    }
    throw "Unknown bitmap used : Cannot unpack";
  }

  unsigned int *group_refs, *group_widths, *group_lengths, *udata;
  group_refs = group_widths = group_lengths = udata = 0;
  // allocate group widths and group lengths
  group_refs = new unsigned int[pkci.ngroups];
  group_widths = new unsigned int[pkci.ngroups];
  group_lengths = new unsigned int[pkci.ngroups];
  udata = new unsigned int[pkci.npnts];
  if (group_refs == 0    || group_widths == 0 ||
      group_lengths == 0 || udata == 0)
    throw "Complex Unpacking fatal error";

  // read any extra values
  unsigned char *d = sec7+5;
  int extra_vals[2];
  int min_val = 0;
  extra_vals[0] = 0;
  extra_vals[1] = 0;
  if (pkci.extra_octets)
  {
    extra_vals[0] = uint_n(d, pkci.extra_octets);
    d +=  pkci.extra_octets;
    if (pkci.ctable_5_6 == 2)
    {
      extra_vals[1] = uint_n(d, pkci.extra_octets);
      d +=  pkci.extra_octets;
    }
    min_val = int_n(d, pkci.extra_octets);
    d +=  pkci.extra_octets;
  }

  try
  {
    unsigned int i, j, k;

    // read the group reference values
    rd_bitstream(d, group_refs, pkci.nbits, pkci.ngroups);
    d += (pkci.nbits*pkci.ngroups + 7)/8;

    // read the group widths
    rd_bitstream(d,group_widths,pkci.nbit_group_width,pkci.ngroups);
    for (i = 0; i < pkci.ngroups; i++)
      group_widths[i] += pkci.ref_group_width;
    d += (pkci.ngroups * pkci.nbit_group_width + 7) / 8;

    // read the group lengths
    if (pkci.ctable_5_4 == 1)
    {
      rd_bitstream(d, group_lengths, pkci.nbits_group_len, pkci.ngroups-1);
      d += (pkci.ngroups * pkci.nbits_group_len + 7) / 8;
      for (i = 0; i < pkci.ngroups-1; i++)
        group_lengths[i] = group_lengths[i] * pkci.group_length_factor +
                           pkci.ref_group_length;
      group_lengths[pkci.ngroups-1] = pkci.len_last;
    }
    // do a check for number of grid points and size
    for (i = j = k = 0; i < pkci.ngroups; i++) {
        j += group_lengths[i];
        k += group_lengths[i]*group_widths[i];
    }
    if (j != (unsigned int) pkci.npnts) throw "bad complex packing";
    if (d + (k+7)/8 - sec7 != uint4(sec7))
        throw "complex unpacking size mismatch";

    // set nbits of each data point
    for (i = j = 0; i < pkci.ngroups; i++) {
      for (k = 0; k < group_lengths[i]; k++) {
        udata[j++] = group_widths[i];
      }
    }
    // get data points values (uint)
    rd_var_len_bitstream(d, udata, pkci.npnts);

    // handle substitute, missing values and reference value
    if (pkci.n_sub_missing == 0) {
      for (i = j = 0; i < pkci.ngroups; i++) {
        for (k = 0; k < group_lengths[i]; k++) {
          udata[j++] += group_refs[i];
        }
      }
    }
    else if (pkci.n_sub_missing == 1) {
      for (i = j = 0; i < pkci.ngroups; i++) {
        if (group_widths[i] == 0) {
          unsigned int m1 = (1 << pkci.nbits) - 1;
          if (m1 == group_refs[i]) {
            for (k = 0; k < group_lengths[i]; k++)
              udata[j++] = INT_MAX;
          }
          else {
             for (k = 0; k < group_lengths[i]; k++)
               udata[j++] += group_refs[i];
          }
        }
        else {
          unsigned int m1 = (1 << group_widths[i]) - 1;
          for (k = 0; k < group_lengths[i]; k++) {
            if (udata[j] == m1) udata[j] = INT_MAX;
            else udata[j] += group_refs[i];
            j++;
          }
        }
      }
    }
    else if (pkci.n_sub_missing == 2)
    {
      for (i = j = 0; i < pkci.ngroups; i++) {
        if (group_widths[i] == 0) {
          unsigned int m1 = (1 << pkci.nbits) - 1;
          unsigned int m2 = m1 - 1;
          if (m1 == group_refs[i] || m2 == group_refs[i]) {
            for (k = 0; k < group_lengths[i]; k++)
              udata[j++] = INT_MAX;
          }
          else {
            for (k = 0; k < group_lengths[i]; k++)
               udata[j++] += group_refs[i];
          }
        }
        else {
          unsigned int m1 = (1 << group_widths[i]) - 1;
          unsigned int m2 = m1 - 1;
          for (k = 0; k < group_lengths[i]; k++) {
            if (udata[j] == m1 || udata[j] == m2) udata[j] = INT_MAX;
            else udata[j] += group_refs[i];
            j++;
          }
        }
      }
    }
    if (pkci.packing == 3)
    {
      if (pkci.ctable_5_6 == 1) {
        unsigned int last = extra_vals[0];
        long i = 0;
        while (i < pkci.npnts) {
          if (udata[i] == INT_MAX) i++;
          else {
            udata[i++] = extra_vals[0];
            break;
          }
        }
        while (i < pkci.npnts) {
          if (udata[i] == INT_MAX) i++;
          else {
            udata[i] += last + min_val;
            last = udata[i++];
          }
        }
      }
      else if (pkci.ctable_5_6 == 2) {
        unsigned int penultimate = extra_vals[0];
        unsigned int last = extra_vals[1];
        long i = 0;
        while (i < pkci.npnts) {
          if (udata[i] == INT_MAX) i++;
          else {
            udata[i++] = extra_vals[0];
            break;
          }
        }
        while (i < pkci.npnts) {
          if (udata[i] == INT_MAX) i++;
          else {
            udata[i++] = extra_vals[1];
            break;
          }
        }
        while (i < pkci.npnts) {
          if (udata[i] == INT_MAX) i++;
          else {
            udata[i] =  udata[i] + min_val + last + last - penultimate;
            penultimate = last;
            last = udata[i++];
          }
        }
      }
      else
        throw "Unsupported code table 5.6";
    }
    // convert to float

    if (pkci.bitmap_flag == 255) {
      for (i = 0; i < pkci.ndata; i++) {
        data[i] = (udata[i] == INT_MAX) ? UNDEFINED :
                  ref_val + udata[i] * factor;
      }
    }
    else if (pkci.bitmap_flag == 0 || pkci.bitmap_flag == 254)
    {
      int k = 0;
      unsigned char mask = 0;
      for (i = 0; i < pkci.ndata; i++) {
        if ((i & 7) == 0) mask = *mask_pointer++;
        if (mask & 128) {
          if (udata[k] == INT_MAX) data[i] = UNDEFINED;
          else data[i] = ref_val + udata[k] * factor;
          k++;
        }
        else data[i] = UNDEFINED;
        mask <<= 1;
      }
    }
    else
      throw "Unknown bitmap!";
  }
  catch (char const* str)
  {
    std::cerr << str << std::endl;
    delete [] group_refs;
    delete [] group_widths;
    delete [] group_lengths;
    delete [] udata;
    throw "Error in Unpacking complex";
  }

  delete [] group_refs;
  delete [] group_widths;
  delete [] group_lengths;
  delete [] udata;

  return 0;
}

int unpk_run_length(float *data, t_g2_pinfo &pkci, unsigned char *sec7)
{
  const unsigned int mask[] = {128,64,32,16,8,4,2,1};

  unsigned int size_compressed_data = uint4(sec7)-5;
  unsigned int nvals = ( size_compressed_data * 8) / pkci.nbits;

  double *levels = pkci.rlv;
  unsigned int *vals = new unsigned int[nvals];

  try
  {
    rd_bitstream(sec7+5, vals, pkci.nbits, nvals);

    unsigned int i, j, ncheck;
    i = j = ncheck = 0;

    long range = (1 << pkci.nbits) - 1 - pkci.mv;
    if (range <= 0) throw "test rlp range error in unpk_run_length";
    unsigned char *mask_pointer = pkci.mask;
    if (pkci.bitmap_flag != 0   &&
        pkci.bitmap_flag != 254 &&
        pkci.bitmap_flag != 255)
      throw "Unknown bitmap used : Cannot unpack";
    while (i < nvals)
    {
      if (vals[i] > pkci.mv) throw "RLP val greater than Max value ?";
      unsigned int v = vals[i++];
      int n = 1;
      int factor = 1;
      while (vals[i] > pkci.mv && i < nvals) {
        n += factor * (vals[i]-pkci.mv-1);
        factor = factor * range;
        i++;
      }
      ncheck += n;
      if (ncheck > pkci.ndata)
        throw "More points than needed in unpk_run_length";
      if (pkci.bitmap_flag == 255) {
        for (int k = 0; k < n; k++) data[j++] = levels[v];
      }
      else 
      {
        for (int k = 0; k < n; k++) {
          while (mask_pointer[j >> 3] & mask[j & 7]) {
            data[j++] = UNDEFINED;
          }
          data[j++] = levels[v];
        }
      }
      if (j != pkci.ndata) throw "Bitmap problem in unpk_run_length";
    }
  }
  catch (char const* str)
  {
    std::cerr << str << std::endl;
    delete [] vals;
    throw "Fatal error in unpk_run_length";
  }

  delete [] vals; 
  return 0;
}

struct png_stream {
   unsigned char *stream_ptr;     /*  location to write PNG stream  */
   long stream_len;               /*  number of bytes written       */
};
typedef struct png_stream png_stream;

void user_read_data(png_structp png_ptr, png_bytep data, png_uint_32 length)
{
  char *ptr;
  long offset;
  png_stream *mem;

  mem=(png_stream *)png_get_io_ptr(png_ptr);
  ptr= (char *)((void *)mem->stream_ptr);
  offset=mem->stream_len;
  memcpy(data,ptr+offset,length);
  mem->stream_len += length;
  return;
}

int dec_png(unsigned char *image, long *w, long *h, unsigned char *data)
{
  int interlace,color,compres,filter,bit_depth;
  long j,k,n,bytes,clen;
  png_structp png_ptr;
  png_infop info_ptr,end_info;
  png_bytepp row_pointers;
  png_stream read_io_ptr;
  png_uint_32 h32, w32;

  if (png_sig_cmp(image,0,8) != 0)
    return (-3); // Invalid stream

  png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING,
                     (png_voidp)NULL, NULL, NULL);
  if (!png_ptr)
    return (-1);

  info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr)
  {
    png_destroy_read_struct(&png_ptr,(png_infopp)NULL,(png_infopp)NULL);
    return (-2);
  }
  end_info = png_create_info_struct(png_ptr);
  if (!end_info)
  {
    png_destroy_read_struct(&png_ptr,(png_infopp)info_ptr,(png_infopp)NULL);
    return (-2);
  }
  if (setjmp(png_jmpbuf(png_ptr)))
  {
    png_destroy_read_struct(&png_ptr, &info_ptr,&end_info);
    return (-3);
  }
  read_io_ptr.stream_ptr=(unsigned char *) ((png_voidp) image);
  read_io_ptr.stream_len=0;
  png_set_read_fn(png_ptr,(voidp)&read_io_ptr,(png_rw_ptr)user_read_data);
  //  Read and decode PNG stream
  png_read_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);
  // Get pointer to each row of image data
  row_pointers = png_get_rows(png_ptr, info_ptr);
  (void)png_get_IHDR(png_ptr, info_ptr, &w32, &h32,
               &bit_depth, &color, &interlace, &compres, &filter);
  *h = h32;
  *w = w32;
  if ( color == PNG_COLOR_TYPE_RGB ) {
    bit_depth=24;
  }
  else if ( color == PNG_COLOR_TYPE_RGB_ALPHA ) {
    bit_depth=32;
  }

  // Copy image data to output string
  n=0;
  bytes=bit_depth/8;
  clen=(w32)*bytes;
  for (j=0;j<(int)h32;j++) {
    for (k=0;k<clen;k++) {
      data[n]=*(row_pointers[j]+k);
      n++;
    }
  }
  png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
  return 0;
}

int unpk_grib2(float *data, t_g2_pinfo &pkci, unsigned char *sec7)
{
  if (pkci.bitmap_flag != 0   &&
      pkci.bitmap_flag != 254 &&
      pkci.bitmap_flag != 255)
    throw "Unknown bitmap used : Cannot unpack";

  unsigned char *mask_pointer = pkci.mask;

  if (pkci.packing == 4)
  {
    if (pkci.precision != 1)
      throw "Unsupported IEEE packing precision";
    if (pkci.bitmap_flag == 255)
    {
      for (int i = 0; i < pkci.ndata; i++)
        data[i] = ieee2flt_nan(sec7+5+i*4);
      return 0;
    }
    if (pkci.bitmap_flag == 0 || pkci.bitmap_flag == 254)
    {
      unsigned char *ieee = sec7+5;
      unsigned char mask = 0;
      for (int i = 0; i < pkci.ndata; i++)
      {
        if ((i & 7) == 0) mask = *mask_pointer++;
        if (mask & 128)
        {
          data[i] = ieee2flt_nan(ieee);
          ieee += 4;
        }
        else
          data[i] = UNDEFINED;
        mask <<= 1;
      }
      return 0;
    }
    throw "Unknown bitmap used : Cannot unpack";
  }
  else if (pkci.packing == 0)
  {
    float tmp;
    if (pkci.nbits == 0)
    {
      tmp = pkci.reference*pkci.dec_scale;
      if (pkci.bitmap_flag == 255)
      {
        for (int i = 0; i < pkci.ndata; i++)
          data[i] = tmp;
        return 0;
      }
      if (pkci.bitmap_flag == 0 || pkci.bitmap_flag == 254)
      {
        unsigned char mask = 0;
        for (int i = 0; i < pkci.ndata; i++)
        {
          if ((i & 7) == 0) mask = *mask_pointer++;
          data[i] = (mask & 128) ?  tmp : UNDEFINED;
          mask <<= 1;
        }
        return 0;
      }
    }
    unpk_grib1(data, sec7+5, mask_pointer, pkci.nbits, pkci.ndata,
           pkci.reference*pkci.dec_scale, pkci.bin_scale*pkci.dec_scale);
    return 0;
  }
  else if (pkci.packing == 2 || pkci.packing == 3)
    return unpk_complex(data, pkci, sec7);
  else if (pkci.packing == 200)
    return unpk_run_length(data, pkci, sec7);
  else if (pkci.packing == 40 || pkci.packing == 40000)
  {
    if (pkci.nbits == 0)
    {
      double tmp = pkci.reference*pkci.dec_scale;
      if (pkci.bitmap_flag == 255) {
        for (int i = 0; i < pkci.ndata; i++) {
          data[i] = tmp;
        }
        return 0;
      }
      else if (pkci.bitmap_flag == 0 || pkci.bitmap_flag == 254)
      {
        unsigned char mask = 0;
        for (int i = 0; i < pkci.ndata; i++) {
          if ((i & 7) == 0) mask = *mask_pointer++;
          data[i] = (mask & 128) ? tmp : UNDEFINED;
          mask <<= 1;
        }
        return 0;
      }
      throw "Unknown bitmap";
    }
    // decode jpeg2000
    jas_image_t *image = 0;
    char *opts = 0;
    jas_stream_t *jpcstream = 0;
    jas_image_cmpt_t *pcmpt = 0;
    jas_matrix_t *jas_data = 0;
    try
    {
      jpcstream=jas_stream_memopen((char*) sec7+5, uint4(sec7)-5);
      image = jpc_decode(jpcstream, opts);
      if (image == NULL) throw "Jpeg2000 decoding error";
      pcmpt = image->cmpts_[0];
      if (image->numcmpts_ != 1) throw "Color image found instead of grayscale";
      jas_data=jas_matrix_create(jas_image_height(image),
                                 jas_image_width(image));
      jas_image_readcmpt(image,0,0,0,
               jas_image_width(image), jas_image_height(image),jas_data);
      // transfer data

      int k = pkci.ndata - pcmpt->height_ * pcmpt->width_;
      for (int i=0; i < pcmpt->height_; i++) {
        for (int j=0; j < pcmpt->width_; j++) {
          data[k++] = (((jas_data->rows_[i][j])*
               pkci.bin_scale)+pkci.reference)*pkci.dec_scale;
        }
      }
      if (pkci.bitmap_flag == 0 || pkci.bitmap_flag == 254) {
        k = pkci.ndata - pcmpt->height_ * pcmpt->width_;
        unsigned char mask = 0;
        for (int i = 0; i < pkci.ndata; i++) {
          if ((i & 7) == 0) mask = *mask_pointer++;
          data[i] = (mask & 128) ? data[k++] : UNDEFINED;
          mask <<= 1;
        }
      }
    }
    catch (char const* str)
    {
      jas_matrix_destroy(jas_data);
      jas_stream_close(jpcstream);
      jas_image_destroy(image);
      std::cerr << str << std::endl;
      throw "JPEG 2000 unpack error";
    }
    jas_matrix_destroy(jas_data);
    jas_stream_close(jpcstream);
    jas_image_destroy(image);
    return 0;
  }
  else if (pkci.packing == 41)
  {
    if (pkci.nbits == 0) {
      double tmp = pkci.reference*pkci.dec_scale;
      if (pkci.bitmap_flag == 255) {
        for (int i = 0; i < pkci.ndata; i++) {
          data[i] = tmp;
        }
        return 0;
      }
      if (pkci.bitmap_flag == 0 || pkci.bitmap_flag == 254) {
        unsigned char mask = 0;
        for (int i = 0; i < pkci.ndata; i++) {
          if ((i & 7) == 0) mask = *mask_pointer++;
          data[i] = (mask & 128) ?  tmp : UNDEFINED;
          mask <<= 1;
        }
        return 0;
      }
      throw "Unsupported bitmap found";
    }

    unsigned char *c = new unsigned char[pkci.ndata*4];
    if (c == 0) throw "Memory error in PNG decoding";
    try
    {
      long width,height;
      if (dec_png(sec7+5, &width, &height, c) != 0)
        throw "Error in decoding PNG coded data";
      // Check if correct size
      if (pkci.mask == 0)
      {
        if (pkci.ndata != width*height)
          throw "PNG image size mismatch";
      }
      else
      {
        if (pkci.ndata != width*height +
                missing_points(mask_pointer, pkci.ndata))
          throw "PNG image size mismatch";
      }

      // Unpack coded GRIB1 data
      unpk_grib1(data, c, mask_pointer, pkci.nbits, pkci.ndata,
             pkci.reference*pkci.dec_scale, pkci.bin_scale*pkci.dec_scale);
    }
    catch (char const* str)
    {
      std::cerr << str << std::endl;
      delete [] c;
      throw "PNG unpack error";
    }
    delete [] c;
    return 0;
  }
  throw "Unsupported packing";
  return 0;
}
