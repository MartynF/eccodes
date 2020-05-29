/*
 * (C) Copyright 2005- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
 * virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
 */

/***************************************************************************
 *   Enrico Fucile  - 19.06.2007                                           *
 *                                                                         *
 ***************************************************************************/

#include <stdint.h>
#ifdef __GNUC__
//define restrict will work with GCC and ICC for c89 support of restrict
#define RESTRICT __restrict
#else
#warn "no support for restrict keyword"
#define RESTRICT
#endif

int grib_decode_long_array(const unsigned char* p, long* bitp, long bitsPerValue,
                           size_t n_vals, long* val)
{
    long i               = 0;
    unsigned long lvalue = 0;

    if (bitsPerValue % 8) {
        int j = 0;
        for (i = 0; i < n_vals; i++) {
            lvalue = 0;
            for (j = 0; j < bitsPerValue; j++) {
                lvalue <<= 1;
                if (grib_get_bit(p, *bitp))
                    lvalue += 1;
                *bitp += 1;
            }
            val[i] = lvalue;
        }
    }
    else {
        int bc;
        int l    = bitsPerValue / 8;
        size_t o = *bitp / 8;

        for (i = 0; i < n_vals; i++) {
            lvalue = 0;
            lvalue <<= 8;
            lvalue |= p[o++];

            for (bc = 1; bc < l; bc++) {
                lvalue <<= 8;
                lvalue |= p[o++];
            }
            val[i] = lvalue;
        }
        *bitp += bitsPerValue * n_vals;
    }

    return 0;
}


int grib_decode_double_array(const unsigned char* p, long* bitp, long bitsPerValue,
                             double reference_value, double s, double d,
                             size_t n_vals, double* val)
{
    size_t i             = 0;
    unsigned long lvalue = 0;

    if (bitsPerValue % 8) {
        int j = 0;
        for (i = 0; i < n_vals; i++) {
            lvalue = 0;
            for (j = 0; j < bitsPerValue; j++) {
                lvalue <<= 1;
                if (grib_get_bit(p, *bitp))
                    lvalue += 1;
                *bitp += 1;
            }
            val[i] = (double)(((lvalue * s) + reference_value) * d);
        }
    }
    else {
        int bc;
        int l    = bitsPerValue / 8;
        size_t o = 0;

        for (i = 0; i < n_vals; i++) {
            lvalue = 0;
            lvalue <<= 8;
            lvalue |= p[o++];

            for (bc = 1; bc < l; bc++) {
                lvalue <<= 8;
                lvalue |= p[o++];
            }
            val[i] = (double)(((lvalue * s) + reference_value) * d);
        }
    }

    return 0;
}

int grib_decode_double_array_complex(const unsigned char* p, long* bitp, long nbits, double reference_value, double s, double* d, size_t size, double* val)
{
    return GRIB_NOT_IMPLEMENTED;
}

static void grib_encode_double_array_opt(size_t n_vals,const double* RESTRICT val,long bits_per_value,
  const double reference_value,const double d,const double divisor,unsigned char* p,long *off)
{
  const size_t block=64;
  size_t ib=0;
  const int bytes_to_copy=bits_per_value/8;

  while(ib<n_vals)
  {
    size_t chunk=(n_vals-ib)>block?block:(n_vals-ib);

    double tmp_flt[block];
    size_t i;

    //#pragma vector always
    for( i=0;i<chunk;i++)
    {
      tmp_flt[i]= ((((val[i+ib]*d)-reference_value)*divisor)+0.5);
    }

    uint64_t tmp_uint[block];
    //compute integers in a vectorizable way

    for(i=0;i<chunk;i++)
    {
      tmp_uint[i]=(uint64_t) (tmp_flt[i]);
    }

    //for each integer, byte swap it

    for(i=0;i<chunk;i++)
    {
      tmp_uint[i]=__builtin_bswap64(tmp_uint[i]);
    }
    switch (bytes_to_copy)
    {
      case 1:
        for(i=0;i<chunk;i++)
        {
          *p++=(tmp_uint[i]>>56)&0xff;
        }
      break;
      case 2:
        for(i=0;i<chunk;i++)
        {
          *((uint16_t*)p)=(tmp_uint[i]>>48)&0xffff;
          p+=2;
        }
      break;
      case 3:
        for(i=0;i<chunk;i++)
        {
          *p++=(tmp_uint[i]>>40)&0xff;
          *((uint16_t*)p)=(tmp_uint[i]>>48)&0xffff;
          p+=2;
        }
      break;
      case 4:
        for(i=0;i<chunk;i++)
        {
          *((uint32_t*)p)=(tmp_uint[i]>>32);
          p+=4;
        }
      break;
    default:
      for(i=0;i<chunk;i++)
      {
        unsigned char* s=&tmp_uint[i];
        int j;
        for(j=0;j<bytes_to_copy;j++)
        {
          *p++=s[8-bytes_to_copy+j];
        }
      }
    }

    ib+=chunk;
  }

  *off+=bits_per_value*n_vals;
}

int grib_encode_double_array(size_t n_vals, const double* val, long bits_per_value, double reference_value, double d, double divisor, unsigned char* p, long* off)
{
    size_t i                   = 0;
    unsigned long unsigned_val = 0;
    unsigned char* encoded     = p;
    if (bits_per_value % 8) {
      for (i = 0; i < n_vals; i++) {
	unsigned_val = (unsigned long)((((val[i] * d) - reference_value) * divisor) + 0.5);
	grib_encode_unsigned_longb(encoded, unsigned_val, off, bits_per_value);
      }
    }
    else  if(bits_per_value==8 || bits_per_value==16 
	     || bits_per_value==24 || bits_per_value==32) {
      grib_encode_double_array_opt(n_vals,val,bits_per_value,reference_value,d,divisor,p,off);
    } else {	  
      for (i = 0; i < n_vals; i++) {
	int blen     = 0;
	blen         = bits_per_value;
	unsigned_val = (unsigned long)((((val[i] * d) - reference_value) * divisor) + 0.5);
	while (blen >= 8) {
	  blen -= 8;
	  *encoded = (unsigned_val >> blen);
	  encoded++;
	  *off += 8;
	}
      }
    }
    return GRIB_SUCCESS;
}

int grib_encode_double_array_complex(size_t n_vals, double* val, long nbits, double reference_value,
                                     double* scal, double d, double divisor, unsigned char* p, long* bitp)
{
    return GRIB_NOT_IMPLEMENTED;
}
