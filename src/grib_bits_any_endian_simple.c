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
#include <stdlib.h>
#include <stdint.h>
#ifdef __GNUC__
//define restrict will work with GCC and ICC for c89 support of restrict
#define RESTRICT __restrict
#else
#warn "no support for restrict keyword"
#define RESTRICT
#endif

static void encode_double_array_opt(size_t n_vals,const double* RESTRICT val,long bits_per_value,
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

    //now copy bytes per bytes (not optimal)
    /*for(size_t i=0;i<chunk;i++)
    {
      unsigned char* s=&tmp_uint[i];
      for(int i=0;i<bytes_to_copy;i++)
      {
        *p++=s[8-bytes_to_copy+i];
      }
    }*/
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
        #if 0
        for(i=0;i<chunk;i++)
        {
          *p++=(tmp_uint[i]>>40)&0xff;
          *((uint16_t*)p)=(tmp_uint[i]>>48)&0xffff;
          p+=2;
        }
        #else
        for(i=0;i<(chunk-1);i++)
        {
          uint32_t dw=tmp_uint[i]>>40&0xffffff;
          *((uint32_t*)p)=dw;
          p+=3;
        }
        *p++=(tmp_uint[chunk-1]>>40)&0xff;
        *((uint16_t*)p)=(tmp_uint[chunk-1]>>48)&0xffff;
        p+=2;
        #endif
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

/* A mask with x least-significant bits set, possibly 0 or >=32 */
/* -1UL is 1111111... in every bit in binary representation */
#define BIT_MASK1(x) \
    (((x) >= max_nbits) ? (unsigned long)-1UL : (1UL << (x)) - 1)

/**
 * decode an array of n_vals values from a octet-stream
 *
 * @param p input bitstream, for technical reasons put into octets
 * @param bitp current position in the bitstream
 * @param bitsPerValue number of bits needed to build a number (e.g. 8=byte, 16=short, 32=int, but also other sizes allowed)
 * @param n_vals number of values to decode
 * @param val output, values encoded as 32/64bit numbers
 */
int grib_decode_long_array(const unsigned char* p, long* bitp, long bitsPerValue,
                           size_t n_vals, long* val)
{
    unsigned long mask = BIT_MASK1(bitsPerValue);

    /* pi: position of bitp in p[]. >>3 == /8 */
    long pi = *bitp / 8;
    size_t i;
    /* number of useful bits in current byte */
    int usefulBitsInByte = 8 - (*bitp & 7);
    for (i = 0; i < n_vals; i++) {
        /* read at least enough bits (byte by byte) from input */
        long bitsToRead = bitsPerValue;
        long ret        = 0;
        while (bitsToRead > 0) {
            ret <<= 8;
            /*   ret += p[pi];         */
            /*   Assert( (ret & p[pi]) == 0 ); */
            ret = ret | p[pi];
            pi++;
            bitsToRead -= usefulBitsInByte;
            usefulBitsInByte = 8;
        }
        *bitp += bitsPerValue;
        /*fprintf(stderr, "%d %d %d %d %d\n", bitsPerValue, *bitp, pi, ret, bitsToRead);*/
        /* bitsToRead might now be negative (too many bits read) */
        /* remove those which are too much */
        ret >>= -1 * bitsToRead;
        /* remove leading bits (from previous value) */
        ret &= mask;
        val[i] = ret;

        usefulBitsInByte = -1 * bitsToRead; /* prepare for next round */
        if (usefulBitsInByte > 0) {
            pi--; /* reread the current byte */
        }
        else {
            usefulBitsInByte = 8; /* start with next full byte */
        }
    }
    return 0;
}

/**
 * decode an array of n_vals values from an octet-bitstream to double-representation
 *
 * @param p input bitstream, for technical reasons put into octets
 * @param bitp current position in the bitstream
 * @param bitsPerValue number of bits needed to build a number (e.g. 8=byte, 16=short, 32=int, but also other sizes allowed)
 * @param n_vals number of values to decode
 * @param val output, values encoded as 32/64bit numbers
 */
int grib_decode_double_array(const unsigned char* p, long* bitp, long bitsPerValue,
                             double reference_value, double s, double d,
                             size_t n_vals, double* val)
{
    long i               = 0;
    unsigned long lvalue = 0;
    double x;

#if 0
    /* slow reference code */
    int j=0;
    for(i=0;i < n_vals;i++) {
        lvalue=0;
        for(j=0; j< bitsPerValue;j++){
            lvalue <<= 1;
            if(grib_get_bit( p, *bitp)) lvalue += 1;
            *bitp += 1;
        }
        x=((lvalue*s)+reference_value)*d;
        val[i] = (double)x;
    }
#endif
    if (bitsPerValue % 8 == 0) {
        /* See ECC-386 */
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
            x      = ((lvalue * s) + reference_value) * d;
            val[i] = (double)x;
            /*  *bitp += bitsPerValue * n_vals; */
        }
    }
    else {
        unsigned long mask = BIT_MASK1(bitsPerValue);

        /* pi: position of bitp in p[]. >>3 == /8 */
        long pi = *bitp / 8;
        /* some bits might of the current byte at pi might be used */
        /* by the previous number usefulBitsInByte gives remaining unused bits */
        /* number of useful bits in current byte */
        int usefulBitsInByte = 8 - (*bitp & 7);
        for (i = 0; i < n_vals; i++) {
            /* value read as long */
            long bitsToRead = 0;
            lvalue          = 0;
            bitsToRead      = bitsPerValue;
            /* read one byte after the other to lvalue until >= bitsPerValue are read */
            while (bitsToRead > 0) {
                lvalue <<= 8;
                lvalue += p[pi];
                pi++;
                bitsToRead -= usefulBitsInByte;
                usefulBitsInByte = 8;
            }
            *bitp += bitsPerValue;
            /* bitsToRead is now <= 0, remove the last bits */
            lvalue >>= -1 * bitsToRead;
            /* set leading bits to 0 - removing bits used for previous number */
            lvalue &= mask;

            usefulBitsInByte = -1 * bitsToRead; /* prepare for next round */
            if (usefulBitsInByte > 0) {
                pi--; /* reread the current byte */
            }
            else {
                usefulBitsInByte = 8; /* start with next full byte */
            }
            /* scaling and move value to output */
            x      = ((lvalue * s) + reference_value) * d;
            val[i] = (double)x;
        }
    }
    return 0;
}

int grib_decode_double_array_complex(const unsigned char* p, long* bitp, long nbits, double reference_value, double s, double* d, size_t size, double* val)
{
    return GRIB_NOT_IMPLEMENTED;
}

int grib_encode_long_array(size_t n_vals, const long* val, long bits_per_value, unsigned char* p, long* off)
{
    size_t i                   = 0;
    unsigned long unsigned_val = 0;
    unsigned char* encoded     = p;
    if (bits_per_value % 8) {
        for (i = 0; i < n_vals; i++) {
            unsigned_val = val[i];
            grib_encode_unsigned_longb(encoded, unsigned_val, off, bits_per_value);
        }
    }
    else {
        for (i = 0; i < n_vals; i++) {
            int blen     = 0;
            blen         = bits_per_value;
            unsigned_val = val[i];
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

int grib_encode_double_array(size_t n_vals, const double* val, long bits_per_value, double reference_value, double d, double divisor, unsigned char* p, long* off)
{
    size_t i                   = 0;
    unsigned long unsigned_val = 0;
    unsigned char* encoded     = p;
    double x;
    if (bits_per_value % 8) {
        for (i = 0; i < n_vals; i++) {
            x            = (((val[i] * d) - reference_value) * divisor) + 0.5;
            unsigned_val = (unsigned long)x;
            grib_encode_unsigned_longb(encoded, unsigned_val, off, bits_per_value);
        }
    }
    else if(bits_per_value==8 || bits_per_value==16 
	    || bits_per_value==24 || bits_per_value==32) {
      encode_double_array_opt(n_vals,val,bits_per_value,reference_value,d,divisor,p,off);
    } else {
      for (i = 0; i < n_vals; i++) {
	int blen     = 0;
	blen         = bits_per_value;
	x            = ((((val[i] * d) - reference_value) * divisor) + 0.5);
	unsigned_val = (unsigned long)x;
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
