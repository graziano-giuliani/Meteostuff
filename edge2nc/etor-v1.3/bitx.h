/************************** Start of BITIO.H *************************/

#ifndef BITIO_H
#define BITIO_H

#include <stdio.h>

typedef struct bit_file {
    unsigned char *addr;
    unsigned char mask;
    int rack;
    int len;
    int count;
    int pacifier_counter;
} BIT_STRM;

#ifdef __STDC__

BIT_STRM     *OpenInputBitStream( unsigned char *addr, int len );
BIT_STRM     *OpenOutputBitStream( unsigned char *addr, int len );
void          StreamOutputBit( BIT_STRM *bit_file, int bit );
void          StreamOutputBits( BIT_STRM *bit_file,
                          unsigned long code, int count );
int           StreamInputBit( BIT_STRM *bit_file );
unsigned long StreamInputBits( BIT_STRM *bit_file, int bit_count );
int          CloseInputBitStream( BIT_STRM *bit_file );
int          CloseOutputBitStream( BIT_STRM *bit_file );

#else   /* __STDC__ */

BIT_STRM     *OpenInputBitStream();
BIT_STRM     *OpenOutputBitStream();
void          StreamOutputBit();
void          StreamOutputBits();
int           StreamInputBit();
unsigned long StreamInputBits();
int          CloseInputBitStream();
int          CloseOutputBitStream();

#endif  /* __STDC__ */

#endif  /* BITIO_H */

/*************************** End of BITIO.H **************************/

