/*************************************************************************
**
** bitx -- bit-wise input and output routines
**
**************************************************************************
**
** HEADER INFORMATION:
**
**	Software Suite		- EDGE
**	Package			- Util
**
**	Reference number	- SP1/PGM/
**	Revision number		- $Revision: 1.1.1.1 $
**	Release state		- $State: Exp $
**	Author, designer	- Mark Nelson, Greenleaf Software
**	Source file		- bitx.c
**	Modification date	- $Date: 1999/03/15 16:26:02 $
**	Modified by		- $Author: merritt $
** 
** 
** COPYRIGHT NOTICE
**
**	Copyright (c) 1997 by Enterprise Electronics Corporation
**	All Rights Reserved
** 
** This program is copyright by Enterprise Electronics  Corporation,
** Enterprise,  Alabama,  USA  36330 (334) 347-3478.  It is licensed
** for use on a specific CPU and is not to be transferred or  other-
** wise divulged. Copies or modifications of this program must carry
** this copyright notice.
** 
** USAGE:
** void dump_bit_file(bf)
** BIT_STRM *bf;
**
** BIT_STRM *OpenOutputBitStream( name, len )
** unsigned char *name;
** int len;
**
** BIT_STRM *OpenInputBitStream( name, len )
** unsigned char *name;
** int len;
**
** int nextc(bit_strm)
** BIT_STRM *bit_strm;
**
** int outc(c, bit_strm)
** BIT_STRM *bit_strm;
** unsigned char c;
**
** int CloseOutputBitStream( bit_strm )
** BIT_STRM *bit_strm;
**
** int CloseInputBitStream( bit_strm )
** BIT_STRM *bit_strm;
**
** void StreamOutputBits( bit_strm, code, count )
** BIT_STRM *bit_strm;
** unsigned long code;
** int count;
**
** int StreamInputBit( bit_strm )
** BIT_STRM *bit_strm;
**
** unsigned long StreamInputBits( bit_strm, bit_count )
** BIT_STRM *bit_strm;
** int bit_count;
** 
** PROCESSING:
**
** This utility file contains all of the routines needed to impement
** bit oriented routines under either ANSI or K&R C.  It needs to be
** linked with every program used in the entire book.
**
** MODIFICATION RECORD:
**
** 10 Sep 92	KenB	From "The Data Compression Book", Marek Nelson
** 10 Sep 92	KenB	Modified to input from array, not files
**
**	$Log: bitx.c,v $
**	Revision 1.1.1.1  1999/03/15 16:26:02  merritt
**	Original import.
**	
 * Revision 1.2  1997/06/03  12:50:08  bobstaff
 * Fixed a pointer problem where a pointer was being used after being freed
 *
 * Revision 1.1  1995/12/04  16:56:21  edge
 * Initial revision
 *
 * Revision 1.1  1995/12/04  16:56:21  edge
 * Initial revision
 *
 * Revision 1.1  1994/09/09  01:43:59  stafford
 * Initial revision
 *
 * Revision 1.1  92/12/08  08:39:32  kenb
 * Initial revision
 * 
** 
------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* Includes:								*/
/*----------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "bitx.h"

/*----------------------------------------------------------------------*/
/* Defines:								*/
/*----------------------------------------------------------------------*/

#define PACIFIER_COUNT 2047

/*----------------------------------------------------------------------*/
/* Global Variables:							*/
/*----------------------------------------------------------------------*/
int lzw_verbose;

/*----------------------------------------------------------------------*/
void dump_bit_file(BIT_STRM *bf)
{
    fprintf(stderr, "\nptr         :%08lx\n", bf);
    fprintf(stderr, "addr        :%08lx\n", bf->addr);
    fprintf(stderr, "mask        :%02x\n", bf->mask);
    fprintf(stderr, "rack        :%08x\n", bf->rack);
    fprintf(stderr, "len         :%d\n", bf->len);
    fprintf(stderr, "count       :%d\n", bf->count);
    fprintf(stderr, "pacifier_cou:%d\n", bf->pacifier_counter);
}


/*----------------------------------------------------------------------*/
BIT_STRM *OpenOutputBitStream( unsigned char *name, int len )
{
	BIT_STRM *bit_strm;

	bit_strm = (BIT_STRM *) calloc( 1, sizeof( BIT_STRM ) );
	if ( bit_strm == NULL )
		return( bit_strm );
	bit_strm->addr = name;
	bit_strm->rack = 0;
	bit_strm->mask = 0x80;
	bit_strm->len = len;
	bit_strm->count = 0;
	bit_strm->pacifier_counter = 0;
	if(lzw_verbose > 2)
	{
		fprintf(stderr, "\nOpenOutputBitStream:");
		dump_bit_file(bit_strm);
	}
	return( bit_strm );
}

/*----------------------------------------------------------------------*/
BIT_STRM *OpenInputBitStream( unsigned char *name, int len )
{
	BIT_STRM *bit_strm;

	bit_strm = (BIT_STRM *) calloc( 1, sizeof( BIT_STRM ) );
	if ( bit_strm == NULL )
		return( bit_strm );
	bit_strm->addr = name;
	bit_strm->rack = 0;
	bit_strm->count = 0;
	bit_strm->mask = 0x80;
	bit_strm->len = len;
	bit_strm->pacifier_counter = 0;
	if(lzw_verbose > 2)
	{
		fprintf(stderr, "\nOpenInputBitStream:");
		dump_bit_file(bit_strm);
	}
	return( bit_strm );
}

/*----------------------------------------------------------------------*/
int nextc(BIT_STRM *bit_strm)
{
	if(bit_strm->len == bit_strm->count)
		return EOF;

	bit_strm->count++;

	return (int)*bit_strm->addr++;
}

/*----------------------------------------------------------------------*/
int outc(unsigned char c, BIT_STRM *bit_strm)
{
	if(bit_strm->len == bit_strm->count)
	{
		printf("buffer_full %d\n",bit_strm->len);
		return EOF;
	}

	bit_strm->count++;

	return (int)(*bit_strm->addr++ = c);
}


/*----------------------------------------------------------------------*/
int CloseOutputBitStream( BIT_STRM *bit_strm )
{
	int ret;

	if ( bit_strm->mask != 0x80 )
		if ( outc( bit_strm->rack, bit_strm ) != bit_strm->rack )
		{
			fprintf(stderr, "Fatal error in CloseBitStream!\n" );
			exit(1);
		}
	if(lzw_verbose > 2)
	{
		fprintf(stderr, "\nCloseOutputBitStream:");
		dump_bit_file(bit_strm);
	}
	ret=bit_strm->count;
	free( (char *) bit_strm );
	return ret;
}

/*----------------------------------------------------------------------*/
int CloseInputBitStream( BIT_STRM *bit_strm )
{
	int ret;

	if(lzw_verbose > 2)
	{
		fprintf(stderr, "\nCloseInputBitStream:");
		dump_bit_file(bit_strm);
	}
	ret=bit_strm->count;
	free( (char *) bit_strm );
	return ret;
}

/*----------------------------------------------------------------------*/
void StreamOutputBit( BIT_STRM *bit_strm, int bit )
{
	if ( bit )
		bit_strm->rack |= bit_strm->mask;
	bit_strm->mask >>= 1;
	if ( bit_strm->mask == 0 )
	{
		if ( outc( bit_strm->rack, bit_strm ) != bit_strm->rack )
		{
			fprintf(stderr, "Fatal error in StreamOutputBit!\n" );
			exit(1);
		}
		else 
			if(((bit_strm->pacifier_counter++ & PACIFIER_COUNT)
				== 0) && lzw_verbose)
				putc( '.', stderr );
		bit_strm->rack = 0;
		bit_strm->mask = 0x80;
	}
}

/*----------------------------------------------------------------------*/
void StreamOutputBits( BIT_STRM *bit_strm, unsigned long code, int count )
{
	unsigned long mask;

	mask = 1L << ( count - 1 );
	while ( mask != 0) {
		if ( mask & code )
			bit_strm->rack |= bit_strm->mask;
		bit_strm->mask >>= 1;
		if ( bit_strm->mask == 0 )
		{
			if( outc( bit_strm->rack, bit_strm ) != bit_strm->rack )
			{
				fprintf(stderr, 
				    "Fatal error in StreamOutputBits!\n" );
				exit(1);
			}
			else if(((bit_strm->pacifier_counter++&PACIFIER_COUNT)
					==0) && lzw_verbose)
				putc( '.', stderr );
			bit_strm->rack = 0;
			bit_strm->mask = 0x80;
		}
		mask >>= 1;
	}
}

/*----------------------------------------------------------------------*/
int StreamInputBit( BIT_STRM *bit_strm )
{
	int value;

	if ( bit_strm->mask == 0x80 ) {
		bit_strm->rack = nextc( bit_strm );
		if ( bit_strm->rack == EOF )
		{
			fprintf(stderr, "Fatal error in StreamInputBit!\n" );
			exit(1);
		}
		if (( ( bit_strm->pacifier_counter++ & PACIFIER_COUNT ) == 0 )
			&& lzw_verbose)
			putc( '.', stderr );
	}
	value = bit_strm->rack & bit_strm->mask;
	bit_strm->mask >>= 1;
	if ( bit_strm->mask == 0 )
		bit_strm->mask = 0x80;
	return( value ? 1 : 0 );
}

/*----------------------------------------------------------------------*/
unsigned long StreamInputBits( BIT_STRM *bit_strm, int bit_count )
{
	unsigned long mask;
	unsigned long return_value;

	mask = 1L << ( bit_count - 1 );
	return_value = 0;
	while ( mask != 0) {
		if ( bit_strm->mask == 0x80 )
		{
			bit_strm->rack = nextc( bit_strm );
			if ( bit_strm->rack == EOF )
			{
				fprintf(stderr, 
					"Fatal error in StreamInputBit!\n" );
				exit(1);
			}
			if(( ( bit_strm->pacifier_counter++ & PACIFIER_COUNT ) 
				== 0 ) && lzw_verbose)
				putc( '.', stderr );
		}
		if ( bit_strm->rack & bit_strm->mask )
			return_value |= mask;
		mask >>= 1;
		bit_strm->mask >>= 1;
		if ( bit_strm->mask == 0 )
			bit_strm->mask = 0x80;
	}
	return( return_value );
}
/*************************** End of BITIO.C **************************/
/*-END OF MODULE--------------------------------------------------------*/
