/*----------------------------------------------------------------------**
**
** crc.c - Cyclic Redundancy Check
**
**----------------------------------------------------------------------**
**
** DESCRIPTION
**
** Cyclic Redundancy Check .  Classical.
**
**
**
** USAGE:
**
** unsigned int crc( 		-* CRC is returned			*-
** 	void *buffer, 		-* Buffer to check			*-
** 	unsigned int count)	-* Size of buffer			*-
**
**
** PROCESSING:
**
** This is the routine used to calculate the 32 bit CRC of a block of data.
** This is done by processing the input buffer using the coefficient table
** that was created when the program was initialized.  This routine takes
** an input value as a seed, so that a running calculation of the CRC can
** be used as blocks are read and written.
**
** COPYRIGHT NOTICE
**
**      Copyright (c) 1993 by Lassen Research
**      All Rights Reserved
**
**      This program is copyright by Lassen Research, Chico, California,
**      95928, (916) 343-6421. It is licensed for use on a specific CPU
**      and is not to be transferred or otherwise divulged.  Copies or
**      modifications of this program must carry this copyright notice.
**
** HEADER INFOMATION
**
**      Software Suite          - EDGE
**      Package                 -
**      Reference number                - SP1/PGM/
**      Revision number         - $Revision: 1.1.1.1 $
**      Release State           - $State: Exp $
**      Author, designer                -
** Modification Date            - $Date: 1999/03/15 16:26:02 $
** Modified by                  - $Author: merritt $
** $Source: /nfs/trmm/src/CVS/etor/crc.c,v $
**
** MODIFICATION RECORD
**
** $Log: crc.c,v $
** Revision 1.1.1.1  1999/03/15 16:26:02  merritt
** Original import.
**
 * Revision 1.1  1995/12/04  16:56:21  edge
 * Initial revision
 *
 * Revision 1.1  1995/12/04  16:56:21  edge
 * Initial revision
 *
 * Revision 1.1  1994/09/09  01:43:59  stafford
 * Initial revision
 *
**
**----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* Feature Test Switches                                                */
/*----------------------------------------------------------------------*/
#define _POSIX_SOURCE   1

/*----------------------------------------------------------------------*/
/* System Headers                                                       */
/*----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

/*----------------------------------------------------------------------*/
/* Macros                                                               */
/*----------------------------------------------------------------------*/
#define CRC_MASK           0xFFFFFFFFL
#define CRC32_POLYNOMIAL   0xEDB88320L

/*----------------------------------------------------------------------*/
/* Local (Static) Variables                                             */
/*----------------------------------------------------------------------*/
static unsigned int crc_table[ 256 ];    /* This array holds the CRC	*/
				      /* table used to calculate the 32 */
				      /* bit CRC values.                */

/*----------------------------------------------------------------------*/
/* Local Function                                                       */
/*----------------------------------------------------------------------*/
/*
 * This routine simply builds the coefficient table used to calculate
 * 32 bit CRC values throughout this program.  The 256 int word table
 * has to be set up once when the program starts.  Alternatively, the
 * values could be hard coded in, which would offer a miniscule improvement
 * in overall performance of the program.
 */

static void build_crc_table()
{
	int i;
    	int j;
	unsigned int value;

	for ( i = 0; i <= 255 ; i++ ) 
	{
		value = i;
		for ( j = 8 ; j > 0; j-- ) 
		{
			if ( value & 1 )

				value = ( value >> 1 ) ^ CRC32_POLYNOMIAL;
	    		else
				value >>= 1;
		}
		crc_table[ i ] = value;
    	}
}

/*
 * This is the routine used to calculate the 32 bit CRC of a block of data.
 * This is done by processing the input buffer using the coefficient table
 * that was created when the program was initialized.  This routine takes
 * an input value as a seed, so that a running calculation of the CRC can
 * be used as blocks are read and written.
 */
/*----------------------------------------------------------------------*/
/* Main Function                                                        */
/*----------------------------------------------------------------------*/

unsigned int crc( 		/* CRC is returned			*/
	void *buffer, 		/* Buffer to check			*/
	unsigned int count)	/* Size of buffer			*/
{
    	unsigned char *p = (unsigned char *) buffer;
    	unsigned int temp1;
    	unsigned int temp2;
    	unsigned int crc=0; 

	if (crc_table[1]==0)
		build_crc_table();

    	while ( count-- != 0 ) 
	{
		temp1 = ( crc >> 8 ) & 0x00FFFFFFL;
		temp2 = crc_table[ ( (int) crc ^ *p++ ) & 0xff ];
		crc = temp1 ^ temp2;
    	}
    	return( crc );
}

/*-END OF MODULE--------------------------------------------------------*/
