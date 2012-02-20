/*----------------------------------------------------------------------**
**
** run_len_dec.c -- Decode run-length encoded data
**
**----------------------------------------------------------------------**
**
** DESCRIPTION
**
** Decodes run-length encoded data
**
**
** USAGE:
**
** int run_len_dec(
** 	unsigned char *frombuf, int csize, 	-* Source buffer & size *-
** 	unsigned char *tobuf, int usize)	-* Dest. buffer & size  *-
**
**
** PROCESSING:
**
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
** $Source: /nfs/trmm/src/CVS/etor/run_len_dec.c,v $
**
** MODIFICATION RECORD
**
** $Log: run_len_dec.c,v $
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

/*----------------------------------------------------------------------*/
/* Application Headers                                                  */
/*----------------------------------------------------------------------*/
#include "run_enc.h"

/*----------------------------------------------------------------------*/
/* Local (Static) Variables                                             */
/*----------------------------------------------------------------------*/
static unsigned char tbl[32];
static int first=1;

/*----------------------------------------------------------------------*/
/* Main Function                                                        */
/*----------------------------------------------------------------------*/

int run_len_dec(
	unsigned char *frombuf, int csize, 	/* Source buffer & size */
	unsigned char *tobuf, int usize)	/* Dest. buffer & size  */
{
	register int i, pos=0;
	register unsigned char c;
	unsigned char c1;

	/*
	** Set up the table the first time through.
	** This table is used for a fast table lookup so that incoming
	** data can be used as an index into this table for the actual
	** data value
	*/
	if(first)
	{
		tbl[1]=0;
		for(i=2;i<26;i++)
			tbl[i]=(unsigned char)((10.63*(i-1))-5.0);

		first=0;
	}

	i=0;
	while(1)
	{
		register int len, to;

		/*
		** Get the character to work with
		*/
		c=frombuf[i++];

		/*
		** If it is header information then something is wrong
		*/
		if(c==TX_SOT||c==TX_EOT||c==TX_SOR||c==TX_SOH||c==TX_EOH)
		{
			fprintf(stderr, "found bad data in ray!\n");
			return(-1);
		}
		/*
		** A NULL signifies the end of a ray
		*/
		else if(c==TX_NULL)
		{
			/*
			** Clear out the rest of a short ray (the rest of the
			** ray must have been 0)
			*/
			for(;pos<usize; pos++)
				tobuf[pos]=0;

			break;
		}
		else
		{
			/*
			** A TX_2BBYTE signifies that this compression is
			** 9 or more bytes, so we need to get the next data
			** byte for a 6 bit compression value instead of 3
			*/
			if((c&0x1f)==TX_2BYTE)
			{
				register unsigned char value;

				c1=frombuf[i++];

				/*
				** Make sure second value is not a header
				** mark
				*/
				if(c1==TX_SOT||c1==TX_EOT||c1==TX_SOR||
					c1==TX_SOH||c1==TX_EOH||c1==TX_NULL)
				{
					fprintf(stderr,
						"found bad data in ray!\n");
					return(-1);
				}

				/*
				** length is 9-72 bytes
				*/
				len=(((c>>5)&0x7)+((c1&0xe0)>>2))+9;
				to=pos+len;
				if(pos+len>usize) to=usize;
				value=tbl[c1&0x1f];

				/*
				** Fill 'len' bytes with the value in the
				** table corresponding to the data value
				** in c1.
				*/
				for(;pos<to;pos++)
					tobuf[pos]=(char)value;
			}
			else
			{
				register unsigned char value;

				/*
				** Length is 1-8 bytes
				*/
				len=(((c>>5)&0x7)+1);
				to=pos+len;
				if(pos+len>usize) to=usize;
				value=tbl[c&0x1f];

				/*
				** Fill 'len' bytes with the value in the
				** table corresponding to the data value
				** in c.
				*/
				for(;pos<to;pos++)
					tobuf[pos]=(char)value;
			}
		}
	}
	return(pos);
}
/******************************* MODULE_END ********************************/
/*-END OF MODULE--------------------------------------------------------*/
