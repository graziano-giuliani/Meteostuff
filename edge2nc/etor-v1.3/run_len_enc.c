/*----------------------------------------------------------------------**
**
** run_len_enc.c - Compress data run-length
**
**----------------------------------------------------------------------**
**
** DESCRIPTION
**
** Compress data run-length.
**
**
**
** USAGE:
**
** int run_len_enc(unsigned char *frombuf, int usize, -* Input buffer & size *-
** 		unsigned char *tobuf, int csize)   -* Output buffer & size *-
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
** $Source: /nfs/trmm/src/CVS/etor/run_len_enc.c,v $
**
** MODIFICATION RECORD
**
** $Log: run_len_enc.c,v $
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
#include "stdio.h"


/*----------------------------------------------------------------------*/
/* Application Headers                                                  */
/*----------------------------------------------------------------------*/
#include "run_enc.h"


/*----------------------------------------------------------------------*/
/* Local (Static) Variables                                             */
/*----------------------------------------------------------------------*/
static int first=1;
static unsigned char table[256];

/*----------------------------------------------------------------------*/
/* Main Function                                                        */
/*----------------------------------------------------------------------*/

int run_len_enc(unsigned char *frombuf, int usize, /* Input buffer & size */
		unsigned char *tobuf, int csize)   /* Output buffer & size */
{
	register int i, j, pos=0;
	unsigned char tmptot, tmpval;

	if(first)
	{
		table[0]=1;
		for(i=1;i<256;i++)
			table[i]=(unsigned char)((i/10.63)+2.0);

		first=0;
	}

	i=1;
	tmptot=1;
	tmpval=frombuf[0];

	while(i<usize)
	{
		if(tmpval==0 && i>1)
		{
			int done=0, doit=1;

			j=i;
			while(j<usize&&!done)
			{
				if(frombuf[j++]!=0x00)
				{
					done=1;
					doit=0;
				}
			}

			if(doit)
				break;
		}

		tmpval=table[tmpval];

		while((tmpval==table[ frombuf[i] ])&&tmptot<72&&i<usize)
		{
			tmptot++;
			i++;
		}

		if(tmptot<=8)
		{
			tmptot--;
			tobuf[pos++]= (unsigned char)(((tmptot<<5)&0xe0)
				|(tmpval));
		}
		else
		{
			tmptot-=9;
			tobuf[pos++]=(unsigned char)(0x1f|((tmptot&0x7)<<5));
			tobuf[pos++]=(unsigned char)(((tmptot&0x38)<<2)|
				(tmpval));
		}

		if(i<usize)
			tmpval=frombuf[i++];
		tmptot=1;
	}
	tobuf[pos++]=TX_NULL;

	return(pos);
}
/******************************* MODULE_END ********************************/
/*-END OF MODULE--------------------------------------------------------*/
