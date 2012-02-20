/*----------------------------------------------------------------------**
**
** copy_buf.c -- Use memcpy to move data
**
**----------------------------------------------------------------------**
**
** DESCRIPTION
**
** Calls memcpy(3) to move data.  Used in compress jump 
** vector for no compression
**
**
**
** USAGE:
**
** int copy_buff(  unsigned char *inbuf,	-* Source buffer	*-
** 		int inbufsize, 			-* Source buffer size	*-
** 		unsigned char *outbuf, 		-* Dest. buffer		*-
** 		int outbufsize)			-* Dest. buffer	size	*-
**
** PROCESSING:
**
** memcpy(3) is called
**
** COPYRIGHT NOTICE
**
**	Copyright (c) 1993 by Lassen Research
**	All Rights Reserved
**
**	This program is copyright by Lassen Research, Chico, California,
**	95928, (916) 343-6421. It is licensed for use on a specific CPU
**	and is not to be transferred or otherwise divulged.  Copies or
**	modifications of this program must carry this copyright notice.
**
** HEADER INFOMATION
**
**	Software Suite 		- EDGE
**	Package			-
**	Reference number		- SP1/PGM/
**	Revision number		- $Revision: 1.1.1.1 $
**	Release State		- $State: Exp $
**	Author, designer		-
** Modification Date		- $Date: 1999/03/15 16:26:02 $
** Modified by			- $Author: merritt $
** $Source: /nfs/trmm/src/CVS/etor/copy_buff.c,v $
**
** MODIFICATION RECORD
**
** $Log: copy_buff.c,v $
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
#define _POSIX_SOURCE	1

/*----------------------------------------------------------------------*/
/* System Headers                                                       */
/*----------------------------------------------------------------------*/
#include	<string.h>	/* String functions			*/

/*----------------------------------------------------------------------*/
/* Application Headers                                                  */
/*----------------------------------------------------------------------*/
#include "err.h"
#include "serror.h"
#include "radar.h"

/*----------------------------------------------------------------------*/
/* External (Import) Variables                                          */
/*----------------------------------------------------------------------*/
extern struct radar_struct *rad;

/*----------------------------------------------------------------------*/
/* Main Function                                                        */
/*----------------------------------------------------------------------*/
int copy_buff(  unsigned char *inbuf, 		/* Source buffer	*/
		int inbufsize, 			/* Source buffer size	*/
		unsigned char *outbuf, 		/* Dest. buffer		*/
		int outbufsize)			/* Dest. buffer	size	*/
{
	if (inbufsize > outbufsize)
		return seterr(ARG_ERR,"in buffer larger than out buffer");

	memcpy(outbuf,inbuf,inbufsize);

	return inbufsize;
}
/*-END OF MODULE--------------------------------------------------------*/
