/*----------------------------------------------------------------------**
**
** culprit.c - Set error info pointer
**
**----------------------------------------------------------------------**
**
** DESCRIPTION
**
** Copy parameter passed to function into the global variable error_info
** which is printed as added information after error messages.  For example
**
**	if(open(culprit("nofile.fil"),0) < 0)
**		return error( "func(), open");
**	return OK;
**
** Could produce the following error message:
**
** SYS:*:20Sep93:18:06:23:err/err.c(31):No such file or directory:nofile.fil
**
** The "nofile.fil" part was supplied by culprit in the open.  Returning its
** argument makes culprit a bit more usable.
**
**
** USAGE:
**
** char *culprit( error_aux )
** char *error_aux			-* Auxillary error information *-
**
** PROCESSING:
**
** This routine simply copies the value passed to it into the global variable
** error_info which is a pointer to a character declared in serror.c.
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
**	Package			- ERR
**	Reference number	- SP1/PGM/
**	Revision number		- $Revision: 1.1.1.1 $
**	Release State		- $State: Exp $
**	Author, designer	- KenB
** Modification Date		- $Date: 1999/03/15 16:26:02 $
** Modified by			- $Author: merritt $
** $Source: /nfs/trmm/src/CVS/etor/culprit.c,v $
**
** MODIFICATION RECORD
**
** $LOG:$
**
**----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* Feature Test Switches                                                */
/*----------------------------------------------------------------------*/
#define _POSIX_SOURCE	1

/*----------------------------------------------------------------------*/
/* System Headers                                                       */
/*----------------------------------------------------------------------*/
#include	<stdio.h>	/* stdio library			*/
#include	<string.h>

/*----------------------------------------------------------------------*/
/* External (Import) Variables                                          */
/*----------------------------------------------------------------------*/
extern char *error_info;	/* Printed as extra info with errors	*/
static char error_buff[128];

/*----------------------------------------------------------------------*/
/* Main Function                                                        */
/*----------------------------------------------------------------------*/

char *culprit( 
char *error_aux		/* Some auxillary error information	*/
)
{
#ifdef DUMP_ARGS
	fprintf(stderr, "name(error_aux:'%s')\n", error_aux);
#endif /* DUMP_ARGS */

	if (error_aux)
	{
		strncpy(error_buff,error_aux,128);
		return error_info = error_buff;
	}
	else
	{
		return error_info = error_aux;
	}
}

/*-END OF MODULE--------------------------------------------------------*/



