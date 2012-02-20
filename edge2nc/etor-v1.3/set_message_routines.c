/*----------------------------------------------------------------------**
**
** set_message_routines.c - Set the pointers to the pre and post message
**			    routines.
**
**----------------------------------------------------------------------**
**
** DESCRIPTION
**
** The error handling routines contain two pointers pre_message and 
** post_message they are initially zero but can be set to point to 
** used defined routines to be called before and after error messages
** are printed.
**
**
** USAGE:
**
** void set_message_routines( 
** 	void (*pre)(void), 	-* Function to call before errors printed *-
** 	void (*post)(int, char*)-* Function to call after errors printed  *-
** )				-* Passed the error level and message     *-
**
** PROCESSING:
**
** This function copies the paramters (pre and post) into the global
** variables pre_message and post_message
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
** Modification Date		- $Date: 1999/03/15 16:26:03 $
** Modified by			- $Author: merritt $
** $Source: /nfs/trmm/src/CVS/etor/set_message_routines.c,v $
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

/*----------------------------------------------------------------------*/
/* Application Headers                                                  */
/*----------------------------------------------------------------------*/
#include "serror.h"

/*----------------------------------------------------------------------*/
/* External (Import) Variables                                          */
/*----------------------------------------------------------------------*/
extern char *edge_errlist[];
extern int edge_errno;

/*----------------------------------------------------------------------*/
/* Global (Export) Variables                                            */
/*----------------------------------------------------------------------*/
void (*pre_message)(void) = (void *)0;
void (*post_message)(int, char *) = (void *)0;

/*----------------------------------------------------------------------*/
/* Main Function                                                        */
/*----------------------------------------------------------------------*/

void set_message_routines( 
	void (*pre)(void), 	 /* Function to call before errors printed */
	void (*post)(int, char *)/* Function to call after errors printed  */
)				 /* Passed the error level and message     */
{
#ifdef DUMP_ARGS
	fprintf(stderr, "set_message_routines(pre:0x%08x, post:0x%08x )\n",
		pre, post);
#endif /* DUMP_ARGS */

	pre_message = pre;
	post_message = post;
}

/*-END OF MODULE--------------------------------------------------------*/
