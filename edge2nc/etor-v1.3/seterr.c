/*----------------------------------------------------------------------**
**
** seterr.c - set an error condition.
**
**----------------------------------------------------------------------**
**
** DESCRIPTION
**
** The seterr() macro, which calls _seterr, is used to indicate an error.
** Both system (designated SYS in the message) and application (EDGE) errors
** are handled.  If a system error has occured, it is not necessary to 
** set the error_num argument.  For example, if malloc fails:
**
**	if((abc = malloc(numb)) == NULL)
**		return seterr(0,"Malloc failed me: %d bytes", numb);
**
** Assuming numb = 123, a later call to serror may produce:
**
** SYS:*:19Sep94:19:26:46:err/err.c(98):Not enough core:
** 	err/err.c(83):funcy(), Malloc failed me: 123 bytes
** 
** _seterr is always called from the macro seterr() in serror.h.  Seterr
** calls err_id with the language macros __LINE__ and __FILE__.  {Err_id() 
** stores them in global.}  Seterr()'s parameters are then passed to 
** _seterr.  _seterr works like printf, using a printf style string followed
** by a variable number of arguments; the error number is passed first.
** _seterr puts the error message on a stack which is dumped at a call to
** serror().
**
** USAGE:
**
** int _seterr(
** 	int errnum,		-* Error number				*-
** 	char *info,		-* Informational message, in printf fmt *-
** 	...   			-* Arguments to info, printf style	*-
** )
**
** PROCESSING:
**
** The line and file named stored in global by err_id() are printed into
** a message string along with the printf style string and all arguments.
** If the error stack is not full, this message is pushed on.  If it is full, 
** the top element is REPLACED with the message.  The error numbers are reset
** and an ERROR indication is returned.
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
**	Revision number		- $Revision: 1.2 $
**	Release State		- $State: Exp $
**	Author, designer	- KenB
** Modification Date		- $Date: 1999/03/31 22:35:30 $
** Modified by			- $Author: merritt $
** $Source: /nfs/trmm/src/CVS/etor/seterr.c,v $
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
#include 	<stdarg.h>	/* Var. argument funct. support		*/
#include	<stdio.h>	/* stdio library			*/
#include	<string.h>	/* String functions			*/
#include	<time.h>	/* Time-of-day functs			*/

/*----------------------------------------------------------------------*/
/* Application Headers                                                  */
/*----------------------------------------------------------------------*/
#include "serror.h"

#undef DUMP_ARGS
/*----------------------------------------------------------------------*/
/* External (Import) Variables                                          */
/*----------------------------------------------------------------------*/
extern char *error_stack[MAX_ERROR_STACK];
extern int errc;
extern int edge_errno;

extern char error_message[512]; /* Temp storage for message		*/

/*----------------------------------------------------------------------*/
/* External Functions                                                   */
/*----------------------------------------------------------------------*/

/* extern char *strdup( char * ); */

/*----------------------------------------------------------------------*/
/* Global (Export) Variables                                            */
/*----------------------------------------------------------------------*/
char *err_src;
int err_line;

/*----------------------------------------------------------------------*/
/* Main Function                                                        */
/*----------------------------------------------------------------------*/

int _seterr(
	int errnum,		/* Error number				*/
	char *info,		/* Informational message, in printf fmt */
	...   			/* Arguments to info, printf style	*/
)
{
	va_list args;
	char format[512];	/* Format string 			*/
	char *strdup(const char *s);
#ifdef DUMP_ARGS
	fprintf(stderr, "_seterr(errnum:0x%08x, info:'%s' )\n", errnum, info);
#endif  /* DUMP_ARGS */
	va_start(args, info);

	/*
	** Set EDGE error number, if non zero
	*/
	if(errnum != 0)
		edge_errno = errnum;
	/*
	** Place source file name and line from globals into fmt string
	*/
	sprintf(format, "%s(%d):%s", err_src, err_line, info);

	/*
	** Put user arguments into message
	*/
	vsprintf(error_message, format, args);

	/*
	** If error stack is not full, put this message on, otherwise
	** replace the top message with this one.
	*/
	if(errc < MAX_ERROR_STACK)
	{
		error_stack[errc++] = strdup(error_message);
	}
	else
	{
		errc = MAX_ERROR_STACK-1;
		free(error_stack[errc-1]);
		error_stack[errc-1] = strdup(error_message);
	}
	va_end(args);

	/*
	** Always return ERROR, users will "return seterr"
	*/
	return ERROR;
}

/*-END OF MODULE--------------------------------------------------------*/
