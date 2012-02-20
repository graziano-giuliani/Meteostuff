/*----------------------------------------------------------------------**
**
** error.c -- Register an error returned.
**
**----------------------------------------------------------------------**
**
** DESCRIPTION
**
** When calling a routine which returns an error condition, error() is called
** the macro error() calls _error after passing err_is() the __LINE__ 
** and __FILE__ macros.  _error paramters are printf format strings
** followed with the parameters.
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
** and an ERROR indication is returned.  The processing is just like set_err
** except the error number is not set.
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
** Modification Date		- $Date: 1999/03/31 22:35:29 $
** Modified by			- $Author: merritt $
** $Source: /nfs/trmm/src/CVS/etor/error.c,v $
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
extern char *err_src;
extern int err_line;

/*----------------------------------------------------------------------*/
/* Main Function                                                        */
/*----------------------------------------------------------------------*/

int _error(
	char *info,		/* Informational message, in printf fmt */
	...   			/* Arguments to info, printf style	*/
)
{
	va_list args;
	char format[512];	/* Format string 			*/
	int errnum=0;		/* Error number				*/
	char *strdup(const char *s);
#ifdef DUMP_ARGS
	fprintf(stderr, "_error(errnum:0x%08x, info:'%s' )\n", errnum, info);
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
		free(error_stack[errc]);
		error_stack[errc] = strdup(error_message);
	}
	va_end(args);

	/*
	** Clean up, set globals to keep from printing same values twice.
	*/
	err_line = -1;
	err_src = (char *)0;

	/*
	** Always return ERROR, users will "return seterr"
	*/
	return ERROR;
}

/*-END OF MODULE--------------------------------------------------------*/
