/*----------------------------------------------------------------------**
**
** err_id.c -- Identify the line number and source file of an error
**
**----------------------------------------------------------------------**
**
** DESCRIPTION
**
** This will identify the source file and line number of an error.  It is
** called prior to error() or seterr().  The line number and source file
** will be reported in the error message printed by serror().  The file name 
** is trimmed of leading dots and slashes.  i.e.  ../../../err/error_p.c will
** be trimmed to err/error_p.c.  err_id returns the trimmed source name.
**
**
**
** USAGE:
**
** char *err_id(			-* Identify the line and file   *-
** 		int line, 		-* Line number			*-
** 		char *source		-* Source file name		*-
** )
**
** PROCESSING:
**
** The source file name is passed as a string, the line as an integer.  Each
** character, starting with the first, is examined.  If it is a dot ('.') or
** a slash ('/'), the pointer is incremented.  The resulting pointer and line
** number are assigned to global variables, err_src and err_line,
** defined in seterr().
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
**	Reference number
**	Revision number		- $Revision: 1.1.1.1 $
**	Release State		- $State: Exp $
**	Author, designer	- KenB
** 	Modification Date	- $Date: 1999/03/15 16:26:02 $
** 	Modified by		- $Author: merritt $
** 				  $Source: /nfs/trmm/src/CVS/etor/err_id.c,v $
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
#undef DUMP_ARGS
#undef DEBUG

/*----------------------------------------------------------------------*/
/* System Headers                                                       */
/*----------------------------------------------------------------------*/
#include	<stdio.h>	/* stdio library			*/
#include	<string.h>

/*----------------------------------------------------------------------*/
/* External (Import) Variables                                          */
/*----------------------------------------------------------------------*/
volatile extern char *err_src;
volatile extern int err_line;
/*----------------------------------------------------------------------*/
/* Main Function                                                        */
/*----------------------------------------------------------------------*/
char *err_id(				/* Identify the line and file   */
 		int line, 		/* Line number			*/
		char *source		/* Source file name		*/
)
{
#ifdef DUMP_ARGS
	fprintf(stderr, "void err_id(line:%d,source:'%s')\n", line, source);
#endif /* DUMP_ARGS */

#ifdef TRIM_DOTS

	err_src = source;
	while(*err_src && ((*err_src == '.') || (*err_src == '/')))
		err_src++;
	err_line = line;

#else  /* TRIM_DOTS */
	{
		int i, n;
		int slashes;

		n = strlen(source) - 1;
		slashes = 0;
		for(i=n; i>0; i--)
		{
#ifdef DEBUG
			printf("source[%d]=%c\n",i,source[i]);
#endif /*  DEBUG  */
			if(source[i] == '/')
				slashes++;
			if(slashes > 1)
			{
				i++;
				break;
			}
		}
		err_src = source + i;
		err_line = line;
	}
#endif /* TRIM_DOTS */


	return (char *)err_src;
}
/*-END OF MODULE--------------------------------------------------------*/
