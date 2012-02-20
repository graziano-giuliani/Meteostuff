/*----------------------------------------------------------------------**
**
**int get_path(char *filename, char *pathname) - Separate the pathname from the complete filename
**
**----------------------------------------------------------------------**
**
** DESCRIPTION
**
**
**
** USAGE:
**
** int get_path(char * filename,char *pathname)
** (char *filename	- Filename including complete path
** char *pathname)	- Path separated from filename
**
** PROCESSING:
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
** HEADER INFOMATION
**
**	Software Suite 		- EDGE
**	Package			- UIF
**	Reference number	- SP1/PGM/
**	Revision number		- $Revision: 1.1.1.1 $
**	Release State		- $State: Exp $
**	Author, designer	- Don Burrows
**
** MODIFICATION RECORD
**
**
**
**----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* Feature Test Switches                                                */
/*----------------------------------------------------------------------*/
#define _POSIX_SOURCE	1

/*----------------------------------------------------------------------*/
/* System Headers                                                       */
/*----------------------------------------------------------------------*/
#include	<assert.h>	/* assert() macro			*/
#include 	<ctype.h>	/* Character test fncts			*/
#include	<dirent.h>	/* Dir. entry module			*/
#include 	<errno.h>	/* All error codes			*/
#include	<fcntl.h>	/* File control functs			*/
#include	<float.h>	/* Symbols for floating point		*/
#include	<grp.h>		/* Group database functions		*/
#include 	<limits.h>	/* Implementation limits		*/
#include	<locale.h>	/* Multi-national apps. symbols		*/
#include	<math.h>	/* Standard math functions		*/
#include 	<pwd.h>		/* User database functions		*/
#include 	<setjmp.h>	/* setjmp() / longjmp() macros		*/
#include	<signal.h>	/* Signal's symbols and functs		*/
#include 	<stdarg.h>	/* Var. argument funct. support		*/
#include	<stddef.h>	/* Some popular symbols			*/
#include	<stdio.h>	/* stdio library			*/
#include	<stdlib.h>	/* Some standard funct.			*/
#include	<string.h>	/* String functions			*/
#include 	<sys/stat.h>	/* File stat structure			*/
#include	<sys/times.h>	/* time() structure			*/
#include	<sys/types.h>	/* POSIX data types			*/
#include	<sys/utsname.h> /* uname() structure			*/
#include	<sys/wait.h>	/* wait() functions			*/
#include	<termios.h>	/* Term. manipulation			*/
#include	<time.h>	/* Time-of-day functs			*/
#include	<unistd.h>	/* POSIX symbols definitions		*/
#include 	<utime.h>	/* utime() structure			*/


/*----------------------------------------------------------------------*/
/* Application Headers                                                  */
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* Macros                                                               */
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* External (Import) Variables                                          */
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* External Functions                                                   */
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* Structures and Unions                                       	        */
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* Global (Export) Variables                                            */
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* Local (Static) Variables                                             */
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* Signal Catching Functions                                            */
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* Local Function                                                       */
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* Main Function                                                        */
/*----------------------------------------------------------------------*/

int get_path(const char *filename,char *pathname)
{
	char *name,lastname;
	int filelen,namelen;

	if (pathname != NULL)
	{
		filelen = strlen(filename);
		name = strrchr(filename,'/');
		if (name == '\0') 
		{
			printf("Data path is NULL!\n");
			pathname[0]='\0';
			return 0;
		}
		else
		{
			name++;
		}
		namelen = strlen(name);
		filelen -= namelen;
		strncpy(pathname,filename,filelen);
		pathname[filelen]='\0';
		printf("Data path is %s\n",pathname);
	}
	else 
	{
		fprintf(stderr,"No Data Path storage was provided. pathname is NULL!\n");
		return -1;
	}
	return 0;
}

/*-END OF MODULE--------------------------------------------------------*/
