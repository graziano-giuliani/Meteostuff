/*----------------------------------------------------------------------**
**
** curr_time.c - Convert time to string
**
**----------------------------------------------------------------------**
**
** DESCRIPTION
**
** Either convert time in seconds since the epoch to a string, or
** return the current time as a string. 
**
**
** USAGE:
**
** char *curr_time(int argc, char *argv[])
**
** Example IDL call:
** 
** res=CALL_EXTERNAL('libedge.so','curr_time',0L,/S_VALUE)
** res will be a string containing the current time and date
**
** PROCESSING:
**
** This function is called from IDL and so uses the argc argv parameter passing
** convention. One parameter should be used, argv[0] should be a pointer to a 
** a value of type time_t.
** 
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
**	Package			- RTD
**	Reference number	- SP1/PGM/
**	Revision number		- $Revision: 1.1.1.1 $
**	Release State		- $State: Exp $
**	Author, designer	- Bob Stafford
** Modification Date		- $Date: 1999/03/15 16:26:02 $
** Modified by			- $Author: merritt $
** $Source: /nfs/trmm/src/CVS/etor/curr_time.c,v $
**
** MODIFICATION RECORD
**
** $Log: curr_time.c,v $
** Revision 1.1.1.1  1999/03/15 16:26:02  merritt
** Original import.
**
 * Revision 1.1  1997/02/14  15:53:09  bobstaff
 * Post Pakistan Version
 *
 * Revision 1.1  1994/08/11  00:18:47  stafford
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
#include	<stdio.h>	/* stdio library			*/
#include	<stdlib.h>	/* Some standard funct.			*/
#include	<sys/times.h>	/* time() structure			*/
#include	<time.h>	/* Time-of-day functs			*/


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
struct stng
{
	unsigned short len;
	short res;
	char *str;
};

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

char *curr_time(int argc, char *argv[])

{
	time_t t;
	char *ptr;

#ifdef DUMP_ARGS
	printf("curr_time(), time = %x %x %x\n",argv[0],
		*(time_t *)(argv[0]),time(NULL));
#endif /* DUMP_ARGS */


	t = *(time_t *)(argv[0]);

	if (t==(time_t) 1)
	{
		/*
		** A time of 1 second indicates that no status packet has been
		** received recently
		*/
		return("No Status Packet Received");
	}
	else if (t == (time_t)0)
	{
		/*
		** A time of zero indicates a request for the current time
		*/
		t = time(NULL);
	}

	/*
	** Convert time to string, append a null character and return a pointer
	** to the string
	*/
	if (argc>=2)
	{
		struct stng *idl_str;
		static char buff[32];
		int i;
		char *format;


		if (argc>=3)
		{
			t+=*(int *)(argv[2]);
		}
		if (argc>=4 && (*(int *)(argv[3]))!=0)
		{
			t-=t%*(int *)(argv[3]);
		}
		idl_str = (struct stng *)argv[1];
		format = idl_str->str;
		i=strftime(buff,32,format,localtime(&t));
		return buff;
	}
	else
		return (ptr=ctime(&t),ptr[24]=0,ptr);
}
/*-END OF MODULE--------------------------------------------------------*/
