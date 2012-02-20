/*----------------------------------------------------------------------**
**
** save_data.c -- save data in standard EDGE format.
**
**----------------------------------------------------------------------**
**
** DESCRIPTION
**
** Save data in standard EDGE format.
**
**
**
** USAGE:
**
** int save_data(	
**	FILE *data_file, 	-* File to read from		*-
**	unsigned char *data, 	-* Data	to output		*-
**	unsigned int size, 	-* Size of data			*-
**	int compress,		-* Compression code		*-
**	int struct_size)	-* Size of header structure	*-
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
**	Package			- util
**	Reference number	- SP1/PGM/
**	Revision number		- $Revision: 1.1.1.1 $
**	Release State		- $State: Exp $
**	Author, designer	- Bob Stafford
** Modification Date		- $Date: 1999/03/15 16:26:02 $
** Modified by			- $Author: merritt $
** $Source: /nfs/trmm/src/CVS/etor/save_data.c,v $
**
** MODIFICATION RECORD
**
** $Log: save_data.c,v $
** Revision 1.1.1.1  1999/03/15 16:26:02  merritt
** Original import.
**
 * Revision 1.2  1997/06/03  13:07:27  bobstaff
 * Changed copyright message
 *
 * Revision 1.1  1995/12/04  16:56:21  edge
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


/*----------------------------------------------------------------------*/
/* Application Headers                                                  */
/*----------------------------------------------------------------------*/
#include "serror.h"
#include "vol.h"
#include "radar.h"

/*----------------------------------------------------------------------*/
/* Macros                                                               */
/*----------------------------------------------------------------------*/
#define BLOCKSIZE 131072
#define BUFFSIZE (BLOCKSIZE * 2)

/*----------------------------------------------------------------------*/
/* External (Import) Variables                                          */
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* External Functions                                                   */
/*----------------------------------------------------------------------*/
extern int copy_buff(unsigned char *,int,unsigned char*,int);
extern int lzw_compress(unsigned char *,int,unsigned char*,int);
extern int run_len_enc(unsigned char *,int,unsigned char*,int);

/*----------------------------------------------------------------------*/
/* Structures and Unions                                       	        */
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* Global (Export) Variables                                            */
/*----------------------------------------------------------------------*/
char *comp[]={
"NONE",
"LZW",
"RUN LEN ENC",
NULL
};

int (*comp_func[])(unsigned char *,int,unsigned char *,int)=
{
	copy_buff,
	lzw_compress,
	run_len_enc,
	NULL
};

/*----------------------------------------------------------------------*/
/* Main Function                                                        */
/*----------------------------------------------------------------------*/

int save_data(	FILE *data_file, 	/* File to read from		*/
		unsigned char *data, 	/* Data	to output		*/
		unsigned int size, 	/* Size of data			*/
		int compress,		/* Compression code		*/
		int struct_size)	/* Size of header structure	*/
{
	int i,j;
	unsigned char *ptr;
	unsigned char tempbuf[BUFFSIZE];
	char blk_hdr[32];
	int csize;
	int this_block;

#ifdef CHECK_ARGS
#endif /* CHECK_ARGS */
#ifdef DUMP_ARGS
	printf("save_data(size:%d, compress: %d, struct_size : %d)\n",
		size,compress,struct_size);
#endif /* DUMP_ARGS */

	this_block = struct_size?struct_size:BLOCKSIZE;

	set_lzw_verbose(0);
	ptr = (unsigned char *)data;
	for (j=0,i=size;i>0;j++)
	{
		if ((csize=comp_func[compress](ptr,this_block,tempbuf,BUFFSIZE)) == ERROR)
		{
			return error("save_data, compfunc[]()");
		}
		sprintf(blk_hdr,"%d %d %08x",this_block,
			csize,crc(tempbuf,csize));
		if (fwrite(blk_hdr,32,1,data_file)!=1)
			return error("save_data(), fwrite");

		if (fwrite(tempbuf,csize,1,data_file)!=1)
			return error("save_data(), fwrite");
		ptr += this_block;
		i-=this_block;
		this_block = i>BLOCKSIZE?BLOCKSIZE:i;
	}

	return OK;
}
/*-END OF MODULE--------------------------------------------------------*/
