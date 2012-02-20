/*----------------------------------------------------------------------**
**
** load_data.c - Load a standard EDGE file
**
**----------------------------------------------------------------------**
**
** DESCRIPTION
**
** Load a standard EDGE file
**
** USAGE:
**
** int load_data(
** 	char **data, 			-* Point to read data pinter	*-
** 	char *filename, 		-* File to read			*-
** 	unsigned int filetype )		-* Sort of file, see defs	*-
**
** PROCESSING:
**
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
**	Revision number		- $Revision: 1.3 $
**	Release State		- $State: Exp $
**	Author, designer	- Bob Stafford
** Modification Date		- $Date: 1999/04/02 16:22:10 $
** Modified by			- $Author: merritt $
** $Source: /nfs/trmm/src/CVS/etor/load_data.c,v $
**
** MODIFICATION RECORD
**
** $Log: load_data.c,v $
** Revision 1.3  1999/04/02 16:22:10  merritt
** ready for v1_0
**
** Revision 1.2  1999/03/31 22:35:29  merritt
** the etor lib
**
** Revision 1.1.1.1  1999/03/15 16:26:02  merritt
** Original import.
**
 * Revision 1.4  1997/06/03  13:02:38  bobstaff
 * Updated header
 * include dap.h instead of dap.H
 *
 * Revision 1.3  1997/02/14  15:51:12  bobstaff
 * Post Pakistan Version
 *
 * Revision 1.2  1996/03/15  17:34:11  edge
 * KenB Added global 'compression' to indicated method
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
#include	<string.h>	/* String functions			*/
/*----------------------------------------------------------------------*/
/* Application Headers                                                  */
/*----------------------------------------------------------------------*/
#include "err.h"
#include "serror.h"
#include "vol.h"
#include "util.h"
#include "radar.h"

/*----------------------------------------------------------------------*/
/* Macros                                                               */
/*----------------------------------------------------------------------*/
#define BUFFSIZE 65536

/*----------------------------------------------------------------------*/
/* External (Import) Variables                                          */
/*----------------------------------------------------------------------*/
extern char *comp[];

/*----------------------------------------------------------------------*/
/* External Functions                                                   */
/*----------------------------------------------------------------------*/
#ifndef _UTIL_H_
extern int copy_buff(unsigned char *,int,unsigned char*,int);
extern int lzw_uncompress(unsigned char *,int,unsigned char*,int);
extern int run_len_dec(unsigned char *,int,unsigned char*,int);
#endif

/*----------------------------------------------------------------------*/
/* Global (Export) Variables                                            */
/*----------------------------------------------------------------------*/

int (*uncomp_func[])(unsigned char *,int,unsigned char *,int)=
{
	copy_buff,
	lzw_uncompress,
	run_len_dec,
	NULL
};
int compression;

/*----------------------------------------------------------------------*/
/* Local (Static) Variables                                             */
/*----------------------------------------------------------------------*/
char *env_sym[] = {
	"EDGE_VOL_DIR",
	"EDGE_PROD_DIR",
	NULL
};

/*----------------------------------------------------------------------*/
/* Local Function                                                       */
/*----------------------------------------------------------------------*/
static char temp[256];
unsigned char *get_header()
{
	return (unsigned char *)temp;
}

static int read_header(FILE *data_file, int *compress, int *size, int *checksum)
{
	int i;

	fscanf(data_file,"%*[^\n]\n");
	if (fscanf(data_file,"Compressed  : %10s\n",temp) != 1)
	{
		return error("fscanf failed to read compression");
	}
	fscanf(data_file,"Size        : %08x\n",size);
	fscanf(data_file,"Checksum    : %08x\n",checksum);

	*compress = -1;

	for (i=0;comp[i]!=NULL;i++)
	{
		if (strcmp(comp[i],temp)==0)
		{
			*compress = i;
			break;
		}
	}

	if (*compress == -1)
	{
		return seterr(VOL_READ,"Unknown compression %s",temp);
	}
	/*
	** Read rest of header
	*/
	for(i=0;i<16;i++)
	{
		fgets(temp,255,data_file);
	}
	
	return OK;
}
/*----------------------------------------------------------------------*/
/* Main Function                                                        */
/*----------------------------------------------------------------------*/

int load_data(
	char **data, 			/* Point to read data pinter	*/
	char *filename, 		/* File to read			*/
	unsigned int filetype )		/* Sort of file, see defs	*/
{

	FILE	*data_file;
	int i,j,k;
	unsigned char *ptr;
	unsigned char *tempbuf;
	unsigned int tempbuf_size=BUFFSIZE;
	char blk_hdr[32];
	int temp;
	int compress;
	int size,csize,usize,asize;
	int checksum,bcrc;
	char path[256];
	char *path_ptr;
#ifdef CHECK_ARGS
	if (*data != NULL)
		return seterr(ARG_ERR,"*data != NULL");
	if (filename == NULL)
		return seterr(ARG_ERR,"filename == NULL");
	if ((filetype & 0xff) > PROD_FILE && filetype != -1)
		return seterr(ARG_ERR,"filetype == %d",filetype & 0xff);
#endif /* CHECK_ARGS */
#ifdef DUMP_ARGS
	printf("load_data() filename = %s\n",filename);
	printf("load_data() filetype = %d\n",filetype & 0xff);
#endif /* DUMP_ARGS */

	if ((tempbuf=malloc(tempbuf_size))==NULL)
	{
		return seterr(NO_MEM,"Unable to allocate memory for data");
	}

	if (*data != NULL)
	{
		free(*data);
		*data = NULL;
	}
	if ((data_file = fopen(culprit(filename),"r")) == NULL)
	{
		if (filetype != -1 )
		{
			/*
			** Get default path for data file type
			*/
#ifdef DEBUG
			printf("calling getenv %s\n",env_sym[filetype & 0xff]);
#endif /* DEBUG */
			path_ptr = getenv(env_sym[filetype & 0xff]);
			if (path_ptr == NULL)
			{
				/*
				** No PATH set up so give up
				*/
				free(tempbuf);
				return error("load_data(), Unable to open "
					"file %s",filename);
			}
#ifdef DEBUG
			printf("path = %s\n",path_ptr);
#endif /* DEBUG */

			/*
			** copy path so as not to corrupt enviroment symbol table
			** when calling strcat in a minute (nano-second ?)
			*/
			strcpy(path,path_ptr);
			strcat(path,"/");
			if ((data_file = fopen(culprit(strcat(path,filename)),
				"r")) == NULL)
			{
				if ((data_file = fopen(culprit(strcat(path,".vol")),"r")) == NULL)
				{
			
					free(tempbuf);
					return error("load_data(), fopen");
				}
			}
		}
		else
		{
			free(tempbuf);
			return error("load_data(), fopen");
		}
	}
	
#ifdef DEBUG
	printf("load_data: reading header\n");
#endif /* DEBUG */
	if (read_header(data_file, &compress, &size, &checksum) == ERROR)
	{
		free(tempbuf);
		fclose(data_file);
		return error("load_data(), read_header");
	}
	compression = compress;	/* Set global in case anyone is interested*/

#ifdef DEBUG
	printf("load_data: reading data\n");
#endif /* DEBUG */
	j=0;
	while (!feof(data_file))
	{
		if (fread(blk_hdr,32,1,data_file)!=1)
		{
			break;		/* End of file */
		}
		if (sscanf(blk_hdr,"%d %d %08x",&usize,&csize,&bcrc) == EOF)
		{
			if (*data != NULL)
			{
				free(*data);
				*data=NULL;
			}
			free(tempbuf);
			fclose(data_file);
			return error("load_data(), sscanf");
		}
#ifdef DEBUG
		printf("usize: %d csize : %d crc : %08x\n",usize,csize,bcrc);
#endif /* DEBUG */

		if (csize > tempbuf_size)
		{
			free(tempbuf);
			tempbuf_size=csize;
			if ((tempbuf=malloc(tempbuf_size))==NULL)
			{
				if (*data != NULL)
				{
					free(*data);
					*data=NULL;
				}
				fclose(data_file);
				return seterr(NO_MEM,
					"Unable to allocate memory for data");
			}
		}
		if (j == 0)
		{
			if (filetype & HEADER_ONLY)
			{
#ifdef DEBUG
				printf("header_only\n");
#endif /* DEBUG */
				asize=usize;
			}
			else if (filetype & PARTIAL_VOL)
			{
				asize = NEW_VOL_SIZE;
			}
			else
			{	
				asize = size;
			}
			
#ifdef DEBUG
			printf("mallocing %d bytes\n",asize);
#endif /* DEBUG */
			if ((*data=malloc(asize))==NULL)
			{
				free(tempbuf);
				fclose(data_file);
				return seterr(NO_MEM,
					"Unable to allocate memory for data");
			}
			ptr = (unsigned char *)*data;
		}
		
		if (fread(tempbuf,csize,1,data_file)!=1)
		{
			free(tempbuf);
			fclose(data_file);
			if (*data != NULL)
			{
				free(*data);
				*data=NULL;
			}
			return error("load_data(), fread");
		}
		temp=crc(tempbuf,csize);
		if (temp != bcrc)
		{
			free(tempbuf);
			if (*data != NULL)
			{
				free(*data);
				*data=NULL;
			}
			fclose(data_file);
			return seterr(VOL_READ,"Bad block checksum %08x %08x",
				temp,bcrc);
		}
		if ((temp=uncomp_func[compress](tempbuf,csize,ptr,usize)) == 
		  ERROR)
		{
			free(tempbuf);
			if (*data != NULL)
			{
				free(*data);
				*data=NULL;
			}
			fclose(data_file);
			return error("load_data, uncomp_func[]()");
		}
		ptr+=usize;
		j++;
		if (filetype & HEADER_ONLY)
			break;
	}
	free(tempbuf);
	fclose(data_file);

#ifdef DEBUG
	printf("load_data: checking crc\n");
#endif /* DEBUG */

	if ((filetype & HEADER_ONLY) == 0)
	{
		/*
		** Not header only so check the checksum
		*/
		if ((temp=crc((char *)*data,size)) != checksum)
		{
			if (*data != NULL)
			{
				free(*data);
				*data=NULL;
			}
			return seterr(VOL_READ,"Bad checksum %08x %08x",
				checksum,temp);
		}
	}

#ifdef DEBUG	
	printf("load_data: finished\n");
#endif /* DEBUG */

	return OK;
}
/*-END OF MODULE--------------------------------------------------------*/
