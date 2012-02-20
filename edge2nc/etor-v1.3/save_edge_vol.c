/*----------------------------------------------------------------------**
**
** save_vol.c - Save volume to disk file
**
**----------------------------------------------------------------------**
**
** DESCRIPTION
**
** Save the specified volume to a file with the specified filename using
** the specified compression technique
**
**
**
** USAGE:
**
** int save_vol(
** struct vol_struct *vol, 	-- pointer to volume to save
** char *filename, 		-- volume filename
** int compress			-- compression technique
** 				-- must be one of following macros from vol.h
**				-- NO_COMPRESS
**				-- LZW
**				-- RUN_LEN_ENC
** )
**
** PROCESSING:
**
** This function saves the specified volume to a temporary file on the local
** machine, it is then rcp'ed to each host in the volume destination list.
** If the special destination 'Archive is specified in that list, the volume 
** is also copied into the /edge/data/arch directory on the machine configured
** as the archival host
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
**
** HEADER INFOMATION
**
**	Software Suite 		- EDGE
**	Package			- DAP
**	Reference number	- SP1/PGM/
**	Revision number		- $Revision: 1.1.1.1 $
**	Release State		- $State: Exp $
**	Author, designer	- Bob Stafford
** Modification Date		- $Date: 1999/03/15 16:26:03 $
** Modified by			- $Author: merritt $
** $Source: /nfs/trmm/src/CVS/etor/save_edge_vol.c,v $
**
** MODIFICATION RECORD
**
** $Log: save_edge_vol.c,v $
** Revision 1.1.1.1  1999/03/15 16:26:03  merritt
** Original import.
**
 * Revision 1.8  1997/05/13  13:44:23  bobstaff
 * Use the volume distribution code for disk and tape archival
 *
 * Revision 1.7  1997/03/19  18:12:17  bobstaff
 * Always save the volume if collecting a clutter map
 *
 * Revision 1.6  1997/03/04  15:00:53  bobstaff
 * Dealt with network distribution of volumes, first version
 *
 * Revision 1.5  1997/02/25  19:00:23  kenb
 * Added support for Zdr
 *
 * Revision 1.4  1997/02/07  15:09:40  edge
 * *** empty log message ***
 *
 * Revision 1.3  1995/12/04  16:41:33  edge
 * Pre-Thailand version
 *
 * Revision 1.2  1994/08/09  20:42:13  stafford
 * Change LOG to Log in RCS section of Header so that log appears in file
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
#include	<stdarg.h>	/* Some standard funct.			*/
#include 	<errno.h>	/* All error codes			*/
#include	<fcntl.h>	/* File control functs			*/
#include	<unistd.h>
#include	<string.h>
#include	<time.h>
#include	<sys/param.h>
#include	<sys/types.h>
#include	<sys/stat.h>


/*----------------------------------------------------------------------*/
/* Application Headers                                                  */
/*----------------------------------------------------------------------*/
#include "vol.h"
/* #include "radar.h" */

/*----------------------------------------------------------------------*/
/* Macros                                                               */
/*----------------------------------------------------------------------*/
#define BLOCKSIZE 1024
#define CHECK_ARGS
#define DUMP_ARGS
#define OK 0
#define ERROR -1

/*----------------------------------------------------------------------*/
/* Global (Import) Variables                                            */
/*----------------------------------------------------------------------*/
extern char *comp[];

/*----------------------------------------------------------------------*/
/* Global (Export) Variables                                            */
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* Local (Static) Variables                                             */
/*----------------------------------------------------------------------*/
struct idl_string
{
	unsigned short length;
	short reserved;
	char *s;
};

/*----------------------------------------------------------------------*/
/* Local Function                                                       */
/*----------------------------------------------------------------------*/
/*
** Function to save the file header. This header is a 20 line ASCII section
** containing information about the contents of the file. All edge data files
** should start with such a header
** It may be displayed using the Unix command 
** head -20 <filename>
*/
static int save_header(FILE *vol_file, struct vol_struct *vol, 
	char *filename, int compress)
{
	int i;
	char *pw[]={"0.8","2.0","?","?"};

	fprintf(vol_file,"Volume File : %-40s\n",filename);
	printf("Volume File : %-40s\n",filename);
	fprintf(vol_file,"Compressed  : %-10s\n",comp[compress]);
	printf("Compressed  : %-10s\n",comp[compress]);
	fprintf(vol_file,"Size        : %08x\n",vol->mem_used);
	printf("Size        : %08x\n",vol->mem_used);
	fprintf(vol_file,"Checksum    : %08x\n",crc((unsigned char *)vol,vol->mem_used));
	printf("Checksum    : %08x\n",crc((unsigned char *)vol,vol->mem_used));
        time_t xtime = vol->date;
	fprintf(vol_file,"Date        : %-24s",asctime(localtime(&xtime)));
	printf("Date        : %-24s",asctime(localtime(&xtime)));
	fprintf(vol_file,"Moments     : %c%c%c%c%c\n",
		(vol->moment_enable|ARC_MOMENT)?'A':' ',
		(vol->moment_enable|U_MOMENT)?'U':' ',
		(vol->moment_enable|Z_MOMENT)?'Z':' ',
		(vol->moment_enable|V_MOMENT)?'V':' ',
		(vol->moment_enable|W_MOMENT)?'W':' ',
		(vol->moment_enable|ZDR_MOMENT)?'D':' ');
	printf("Moments     : %c%c%c%c%c\n",
		(vol->moment_enable|ARC_MOMENT)?'A':' ',
		(vol->moment_enable|U_MOMENT)?'U':' ',
		(vol->moment_enable|Z_MOMENT)?'Z':' ',
		(vol->moment_enable|V_MOMENT)?'V':' ',
		(vol->moment_enable|W_MOMENT)?'W':' ',
		(vol->moment_enable|ZDR_MOMENT)?'D':' ');
	fprintf(vol_file,"Scan Type   : %-10s Rng: %d PW:%s\n", 
		SCAN_TYPE(vol->scan_type),
		(int)(vol->sweep[0].max_range),
		pw[vol->sweep[0].rad.pulse_width]);
	printf("Scan Type   : %-10s Rng: %d PW:%s\n", 
		SCAN_TYPE(vol->scan_type),
		(int)(vol->sweep[0].max_range),
		pw[vol->sweep[0].rad.pulse_width]);
	fprintf(vol_file,"Sweeps      : %02d\n", vol->num_sweeps);
	printf("Sweeps      : %02d\n", vol->num_sweeps);
	fprintf(vol_file,"Elevations  :");
	printf("Elevations  :");
	for(i=0;i<vol->num_sweeps;i++)
	{
		fprintf(vol_file," %5.2f",DEGS(vol->sweep[i].rad.el));
		printf(" %5.2f",DEGS(vol->sweep[i].rad.el));
	}
	fprintf(vol_file,"\n");
	printf("\n");
	fprintf(vol_file,"Version     : %5.2f\n", (float)(vol->version)/100.0);
	printf("Version     : %5.2f\n", (float)(vol->version)/100.0);
	fprintf(vol_file,"Proc File   : %s\n", vol->sweep[0].rad.job_name);
	printf("Proc File   : %s\n", vol->sweep[0].rad.job_name);
		
	/*
	** Pad header to 20 lines total
	*/
	for(i=0;i<9;i++)
		fprintf(vol_file,"\n");

	printf("saved header\n");	
	return OK;
}
/*----------------------------------------------------------------------*/
/* Main Function                                                        */
/*----------------------------------------------------------------------*/


int save_edge_vol(
	struct vol_struct *vol, 
	char *filename, 
	int compress
)
{

	FILE	*vol_file;
	char *temp_file="/tmp/vol.vol";
	int i;

	temp_file = filename;

	/*
	** If volume is NULL then return with an error condition
	*/
	if (vol == NULL) 
	{
		fprintf(stderr,"vol == NULL"); 
		return ERROR;
	}
#ifdef CHECK_ARGS
	if (compress<0 || compress>MAX_COMPRESS)
	{
		fprintf(stderr, "compress == %d",compress);
		return ERROR;
	}
#endif /* CHECK_ARGS */

	printf("save_vol() filename = %s compress = %d (%s)\n",filename,
		compress,comp[compress]);

	/*
	** Open volume fime for writing
	*/
	printf("opening %s\n",temp_file);
	if ((vol_file = fopen(temp_file,"w")) == NULL)
	{
		fprintf(stderr, "save_vol(), fopen: %s",temp_file);
		return ERROR;
	}

	/*
	** Save the ascii header at start of file
	*/
	printf("saving header %lx\n",vol_file);
	if (save_header(vol_file, vol, filename , compress) == ERROR)
	{
		fclose(vol_file);
		fprintf(stderr, "save_vol(), save_header");
		return ERROR;
	}

	/*
	** Save the volume data
	*/
	printf("saving data\n");
	if (save_data(vol_file,(unsigned char *)vol,vol->mem_used,
		compress,sizeof(struct vol_struct))==ERROR)
	{
		fclose(vol_file);
		fprintf(stderr, "save_vol(), save_data");
		return ERROR;
	}

	/*
	** Close the file
	*/
	fclose(vol_file);

	printf("volume file saved in temp file\n");
	return OK;

}
