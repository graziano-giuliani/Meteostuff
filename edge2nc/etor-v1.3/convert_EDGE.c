/*----------------------------------------------------------------------**
**
** convert_EDGE.c - Program to convert an EDGE file to RSL radar 
** structure save the structure to disk.
**
**----------------------------------------------------------------------**
**
** DESCRIPTION
**
** This program is meant to provide an interface between the EDGE file format
** and the TRMM Radar Software Library. EDGE volume files would then be 
** accessible to the radar meteorology/hydrology research community worldwide.
**
** USAGE:
**
** int convert_EDGE(char *infile);	-* Output files are named using the -*
**					-* infile name but converting the   -*
**					-* .vol extension to .rad for the   -*
**					-* radar structure and .uf for the  -*
**					-* universal format file.           -*
**
** PROCESSING:
**
** The program reads in an EDGE volume file converts it into an RSL 
** radar structure. The structure is saved to disk. The structure is
** also converted to universal format and saved to disk as a compressed
** universal format file.
**
** COPYRIGHT NOTICE
**
**	Copyright (c) 1997 by Enterprise Electronics Corporation
**	All Rights Reserved
** 
** This program is  copyright  by  Enterprise  Electronics  Corpora-
** tion,    Enterprise,  Alabama,  USA  36330 (334) 347-3478.  It is
** licensed for  use  on  a  specific  CPU   and   is  not    to  be
** transferred  or otherwise divulged.   Copies  or modifications of
** this program must carry this copyright notice.
** 
**
**
** HEADER INFOMATION
**
**	Software Suite 		- EDGE
**	Package			- TRMM RSL Interface
**	Reference number	- SP1/PGM/
**	Revision number		- $Revision: 1.1.1.1 $
**	Release State		- $State: Exp $
**	Author, designer	- Don Burrows
** Modification Date		- $Date: 1999/03/15 16:26:02 $
** Modified by			- $Author: merritt $
** $Source: /nfs/trmm/src/CVS/etor/convert_EDGE.c,v $
**
** MODIFICATION RECORD
**
** $Log: convert_EDGE.c,v $
** Revision 1.1.1.1  1999/03/15 16:26:02  merritt
** Original import.
**
**
**----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* Feature Test Switches                                                */
/*----------------------------------------------------------------------*/
#define _POSIX_SOURCE	1

/*----------------------------------------------------------------------*/
/* System Headers            { full list in stdinc.h }                  */
/*----------------------------------------------------------------------*/
#include	<stdio.h>	/* stdio library			*/
#include	<stddef.h>	/* Some popular symbols			*/
#include	<stdlib.h>	/* Some standard funct.			*/
#include	<unistd.h>	/* POSIX symbols definitions		*/
#include	<string.h>	/* String functions			*/

/*----------------------------------------------------------------------*/
/* Application Headers                                                  */
/*----------------------------------------------------------------------*/
#include "rsl.h"

/*----------------------------------------------------------------------*/
/* Macros                                                               */
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* External (Import) Variables                                          */
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* External Functions                                                   */
/*----------------------------------------------------------------------*/
Radar *RSL_EDGE_to_radar(char *infile);

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

int convert_EDGE(char *infile)	
{
	Radar *rad;
	int dot;
	char *ext;
	char outfile_rad[255],outfile_uf[255];

	printf("Infile: %s\n",infile);
	rad = RSL_EDGE_to_radar(infile);
	printf("Edge converted to radar structure.\n");
	ext = strrchr(infile,'.');
	*ext = 0x00;
	strcpy(outfile_rad,infile);
	strcat(outfile_rad,".rad");
	strcpy(outfile_uf,infile);
	strcat(outfile_uf,".uf");
	printf("Outfile is %s\n",outfile_uf);
	RSL_write_radar_gzip(rad,outfile_rad);
	printf("Radar file written\n");
	RSL_radar_to_uf_gzip(rad,outfile_uf);
	printf("Universal format file written.\n");
	return 0;
}

/*----------------------------------------------------------------------*/
/* Main Function                                                        */
/*----------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
	char *infile;

	if (argc < 2) 
	{
		fprintf(stderr,"convert_EDGE: No input file specified.\n");
		exit(1);
	}
	infile = argv[1];
	exit(convert_EDGE(infile));
}

/*-END OF MODULE--------------------------------------------------------*/
