/*----------------------------------------------------------------------**
**
** RSL_to_EDGE.c
**
**----------------------------------------------------------------------**
**
** DESCRIPTION
**
** Converts an RSL Radar structure to an EDGE volume structure
**
** USAGE:
**
** int RSL_to_EDGE(RSL_rad, data);	-* Usage Section with comments *-
** Radar RSL_rad;			-* RSL Radar structure pointer *-
** char **data;				-* pointer to EDGE volume structure pointer*-
**
** PROCESSING:
**
**	This program accepts an RSL radar data structure copies it into an
**	Edge volume data structure and returns the pointer to the Edge structure.
**	The radar structure is left unaffected. 
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
**	Package			-
**	Reference number	- SP1/PGM/
**	Revision number		- $Revision: 1.2 $
**	Release State		- $State: Exp $
**	Author, designer	- Don Burrows
** Modification Date		- $Date: 1999/03/31 22:35:29 $
** Modified by			- $Author: merritt $
** $Source: /nfs/trmm/src/CVS/etor/RSL_to_EDGE.c,v $
**
** MODIFICATION RECORD
**
** $Log: RSL_to_EDGE.c,v $
** Revision 1.2  1999/03/31 22:35:29  merritt
** the etor lib
**
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
#include	<string.h>
#include	<stddef.h>	/* Some popular symbols			*/
#include	<stdlib.h>	/* Some standard funct.			*/
#include	<unistd.h>	/* POSIX symbols definitions		*/
#include 	<math.h>
#include	<time.h>

/*----------------------------------------------------------------------*/
/* Application Headers                                                  */
/*----------------------------------------------------------------------*/
#include "vol.h"
#include "rsl.h"

/*----------------------------------------------------------------------*/
/* Macros                                                               */
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* External (Import) Variables                                          */
/*----------------------------------------------------------------------*/
extern int radar_verbose_flag;

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



int RSL_to_EDGE(Radar *RSL_rad, char **data)
{
	int num_sweeps;
	int ng,gs,pf;
	float nv,mr;
	int num_gates[MAX_SWEEPS];
	int num_rays[MAX_SWEEPS];
	int gate_size[MAX_SWEEPS];
	int prf[MAX_SWEEPS];
	float max_range[MAX_SWEEPS];
	float nyq_vel[MAX_SWEEPS];
	int offset[MAX_SWEEPS];
	unsigned int mem_alloc;
	int i,j,k,index;
	int nvols;
	Sweep *RSL_sweep;
	Ray *RSL_zray;
	Ray *RSL_vray;
	Ray *RSL_wray;
	struct vol_struct *EDGE_vol;
	unsigned char *EDGE_ray;
	unsigned short *sray;
	struct tm voltime;
	float wavelength,azimuth,elevation;
	float refl,vel,specw;
	float dzgsr,vrgsr,swgsr,zindex;
	int dzindex,vrindex,swindex;

	RSL_radar_verbose_on();
	if (radar_verbose_flag) printf("Month: %d  Day: %d  Year: %d City: %s\n",
		RSL_rad->h.month,RSL_rad->h.day,RSL_rad->h.year,RSL_rad->h.city);
	fprintf(stderr, "WELCOME TO RSL_to_EDGE\n");
	mem_alloc = 262108;
	num_sweeps = 0;
	nvols = 0;
	for (i=0;i<RSL_rad->h.nvolumes;i++)
	{
		if(RSL_rad->v[i] != NULL)
		{
			if (RSL_rad->v[i]->h.nsweeps > num_sweeps)
				num_sweeps = RSL_rad->v[i]->h.nsweeps;
			nvols++;
			if (radar_verbose_flag) 
				printf("volume: %d   num_sweeps = %d\n",i,RSL_rad->v[i]->h.nsweeps);
		}
	}
	if (num_sweeps > MAX_SWEEPS) num_sweeps = MAX_SWEEPS;
	if (radar_verbose_flag) printf("There are %d volumes with data and num sweeps = %d\n",
		nvols,num_sweeps);
	for (i=0;i<num_sweeps;i++)
	{
		gate_size[i] = 999999;
		prf[i] = 0;
		nyq_vel[i] = 0.0;
		max_range[i] = 0.0;
		for(j=0;j<nvols;j++) 
		{
			RSL_sweep = RSL_rad->v[j]->sweep[i];
			if (RSL_sweep == NULL) continue;
			if (RSL_sweep->ray[0] == NULL) continue;
			ng = RSL_sweep->ray[0]->h.nbins;
			gs = RSL_sweep->ray[0]->h.gate_size;
			mr = (float)ng*(float)gs;
			pf = RSL_sweep->ray[0]->h.prf;
			nv = RSL_sweep->ray[0]->h.nyq_vel;
			if (gs < gate_size[i]) gate_size[i] = gs;
			if (pf > prf[i]) prf[i] = pf;
			if (mr > max_range[i]) max_range[i] = mr;
			if (nv > nyq_vel[i]) nyq_vel[i] = nv;
		}
		num_gates[i] = max_range[i]/gate_size[i];
		max_range[i] /= 1000.0;
		if (radar_verbose_flag)
			printf("sweep: %d  gates: %d  max_range: %f\n",
				i,num_gates[i],max_range[i]);
		num_rays[i] = RSL_rad->v[0]->sweep[i]->h.nrays;
		offset[i] = mem_alloc;
		mem_alloc += num_rays[i]*(4*num_gates[i] + 8);
	}
	wavelength = 4.0*nyq_vel[0]/prf[0];
	if (radar_verbose_flag) printf("%d  bytes needed for edge volume\n",mem_alloc);
	if((EDGE_vol = (struct vol_struct *)malloc(mem_alloc)) == NULL)
	{
		fprintf(stderr,"RSL_to_EDGE: Not enough memory for volume structure\n");
		return -1;
	}
	*data = (char *)EDGE_vol;
	memset(EDGE_vol,0x00,mem_alloc);
/*  Now fill the header structures  */
	EDGE_vol->version = 0;
        voltime.tm_sec = (int)RSL_rad->h.sec;
        voltime.tm_min = RSL_rad->h.minute;
        voltime.tm_hour = RSL_rad->h.hour;
        voltime.tm_mday = RSL_rad->h.day;
        voltime.tm_mon = RSL_rad->h.month-1;
        voltime.tm_year = RSL_rad->h.year - 1900;
        voltime.tm_isdst = 0;
	printf("Volume time is %2d:%2d:%2d\n",voltime.tm_hour,voltime.tm_min,voltime.tm_sec);
/* mktime assumes that the time is localtime and converts to GMT by adding 6 hours for CST */
	EDGE_vol->date = mktime(&voltime);
/* Subtract out the Central Time offset to return to GMT */
	EDGE_vol->date -= 21600;

	EDGE_vol->scan_type = 0;
	EDGE_vol->moment_enable = 0x10;
	EDGE_vol->mem_alloc = mem_alloc;
	EDGE_vol->mem_used = mem_alloc;
	EDGE_vol->num_sweeps = num_sweeps;
	EDGE_vol->version = VOL_VERSION;

	for (i=0;i<num_sweeps;i++)
	{
		dzgsr = (float)gate_size[i]/(float)RSL_rad->v[DZ_INDEX]->sweep[i]->ray[0]->h.gate_size;
		vrgsr = (float)gate_size[i]/(float)RSL_rad->v[VR_INDEX]->sweep[i]->ray[0]->h.gate_size;
		swgsr = (float)gate_size[i]/(float)RSL_rad->v[SW_INDEX]->sweep[i]->ray[0]->h.gate_size;
		RSL_sweep = RSL_rad->v[DZ_INDEX]->sweep[i];
		if (RSL_sweep == NULL) continue;
		printf("sweep: %d  Elevation: %f\n",i,RSL_sweep->h.elev);
		EDGE_vol->sweep[i].data_offset = offset[i];
		if (radar_verbose_flag) printf("sweep: %d   num_rays: %d\n",i,RSL_sweep->h.nrays);
		EDGE_vol->sweep[i].num_rays = RSL_sweep->h.nrays;
		EDGE_vol->sweep[i].ray_size = 4*num_gates[i] + 8;
		EDGE_vol->sweep[i].date = EDGE_vol->date;
		EDGE_vol->sweep[i].max_range = num_gates[i]*gate_size[i]/1000.0;
		EDGE_vol->sweep[i].rad.prf1 = prf[i];
		EDGE_vol->sweep[i].rad.prf2 = RSL_sweep->ray[0]->h.prf;
		EDGE_vol->sweep[i].rad.gw1 = gate_size[i];
		EDGE_vol->sweep[i].rad.gw2 = RSL_sweep->ray[0]->h.gate_size;
		EDGE_vol->sweep[i].rad.pulse_width = (unsigned int)(RSL_sweep->ray[0]->h.pulse_width);
		EDGE_vol->sweep[i].rad.pulses = RSL_sweep->ray[0]->h.pulse_count;
		EDGE_vol->sweep[i].rad.gates = num_gates[i];
		EDGE_vol->sweep[i].rad.ray_size = EDGE_vol->sweep[i].ray_size;
		EDGE_vol->sweep[i].rad.fold = 0;
		EDGE_vol->sweep[i].rad.wavelength = wavelength;
		elevation = RSL_sweep->h.elev;
		if (elevation < 0.0) elevation += 360.0;
		elevation = elevation/360.0*65536.0;
		EDGE_vol->sweep[i].rad.el = (unsigned int)elevation;
		strcpy(EDGE_vol->sweep[i].rad.job_name,"UNKNOWN");
		EDGE_vol->sweep[i].rad.long_deg = RSL_rad->h.lond;
		EDGE_vol->sweep[i].rad.long_min = RSL_rad->h.lonm;
		EDGE_vol->sweep[i].rad.long_sec = RSL_rad->h.lons;
		EDGE_vol->sweep[i].rad.lat_deg = RSL_rad->h.latd;
		EDGE_vol->sweep[i].rad.lat_min = RSL_rad->h.latm;
		EDGE_vol->sweep[i].rad.lat_sec = RSL_rad->h.lats;
		EDGE_vol->sweep[i].rad.antenna_height = RSL_rad->h.height;
		strncpy(EDGE_vol->sweep[i].rad.site_name,RSL_rad->h.name,15);
/*
		strncat(EDGE_vol->sweep[i].rad.site_name,RSL_rad->h.radar_name,7);
*/
		strncpy(EDGE_vol->sweep[i].rad.radar_type,RSL_rad->h.radar_type,15);
		EDGE_vol->sweep[i].rad.bright_band = 0x00;
		EDGE_vol->sweep[i].rad.attenuation = 0x00;
		EDGE_vol->sweep[i].rad.occult = 0x00;
		for (j=0;j<RSL_sweep->h.nrays;j++)
		{
			RSL_zray = RSL_rad->v[DZ_INDEX]->sweep[i]->ray[j];
			RSL_vray = RSL_rad->v[VR_INDEX]->sweep[i]->ray[j];
			RSL_wray = RSL_rad->v[SW_INDEX]->sweep[i]->ray[j];
			if (RSL_zray == NULL) continue;
			if (RSL_vray == NULL) continue;
			if (RSL_wray == NULL) continue;
			EDGE_ray = RAY_PTR(EDGE_vol,i,j);
			sray = (unsigned short *)EDGE_ray;
			azimuth = RSL_zray->h.azimuth;
			azimuth = azimuth/360.0*65536.0;
/*
			elevation = RSL_sweep->h.elev;
			if (elevation < 0.0) elevation += 360.0;
			elevation = elevation/360.0*65536.0;
*/
			sray[0] = (unsigned short)azimuth;
			sray[1] = (unsigned short)elevation;
			sray[2] = (unsigned short)azimuth;
			sray[3] = (unsigned short)elevation;
			for (k=0,index=8;k<num_gates[i];index+=4,k++)
			{
				zindex = (float)k*dzgsr+0.5;
				dzindex = (int)zindex;
				vrindex = (int)((float)k*vrgsr+0.5);
				swindex = (int)((float)k*swgsr+0.5);
				if (dzindex <= RSL_zray->h.nbins-1)
				{
					refl = RSL_zray->h.f(RSL_zray->range[dzindex]);
					if (refl > 1000.0) refl = 0.0;
					else
					{
						refl = (refl + 32.0)*2.0;
						if (refl < 1.0) refl = 1.0;
						if (refl > 255.0) refl = 255.0;
					}
					EDGE_ray[index] = (unsigned char)refl;
					EDGE_ray[index+2] = (unsigned char)refl;
				}
				else 
				{
					EDGE_ray[index] = 0xff;
					EDGE_ray[index+2] = 0xff;
				}
				if (vrindex < RSL_vray->h.nbins)
				{
					vel = RSL_vray->h.f(RSL_vray->range[vrindex]);
					if (vel <= -nyq_vel[i]) vel = 0.0;
					else
					{
						vel = vel*128.0/nyq_vel[i]+128.0;
						if (vel < 0.0) vel += 256.0;
						if (vel > 255.0) vel -= 256.0;
						if (vel < 1.0) vel = 1.0;
						if (vel > 255.0) vel = 255.0;
					}
					EDGE_ray[index+1] = (unsigned char)vel;
				}
				else EDGE_ray[index+1] = 0xff;
				if (swindex < RSL_wray->h.nbins)
				{
					specw = RSL_wray->h.f(RSL_wray->range[swindex]);
					
					if (specw <= -nyq_vel[i]) specw = 0.0;
					else
					{
						specw = specw*128.0/nyq_vel[i];
						if (specw < 1.0) specw = 1.0;
						else if (specw > 255.0) specw = 255.0;
					}
					EDGE_ray[index+3] = (unsigned char)specw;
				}
				else EDGE_ray[index+3] = 0xff;
			}
			
		}
	}
	if (radar_verbose_flag) printf("RSL to EDGE conversion complete\n");	
	return 0;
}
/*-END OF MODULE--------------------------------------------------------*/
