/*----------------------------------------------------------------------**
**
** serror.c - Error Message Processing
**
**----------------------------------------------------------------------**
**
** DESCRIPTION
**
** This function is used to report errors detected by the calling program
** It uses the parameters and the value of the variables errno and edge_errno
** to produce a meaningful error message.
**
**
** USAGE:
**
** void _serror(
** 	char *prog,	-* Program name				*-
** 	char *src,	-* Source file name			*-
** 	int line	-* Source file line number		*-
** )
**
** PROCESSING:
**
** The first thing this function does is call the user defined routine pointed
** to by the pre_message variable (if the pointer is no zero). If errno is
** none zero the program name, time source file name, line number and error 
** message are printed to standard error. If edge_errno is none zero a message 
** starting "EDGE:" and containig the program name, time, "Warning" or "Fatal"
** Source file, line number, error number and message is printed. Any error
** messages placed in the error_stack by seterr are printed next. And finally
** the used defined routine pointed to by the pointer post_message is called.
** Post_message is passed the error level and the message, NUll   Terminated.
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
**	Package			- ERR
**	Reference number	- SP1/PGM/
**	Revision number		- $Revision: 1.2 $
**	Release State		- $State: Exp $
**	Author, designer	- KenB
** Modification Date		- $Date: 2001/05/04 19:19:35 $
** Modified by			- $Author: kelley $
** $Source: /nfs/trmm/src/CVS/etor/serror.c,v $
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
#define _OSF_SOURCE	1

/*----------------------------------------------------------------------*/
/* System Headers                                                       */
/*----------------------------------------------------------------------*/
#include 	<stdarg.h>	/* Var. argument funct. support		*/
#include	<stdlib.h>
#include	<stdio.h>	/* stdio library			*/
#include	<string.h>	/* String functions			*/
#include	<time.h>	/* Time-of-day functs			*/
#include 	<errno.h>	/* All error codes			*/
#if defined(__linux)
#include	<sys/time.h>
#else
#include	<sys/timers.h>
#endif

/*----------------------------------------------------------------------*/
/* Application Headers                                                  */
/*----------------------------------------------------------------------*/
#include "serror.h"

/*----------------------------------------------------------------------*/
/* External (Import) Variables                                          */
/*----------------------------------------------------------------------*/
extern void (*pre_message)(void);
extern void (*post_message)(int, char *);
extern char *err_src;
extern int err_line;

/*----------------------------------------------------------------------*/
/* External Functions                                                   */
/*----------------------------------------------------------------------*/
extern char *err_id( int, char * );

/*----------------------------------------------------------------------*/
/* Structures and Unions                                       	        */
/*----------------------------------------------------------------------*/

/*
** A breakdown of the string returned by ctime()
*/
struct ctime_structure
{
	char weekday[3];
	char space;
	char mon[3];
	char space1;
	char day[2];
	char space2;
	char time[8];
	char space3[3];
	char year[2];
	char stuff[2];
};

/*----------------------------------------------------------------------*/
/* Global (Export) Variables                                            */
/*----------------------------------------------------------------------*/

char *error_stack[MAX_ERROR_STACK];
char *error_info = NULL;
int errc = 0;
char error_message[512]; /* Temp storage for message		*/
int serror_action=0;

char *edge_errlist[] =
{
	"NO_ERR, Not an error, errno=0",
	"ARG_ERR, Argument Error",
	"TST_1_ERR, Test fatal 1 error",
	"TST_2_ERR, Test fatal 2 error",
	"WARNING_WAR, Warning message test",
	"NO_VGA_ERR, Can't access VGA",
	"NO_POOL, Pool exausted",
	"SP_READ, SP read truncated",
	"JOB_FULL, Job Table Full",
	"JOB_NULL, Job Table Not Initialised",
	"JOB_INVALID, Invalid Job Id",
	"SCD_ON, Scheduler already running",
	"SCDOFF, Scheduler not running",
	"FOPEN_ERR, fopen() failed",
	"NO_MEM, No memory fo malloc()",
	"NOT_IMP, Feature not yet implemented",
	"VOL_READ, Error reading volume file",
	"RTD_BROAD, Real Time Data broadcast error",
	"DCT_ERR, DCT driver or device error",
	"UNK_TYPE, Unknown file type",
	"ACQ_ERR, Error acquiring raw data",
	"ARCH_ERR, Archival error",
	"ANT_ERR, Antenna error",
	"TSG_FLT, Test Signal Generator Fault",
	"FMT_ERR, Message format error",
	"DIAG_ERR, Radar Unavailable",
	"NOIS_ERR, Unable to get good noise sample",
	"SIGP_ERR, Command not available on current signal processor",
	"TIME_OUT, Time-out on device I/O",
	"SIGP_DIFF, Difference between SOPRM and GPARM",
	"OPT_ERR, Unrecognized or illegal option",
	"PRF_ERR, Illegal PRF value",
	"WAVE_ERR, Illegal WAVE value",
	"PW_ERR, Illegal pulse width index",
	"SAMP_ERR, Illegal Number Samples",
	"REX_ERR, Can not resolve 'rex' IP address",
	"TIM_ERR, Time conversion error",
	"PRD_ERR, Product Script Error",
	"JOB_RUN, Job Already running",
	"VOL_INVALID, Invalid volume number",
	"COMP_SHM, Unable to attach to composite shared memory",
	"SURV_OFF, Surveillance Scan is Off",
	"SURV_ON, Surveillance Scan is On",
	"RAING_ERR, Raing Gauge Error",
	"CSITE_ERR, Not WARN central site",
	"EDGEPC_ERR, Error in EDGE PC configuration",
	"SURV_SIGNAL, Can not signal SURV",
	"STAT_ATT, Can't attach to edge_stat shared memory",
	"NOT_LIC, Feature not licensed",
	"SHM_ERR, Error connecting to shared memory",
	"NOT_DIR, Not a directory",
	"MANY_RAYS, Too many rays in sweep",
	"ANGLE_ERR, Error in antenna angles"
};
char *edge_info[] =
{
	"NO_ERR, Not an error, errno=0",
	"TX_RX_P15F, TX/RX +15Vdc Failure",
	"TX_RX_P15R, TX/RX +15Vdc Recovered",
	"TX_RX_M15F, TX/RX -15Vdc Failure",
	"TX_RX_M15R, TX/RX -15Vdc Recovered",
	"TX_RX_P28F, TX/RX +28Vdc Failure",
	"TX_RX_P28R, TX/RX +28Vdc Recovered",
	"FAULT_CLEARED, BITE shows previous failure cleared",
	"STALO_FAULT, Stalo Unlocked",
	"STALO_REC, Stalo Locked",
	"HV_FAULT, High Voltage Failure",
	"HV_REC, High Voltage Recovered",
	"INLCK_FAULT, Interlock Failure",
	"INLCK_REC, Interlock Recovered",
	"FWD_POW_F, Forward Power Failure",
	"FWD_POW_R, Forward Power Recovered",
	"RVS_POW_F, Reverse Power Failure",
	"RVS_POW_R, Reverse Power Recovered",
	"RAD_FAULT, Radiate Failure",
	"RAD_REC, Radiate Recovered",
	"AIR_FLOW_F, Air Flow Failure",
	"AIR_FLOW_R, Air Flow Recovered",
	"WG_PRESS_F, Waveguide Pressure Failure",
	"WG_PRESS_R, Waveguide Pressure Recovered",
	"MAG_CURR_F, Magnetron Current Failure",
	"MAG_CURR_R, Magnetron Current Recovered",
	"CAB_TEMP_F, Cabinet Temperature Failure",
	"CAB_TEMP_R, Cabinet Temperature Recovered",
	"TX_RX_P5F, TX/RX +5Vdc Failure",
	"TX_RX_P5R, TX/RX +5Vdc Recovered",
	"LO_P15F, Local Oscillator +15Vdc Failure",
	"LO_P15R, Local Oscillator +15Vdc Recovered",
	"SERVO_P15F, Servo +15Vdc Failure",
	"SERVO_P15R, Servo +15Vdc Recovered",
	"SERVO_M15F, Servo -15Vdc Failure",
	"SERVO_M15R, Servo -15Vdc Recovered",
	"SERVO_P24F, Servo +24Vdc Failure",
	"SERVO_P24R, Servo +24Vdc Recovered"
};

int edge_errno = 0;

void clrerr( void )
{
	int i;
	/*
	** Process the error stack, printing each entry, then
	** freeing the space.
	*/
	for(i=errc-1; i>=0; i--)
	{
		free(error_stack[i]);
	}
	err_line = -1;
	err_src = (char *)0;
	edge_errno = 0;
	errno = 0;
	errc = 0;
	error_info = (char *)0;
	error_info = NULL;	/* Printed as extra info with errors	*/
}
char *strerr( int edge_errno )
{
	char msg[128];
	int err;

	err = edge_errno & EDGE_ERR_MASK;
	if((err > MAX_ERRS) || (msg == NULL))
		return "No Message";
	else
		return edge_errlist[err];
}
/*----------------------------------------------------------------------*/
/* Main Function                                                        */
/*----------------------------------------------------------------------*/

void _serror(
	char *prog,	/* Program name				*/
	char *src,	/* Source file name			*/
	int line,	/* Source file line number		*/
	char *txt
)
{
	time_t time_is;	/* A place for time			*/
	struct ctime_structure *ct;
	int err;
	char *msg;
	int i;
	int id;
	int err_lev = 9;
	int exit_num = 0;
	char buff[256];
#if defined(__linux)
	struct timeval tp;
#else
	struct timespec tp;
#endif
	int local_errno;

#ifdef DUMP_ARGS
	fprintf(stderr, "_serror(prog:'%s',src:'%s',line:%d)\n", 
		prog, src, line);
#endif /* DUMP_ARGS */

	/*
	** Use err_id to trim unwanted dots and slashes from
	** the source file name, src.
	*/
	src = err_id(line, src);

	/*
	** If we are printing an error, call pre_message routine
	*/
	if((edge_errno == 0) && (errno == 0))
	{
		errno = 0;
		edge_errno = NO_ERR;
	}
	if(pre_message != (PFV)0)
		pre_message();
	/* 
	** Get system time into structure
	*/
	time_is = time((time_t *)0);
	/*
	** make a copy of errno as ctime may corrupt it (maybe when reading 
	** time zone info)
	*/
	local_errno = errno;
	ct = (struct ctime_structure *)ctime(&time_is);
	errno = local_errno;

#define TIMEOFDAY 1
#if defined(__linux)
	gettimeofday(&tp, NULL);
	id = tp.tv_sec * 10 + (tp.tv_usec/1000);
#else
	getclock(TIMEOFDAY, &tp);
	id = (0x0000ffff & (tp.tv_sec * 10)) + (tp.tv_nsec/100000000);
#endif
	/*
	** If the system errno is set, process the system error.
	*/
	if(errno)
	{
		printf("error_info=%s\n",error_info);
		sprintf(buff, 
		  "!%04x:SYS:%s(%d):%02.2s% 3.3s% 2.2s:% 8.8s:%s(%d):%s:%s\n",
			id,
			prog,getpid(),
			ct->day, ct->mon, ct->year, ct->time,
			src, line, strerror(errno),
			error_info ? error_info : "");
		if (txt != NULL && error_info != NULL)
			strncpy(txt,error_info,128);
		fprintf(stderr,buff);
		err_lev = 8;
		exit_num=errno;
	}
	/*
	** If the edge errno is set, process that one too.
	*/
	if(edge_errno & INFO_BIT)
	{

		err = edge_errno & EDGE_ERR_MASK;
		
		msg = edge_info[err & 0xff];
		if((err > MAX_ERRS) || (msg == NULL))
			msg = "No Message";
		sprintf(buff, 
	"!%04x:INFO:%s(%d):%02.2s% 3.3s% 2.2s:% 8.8s:%s:%s(%d):I%02x:%s:%lx:%s\n", 
			id,
			prog,getpid(),
			ct->day, ct->mon, ct->year, ct->time,
			WARNING(edge_errno) ? "Warning" : "Fatal",
			src, line, err, msg,
			error_info,
			((error_info != NULL) ? error_info : ""));
		if (txt != NULL && error_info != NULL)
			strncpy(txt,error_info,128);
		fprintf(stderr,buff);
		/* printf("err>>8:%x (err>>8)+1:%x\n",err>>8,(err>>8)+1); */
		err_lev=(err>>8)+1;
	}
	else if (edge_errno)
	{

		sprintf(buff, 
	"!%04x:EDGE:%s(%d):%02.2s% 3.3s% 2.2s:% 8.8s:%s:%s(%d):E%02x:%s:%lx:%s\n", 
			id,
			prog,
			getpid(),
			ct->day, ct->mon, ct->year, ct->time,
			WARNING(edge_errno) ? "Warning" : "Fatal",
			src, line, err, strerr(edge_errno),
			error_info,
			((error_info != NULL) ? error_info : ""));
		if (txt != NULL && error_info != NULL)
			strncpy(txt,error_info,128);
		fprintf(stderr,buff);
		err_lev = 8;
		exit_num=errno;
	}
	/*
	** Process the error stack, printing each entry, then
	** freeing the space.
	*/
	if (txt != NULL )
	{
		strncpy(txt,error_stack[errc-1]+
			strcspn(error_stack[errc-1],":")+1,128);
	}
	for(i=errc-1; i>=0; i--)
	{
		sprintf(buff, "%04x:\t%s\n", id, error_stack[i]);
		fprintf(stderr,buff);
		free(error_stack[i]);
	}
	fprintf(stderr, "\n");
	if(post_message != (PFV)0)
		post_message(err_lev, buff);

	if (serror_action != 0 && exit_num != 0)
		exit(exit_num);
	/*
	** Clean up
	*/
	edge_errno = 0;
	errno = 0;
	errc = 0;
	error_info = (char *)0;
	clrerr();
}
/*-END OF MODULE--------------------------------------------------------*/
