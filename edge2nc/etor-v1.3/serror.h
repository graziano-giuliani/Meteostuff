/*
** serror.h -- define error system data structures
*/
#ifndef SERROR_H
#define SERROR_H

#include "err.h"

#ifndef __EXE__
#define __EXE__ "*"
#endif

typedef void (*PFV)();
/*
** All EDGE errors have the high bit set, warnings also have bit 14 set.
*/
#define EDGE_ERR_BASE 	0x8000
#define EDGE_ERR_MASK 	0x1fff
#define WARNING_BIT 	0x4000
#define INFO_BIT 	0x2000
#define MAX_ERRS	0x1fff
/*
** FATAL -- returns true if argument is a EDGE fatal error number
*/
#define FATAL(E)	(((E) & WARNING_BIT) == 0)
/*
** WARNING -- returns true if argument is a EDGE warning message
*/
#define WARNING(E)	(((E) & WARNING_BIT) != 0)
/*
** EDGE_ERROR -- Returns true if argument is a valid EDGE error number
*/
#define EDGE_ERROR(E)	((E) & EDGE_ERR_BASE)
/*
** ERR_NAME, NAME, name -- assign argument to EDGE informational error 
*/
#define ERR_NAME(S)	name(S)
#define NAME(S)		name(S)
/*
** seterr -- set parameters for _seterr
*/
#define error		err_id(__LINE__,__FILE__),_error
#define seterr		err_id(__LINE__,__FILE__),_seterr

#define OK		0		/* Good return value		*/
#define ERROR		-1		/* Bad return value		*/
#define ERR		ERROR		/* Bad return value		*/

/*
** errno definitions above 0x8000, used for EDGE
*/
#define NO_ERR		0x8000
#define ARG_ERR		0x8001
#define TST_1_ERR	0x8002
#define TST_2_ERR	0x8003
#define WARNING_WAR	0xc004
#define NO_VGA_ERR	0x8005
#define NO_POOL		0x8006
#define SP_READ		0x8007
#define JOB_FULL	0x8008
#define JOB_NULL	0x8009
#define JOB_INVALID	0x800a
#define SCD_ON		0x800b
#define SCD_OFF		0x800c
#define FOPEN_ERR	0x800d
#define NO_MEM		0x800e
#define NOT_IMP		0x800f
#define VOL_READ	0x8010
#define RTD_BROAD	0x8011
#define DCT_ERR		0x8012
#define UNK_TYPE	0x8013
#define ACQ_ERR		0x8014
#define ARCH_ERR	0x8015
#define ANT_ERR		0x8016
#define TSG_FLT		0x8017
#define FMT_ERR		0x8018
#define DIAG_ERR	0x8019
#define NOIS_ERR	0x801A
#define SIGP_ERR	0x801B
#define TIME_OUT	0x801C
#define SIGP_DIFF	0x801D
#define OPT_ERR		0x801E
#define PRF_ERR		0x801F
#define WAVE_ERR	0x8020
#define PW_ERR		0x8021
#define SAMP_ERR	0x8022
#define REX_ERR		0x8023
#define TIM_ERR		0x8024
#define PRD_ERR		0x8025
#define JOB_RUN		0x8026
#define VOL_INVALID	0x8027
#define COMP_SHM	0x8028
#define SURV_OFF	0x8029
#define SURV_ON		0x802a
#define RAING_ERR	0x802b
#define CSITE_ERR	0x802c
#define EDGEPC_ERR	0x802d
#define SURV_SIGNAL	0x802e
#define STAT_ATT	0x802f
#define NOT_LIC		0x8030
#define SHM_ERR		0x8031
#define NOT_DIR		0x8032
#define MANY_RAYS	0x8033
#define ANGLE_ERR	0x8034

#define TX_RX_P15F	0x2001
#define TX_RX_P15R	0x2002
#define TX_RX_M15F	0x2003
#define TX_RX_M15R	0x2004
#define TX_RX_P28F	0x2005
#define TX_RX_P28R	0x2006
#define FAULT_CLEARED	0x2007
#define STALO_FAULT	0x2008
#define STALO_REC	0x2009
#define HV_FAULT	0x200a
#define HV_REC		0x200b
#define INLCK_FAULT	0x200c
#define INLCK_REC	0x200d
#define FWD_POW_F	0x200e
#define FWD_POW_R	0x200f
#define RVS_POW_F	0x2010
#define RVS_POW_R	0x2011
#define RAD_FAULT	0x2012
#define RAD_REC		0x2013
#define AIR_FLOW_F	0x2014
#define AIR_FLOW_R	0x2015
#define WG_PRESS_F	0x2016
#define WG_PRESS_R	0x2017
#define MAG_CURR_F	0x2018
#define MAG_CURR_R	0x2019
#define CAB_TEMP_F	0x201a
#define CAB_TEMP_R	0x201b
#define TX_RX_P5F	0x201c
#define TX_RX_P5R	0x201d
#define LO_P15F		0x201e
#define LO_P15R		0x201f
#define SERVO_P15F	0x2020
#define SERVO_P15R	0x2021
#define SERVO_M15F	0x2022
#define SERVO_M15R	0x2023
#define SERVO_P24F	0x2024
#define SERVO_P24R	0x2025

#ifndef errno
#include <errno.h>
#endif

#define MAX_ERROR_STACK 20

#define serror()	_serror(__EXE__, __FILE__, __LINE__,NULL)
#define geterror(T)	_serror(__EXE__, __FILE__, __LINE__,(T))


#endif
