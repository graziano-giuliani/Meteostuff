#ifndef _VOL_H
#define _VOL_H

#include <time.h>
#include <sys/types.h>
#include "radar.h"

#define VOL_VERSION 103
#define MAX_SWEEPS 30


#define NEW_VOL_SIZE	0x1000000	/* start with 8MB */
#define VOL_INC		0x400000	/* Increase in 4MB steps */

struct sweep_struct
{
	int data_offset;
	int num_rays;
	int ray_size;
	int32_t date;
	float max_range;
	struct radar_struct rad;
};

struct vol_struct
{
	int version;
	int32_t date;
	int scan_type;
	char moment_enable;
	unsigned int mem_alloc;
	unsigned int mem_used;
	int num_sweeps;
	struct sweep_struct sweep[MAX_SWEEPS];
};

#define FREE_MEM(v) (v->mem_alloc-v->mem_used)
#define NEXT_FREE(v) ((unsigned char *)(v)+((v)->mem_used))

#define SWEEP(V,T) ((V)->sweep[(T)])
#define RAYS(V,T) ((V)->sweep[(T)].num_rays)
#define RAY_SIZE(V,T) SWEEP(V,T).ray_size
#define RAY_PTR(V,T,R) ((unsigned char *)(V)+SWEEP(V,T).data_offset+ \
	(R)*RAY_SIZE(V,T))
#define AZ(V,T,R) (*(unsigned short *)(RAY_PTR(V,T,R)))
#define EL(V,T,R) (*(unsigned short *)(RAY_PTR(V,T,R)+2))
#define GATES(V,T) (SWEEP(V,T).rad.gates)
#define GATE_SIZE(V,T) (SWEEP(V,T).rad.gw1)
#define END_AZ(P)	(*(unsigned short *)((P)+4))
#define END_EL(P)	(*(unsigned short *)((P)+6))
#define START_AZ(P)	(*(unsigned short *)((P)+0))
#define START_EL(P)	(*(unsigned short *)((P)+2))


#define START_AZ_DEGS(P)        DEGS(START_AZ(P))
#define END_AZ_DEGS(P)          DEGS(END_AZ(P))
#define START_EL_DEGS(P)        DEGS(START_EL(P))
#define END_EL_DEGS(P)          DEGS(END_EL(P))

#define SEP(A,B)        ((A)>(B)?(A)-(B):(B)-(A))
#define SEPARATION(A,B) (SEP(A,B)>180.0?360.0-SEP(A,B):SEP(A,B))
#define RAY_SEPARATION(V,S,R1,R2) SEPARATION(START_AZ_DEGS(RAY_PTR(V,S,R1)),START_AZ_DEGS(RAY_PTR(V,S,R2)))


#define MAX_COMPRESS 2
#define NO_COMPRESS 0
#define LZW 1
#define RUN_LEN_ENC 2

#define ARC_MOMENT 16
#define W_MOMENT 1
#define V_MOMENT 2
#define U_MOMENT 4
#define Z_MOMENT 8
#define SHR_MOMENT 16
#define MM_MOMENT 32
#define R_MOMENT 64  /* Rainfall */
#define D_MOMENT 128 /* Diff. Reflec */
#define ZDR_MOMENT D_MOMENT /* Diff. Reflec */
#define ZDR_ENABLE(V) (((V)->moment_enable&ZDR_MOMENT) != 0 )
#define BYTES_BIN(V) ((ZDR_ENABLE(V)) +4 )
#define C_MOMENT 256 /* Combined Shear */
#define T_MOMENT 512 /* Total Rainfall */
#define DV_MOMENT 1024 /* Dealiased Velocity */
#define H_MOMENT 2048 /* Hydrometeor Value */
#define I_MOMENT 4096
#define Q_MOMENT 8192
#define L_MOMENT 16384
#define S_MOMENT (32768+12)

#define STATUS_MOMENT 0x0020 /* used for transmitting status packets */
#define MSG_MOMENT 0x0043 /* used for transmitting message packets */
#define GPARM_MOMENT 0x0085 /* used for transmitting message packets */
#define GPARM_SIZE	128
#define DCT_MOMENT 0x0086 /* used for transmitting message packets */
#define DCT_SIZE	200
#define IOTESTER_MOMENT	0x0087

#define VOL_FILE 0
#define PROD_FILE 1
#define HEADER_ONLY 0x100
#define PARTIAL_VOL 0x200

	 /* Function prototypes. */
int load_data(
	char **data, 			/* Point to read data pinter	*/
	char *filename, 		/* File to read			*/
	unsigned int filetype );		/* Sort of file, see defs	*/


#endif /* _VOL_H */

