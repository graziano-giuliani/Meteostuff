#ifndef _RADAR_H
#define _RADAR_H

#ifndef _RCU_H
#include "rcu.h"
#endif /* _RCU_H */
/*#include <sys/timers.h> */
#include <sys/times.h>
#include <sys/types.h>
#include <time.h>

#define DATA_RESOL 125 /* metres */
#define GATE_MASK_SIZE	16384
#define GATE_MASK_ARRAY_SIZE (GATE_MASK_SIZE/16)
#define NUM_PW	4

#define speed_of_light (1.5e8)
#define U_RANGE (rad->prf1!=0.0)?(unsigned int)(speed_of_light/(float)(rad->prf1)):0

#define MAX_SITE_NAME_LENGTH 16
#define MAX_RADAR_TYPE_LENGTH 16

#define RVP5 0
#define RVP6 1
#define ESP7 2

#define RC_HAS_NO_PRIM_BITE 1
#define RC_HAS_NO_SEC_BITE 2
#define RC_HAS_NO_BITE 3
#define RC_HAS_NO_TSG  4
#define RC_VEL_INVERT  8

#define RADAR_HAS_PRIM_BITE ((rad->radar_config & RC_HAS_NO_PRIM_BITE) == 0)
#define RADAR_HAS_SEC_BITE ((rad->radar_config & RC_HAS_NO_SEC_BITE) == 0)
#define RADAR_HAS_BITE ((rad->radar_config & RC_HAS_NO_BITE) != 3)
#define RADAR_HAS_TSG ((rad->radar_config & RC_HAS_NO_TSG) == 0)
#define RADAR_VEL_INVERT ((rad->radar_config & RC_VEL_INVERT) != 0)

struct trig_struct
{
	float start;
	float width;
	int polar;
};

/*
**************************************************************************
This structure is included in the volume structure. If you change this 
structure you must change the version numbe in the vol.h file and if
necessary write a translation routine
NOTE:  THIS BABY IS 7516 BYTES IN SIZE! Wed Mar 20 11:20:32 CST 1996
**************************************************************************
*/
struct radar_struct
{

	unsigned int software_sim;	/* 0 off 1 on */
	char sim_file[255];
	unsigned char sim_vel;

	unsigned int prf1,prf2;		/* Hz */

	unsigned int gw1,gw2;		/* metres */
	float gw_partition;		/* Km */
	unsigned char range_avg;
	unsigned int pulse_width;
	unsigned int wave;
	unsigned int pulses;
	/* unsigned char range_norm; */
	unsigned char unused;
	unsigned short clutter_filter;
	float spec_range;
	unsigned short unused2[509];
	unsigned short autopow;

	unsigned short gate_mask[GATE_MASK_ARRAY_SIZE];
	unsigned int gates;
	unsigned int ray_size;
	unsigned short samples;
	unsigned short flags;
	unsigned short sqi;
	unsigned short fold;
	unsigned short log_slope[NUM_PW];
	unsigned short noise_thresh[NUM_PW];
	short clutter_thresh[NUM_PW];
	unsigned short wsp_threshold[NUM_PW];
	short calibration_reflectivity[NUM_PW];
	unsigned short agc_pulses_integrated;
	unsigned short filter_delay;
	unsigned short uncorr_reflec_flags;
	unsigned short corr_reflec_flags;
	unsigned short vel_thresh_flags;
	unsigned short wid_thresh_flags;
	unsigned short tag_invert_ls;
	unsigned short tag_invert_ms;
	unsigned short gas_atten_corr;
	unsigned short zdr_flags;
	short zdr_calib_offset;
	float wavelength;

	int noise_int;
	unsigned int last_noise;



	unsigned char moment_enable;

	unsigned char scan_type;
	int dir;
	unsigned int antenna_speed;		
	unsigned int az,el;			/* ACP's */
	unsigned int rhi_top,rhi_bot;		/* ACP's (sort of) */
	unsigned int sector_start,sector_end;	/* ACP's */
	unsigned int target_az;
	unsigned int target_el;

	int long_deg;
	unsigned int long_min;
	unsigned int long_sec;
	
	int lat_deg;
	unsigned int lat_min;
	unsigned int lat_sec;

	int antenna_height;
	int agcpol;


	int max_prf[4];
	unsigned short calib_data[81];	
	float lin_pow[81];	

	unsigned short gparm[64];
	unsigned short rback[240];

	char site_name[MAX_SITE_NAME_LENGTH];
	char radar_type[MAX_RADAR_TYPE_LENGTH];
	char arch_host[64];

	char job_name[64];

	char proc_description[200];


	unsigned char rcu_status[9];
	unsigned char rcu_bite[2][10];

	struct rcu_control rcu;

	unsigned short az_speed;
	unsigned short el_speed;
	unsigned short last_az;
	unsigned short last_el;
	int min_prf[4];
	int32_t bite_time;
	unsigned int tsg_err;
	unsigned int stat_filt_thold;   /* stat. clut. filt. thold *1000*/
	struct spare_struct
	{
		unsigned short angles[2];
	} spare14;
	float atten_factor;		/* (dB/km) / (mm/hr)		*/
	unsigned char clutter_map; 	/* 0=no map, 1=clutter marked,	*/
					/* 2=clutter interpolated, 	*/
					/* 3=clutter subtracted		*/
	unsigned char bright_band;	/* 0=off, 255=none, 1-254=km*10 */
	unsigned char attenuation;	/* non-zero, has been corrected */
	unsigned char occult;		/* zero = OFF			*/

	/* struct bite_struct bite; */
	char site_code[4];
	unsigned int cmap_flag;
	unsigned short idleang;
	float za_a, za_b;
	float a,b;

	char vdest[1024];
	char pdest[1012];

	float beam_width;
	unsigned int radar_config;

	short sigproc;
	unsigned short interface_type;
	short tags;

	struct trig_struct trigger[4][6];

	/*
	** other good stuff in here
	*/
};

struct stc_struct 
{
	unsigned short stc_wave[8192];
};

struct caldat_struct
{
	unsigned short cal[4][2048];
};

#define RAD_SHM_ID	0xb0ba

#define PPI_SCAN 	0
#define RHI_SCAN 	1
#define SECTOR_SCAN 	2
#define POINT_SCAN 	3

#define ANT_UP		1
#define ANT_DOWN 	-1
#define ANT_CLOCK	1
#define ANT_COUNTER_CLOCK	-1


#define SCAN_TYPE(s) ((s)==PPI_SCAN?"PPI":(s)==RHI_SCAN?"RHI":\
	(s)==SECTOR_SCAN?"SECTOR":"UNKNOWN")

#define SPEEDFACTOR 1.0
#define CVSPEED(s) ((int)((s+SPEEDFACTOR/2.0)/SPEEDFACTOR))
#define DECODE_SPEED(I) ((float)((I))*SPEEDFACTOR)
#define CVELEV(e) ((int)((e)*65536.0/360.0))
#define CVAZ(a) ((int)((a)*65536.0/360.0))

#define GBIT(i) (1<<(((i)&0xf)))
#define GEL(i) ((i)/16)

#ifdef NEVER
#define MASK_ANG(A) ((unsigned short)((A) & 0xffff))
#define ELDEGS(a) (printf("%d %d\n",(a),MASK_ANG(a))),(((MASK_ANG(a))>0x8000?(int)(((a))-65536):((a)))/65536.0*360.0)
#endif
#define MASK_ANG(A) ((unsigned short)((A) & 0xffff))
#define ELDEGS(a) (((signed short)(a))/65536.0*360.0)
#define DEGS(a) (((a))/65536.0*360.0)
#define IDEGS(a) ((int)(((a))/65536.0*360.0))
#define ACPS(a) ((int)(((a))*65536.0/360.0))

#define SUB_AZ(a,b) (((a)-(b))&0xffff)
#define DIFF_AZ(a,b) (SUB_AZ((a),(b))<SUB_AZ((b),(a))?SUB_AZ((a),(b)):SUB_AZ((b),(a)))
#define DIFF_EL(a,b)	DIFF_AZ((a),(b))

#define RCUAZ ((rad->rcu_status[2]<<9)+(rad->rcu_status[1]<<2))
#define RCUEL ((rad->rcu_status[4]<<9)+(rad->rcu_status[3]<<2))

#endif /* _RADAR_H */
