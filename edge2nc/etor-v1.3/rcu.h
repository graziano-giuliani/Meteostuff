#ifndef _RCU_H
#define _RCU_H

#define CONTROL_SEM_ID 0xed9e

struct bite_struct
{
	unsigned char tx_rec_p5v;
	unsigned char tx_rec_m15v;
	unsigned char tx_rec_p15v;
	unsigned char tx_rec_p28v;
	unsigned char lo_p15v;
	unsigned char servo_p15v;
	unsigned char servo_m15v;
	unsigned char servo_p24v;
	unsigned char stalo_lock;
	unsigned char tx_trigger;
	unsigned char tx_hv;
	unsigned char p28v_interlock;
	unsigned char p28v_standby;
	unsigned char fwd_power;
	unsigned char rvs_power;
	unsigned char p28v_radiate;
	unsigned char pw_select;
	unsigned char air_flow;
	unsigned char waveguide;
	unsigned char mag_current;
	unsigned char cabinet_temp;
	unsigned int aux;
};

struct rcu_status
{
	unsigned short angles[2];
	unsigned int low_air;
	unsigned int low_wg_pressure;
	unsigned int servo_power;
	unsigned int radar_power;
	unsigned int interlock	;
	unsigned int standby	;
	unsigned int radiate	;
	unsigned int RCP_down	;
	unsigned int pw		;
	unsigned int local_mode	;
	unsigned int mag_current;
};

struct rcu_control
{
	unsigned short az;
	unsigned short el;
	unsigned int pw		;
	unsigned int cw		;
	unsigned int el_scan	;
	unsigned int az_scan	;
	unsigned int rcp_reset	;
	unsigned int sig_gen_on	;
	unsigned int radiate_off;
	unsigned int radiate_on	;
	unsigned int servo_on	;
	unsigned int radar_power;
	unsigned char sig_gen_level;
	unsigned char tsg[6];
	unsigned char aux;
	char speed;
	unsigned char speed1;
};
#endif /* _RCU_H */
