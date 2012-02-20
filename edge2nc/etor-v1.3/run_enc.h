
#ifndef _RUN_ENC_H
#define _RUN_ENC_H


/*
 * values for first word of transmission
 */
#define TX_PRODUCT 0x88888888
#define TX_ENDOFTX 0xeeeeeeee
#define TX_ECHOTST 0x55555555
#define TX_VOLUME  0x99999999

#define TX_NULL 0x00    /* terminates ray or header, filler character   */
#define TX_SOT  0x1a    /* start of product transmission                */
#define TX_EOT  0x1b    /* end of product transmission                  */
#define TX_SOR  0x1c    /* start of ray                                 */
#define TX_EOH  0x1d    /* end of or header                             */
#define TX_SOH  0x1e    /* start of header                              */
#define TX_2BYTE        0x1f    /* 2-byte format                        */
#define TX_MES	0x3a	/* Start of Message				*/
#define TX_STS	0x3b	/* Start of status packet			*/
#define TX_SOC	0x3c	/* Start of Command				*/
#define TX_RES	0x3d	/* Start of result				*/
#define TX_RESV1        0x3a;   /* 3e is usable 		        */
#define TX_RESV2        0x5a;   /* 5a, 5b, 5c, 5d, 5e are usable        */
#define TX_RESV3        0x7a;   /* 7a, 7b, 7c, 7d, 7e are usable        */
#define TX_RESV4        0x9a;   /* 9a, 9b, 9c, 9d, 9e are usable        */
#define TX_RESV5        0xba;   /* ba, bb, bc, bd, be are usable        */
#define TX_RESV6        0xda;   /* da, db, dc, dd, de are usable        */

#endif 

/******************************* MODULE_END ********************************/
