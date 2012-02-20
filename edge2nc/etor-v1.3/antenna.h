#ifndef ANTENNA_H
#define ANTENNA_H

#undef RAD
#define RAD(A) ((float)(A) / 57.2957795131)

#define Re (8498670.0)		/* 4/3 earth radius			*/
#define Resq (7.222739177E13)	/* 4/3 earth radius squared		*/
/*
** Convert binary angles to nearest integer angle.
*/
#define BIN2IANG(B) (((B) * 45) >> 13)		/* Binary angle to angle */
#define IANG2BIN(B) (((B) << 13) / 45)
/*
** Binary angle to integer angle times 10 & 100
*/
#define BIN2IANG10(B) (((B) * 225) >> 12)
#define BIN2IANG100(B) (((B) * 1125) >> 11)
/*
** Binary angle to signed integer angle and 10 times angle
*/
#define oBINEL2IANG(B) (((B) & 0x8000) ?  (BIN2IANG(B)-360) : BIN2IANG(B))
#define oBINEL2IANG10(B) (((B) & 0x8000) ?  (BIN2IANG10(B)-360) : BIN2IANG10(B))
#define BINEL2IANG(B) (((B) & 0x8000) ?  (-BIN2IANG((unsigned short)-(B))) : BIN2IANG(B))
#define BINEL2IANG10(B) (((B) & 0x8000) ?  (-BIN2IANG10((unsigned short)-(B))) : BIN2IANG10(B))
#define BINEL2IANG100(B) (((B) & 0x8000) ?  (-BIN2IANG100((unsigned short)-(B))) : BIN2IANG100(B))

#endif
