#ifndef _HEADER_POSDIFFERRENTIAL_20150331_
#define _HEADER_POSDIFFERRENTIAL_20150331_
#include "../Common/Common.h"
#include "../Common/GNSSERROR.h"
#include "../Common/BaseFunction.h"
#include "../Common/SVEphemeris.h"
#define NFREQ       3                   /* number of carrier frequencies */


#ifdef GNSSPROCESS_EXPORTS
#define GNSSPROCESS_API _declspec(dllexport)
#else
#define GNSSPROCESS_API _declspec(dllimport)
#endif


typedef struct {        /* observation data record */
	gtime_t time;       /* receiver sampling time (GPST) */
	unsigned char sat,rcv; /* satellite/receiver number */
	double L[NFREQ]; /* observation data carrier-phase (cycle) */
	double P[NFREQ]; /* observation data pseudorange (m) */
} obsd_t;


typedef struct {        /* GPS/QZS/GAL broadcast ephemeris type */
	int sat;            /* satellite number */
	int iode,iodc;      /* IODE,IODC */
	int sva;            /* SV accuracy (URA index) */
	int svh;            /* SV health (0:ok) */
	int week;           /* GPS/QZS: gps week, GAL: galileo week */
	int code;           /* GPS/QZS: code on L2, GAL: data sources */
	int flag;           /* GPS/QZS: L2 P data flag */
	gtime_t toe,toc,ttr; /* Toe,Toc,T_trans */
	double A,e,i0,OMG0,omg,M0,deln,OMGd,idot;
	double crc,crs,cuc,cus,cic,cis;
	double toes;        /* Toe (s) in week */
	double fit;         /* fit interval (h) */
	double f0,f1,f2;    /* SV clock parameters (af0,af1,af2) */
	double tgd[4];      /* group delay parameters */
} eph_t;

GNSSPROCESS_API				int PositionDiff(GNSSDATA Obs, Orb_Clk  m_Orb, double BaseCoor[], double Corr[]);
GNSSPROCESS_API				int Processing(int nObs, obsd_t *Obs, int nEph, eph_t *eph,  int NumbSate, int ObserveSate[], double BaseCoor[], double RoveCoor[]);

#endif