#ifndef _BASEFUNCTION_HEADER_20150119_
#define _BASEFUNCTION_HEADER_20150119_
#include "../Common/Common.h"

#include <string.h>
#include <math.h>
#include <stdlib.h>

#ifdef BASEFUNCTION_EXPORTS
#define BASEFUNCTION_API _declspec(dllexport)
#else
#define BASEFUNCTION_API _declspec(dllimport)
#endif
BASEFUNCTION_API     char *trim(char *line);
BASEFUNCTION_API     char *getstr(char *out,char *line,int ini, int length);
BASEFUNCTION_API     double str2num(const char *s, int i, int n);
BASEFUNCTION_API     int xstrReplace (char szSrc[], char ch2Re, char chRe, unsigned int maxCh);
BASEFUNCTION_API     void xstrStretch (char* pszSrc, unsigned int lenEnd);
BASEFUNCTION_API     void xstrtrim (char* pszDest, char* pszSrc);
BASEFUNCTION_API     void xstrmid(char *szDest, char *szSrc, int nPos, int nCount);

//
BASEFUNCTION_API     int JudgeSys(char strSys, u2 NAVSYS, u2 *SYS);
BASEFUNCTION_API     u2  SatNo(u2 Sat, u2 NAVSYS);
BASEFUNCTION_API     u2 Prn2Sys(u2 Prn, u2 &SYS, u2 &SysNum, char &StrSys);
BASEFUNCTION_API     int Satid2No(const char *id);
//time
BASEFUNCTION_API     double Ttime2GPS(TTime T, int *week);
BASEFUNCTION_API     gtime_t timeadd(gtime_t t, double sec);
BASEFUNCTION_API     double timediff(gtime_t t1, gtime_t t2);
BASEFUNCTION_API     gtime_t epoch2time(const double *ep);
BASEFUNCTION_API     void time2epoch(gtime_t t, double *ep);
BASEFUNCTION_API     gtime_t gpst2time(int week, double sec);
BASEFUNCTION_API     double time2doy(gtime_t t);
BASEFUNCTION_API     double time2gpst(gtime_t t, int *week);
BASEFUNCTION_API	 int MJDN(struct tm *tm);
BASEFUNCTION_API     TTime cal2t (int year, int month, int day, int hour, int minute, double second);
BASEFUNCTION_API     void t2doy(TTime *t, int *year, double *doy);
BASEFUNCTION_API     void t2tm(TTime *t, struct tm *tm, double *seconds);
BASEFUNCTION_API     gtime_t gltime2rtktime(TTime t);
BASEFUNCTION_API     TTime rtktime2gltime(gtime_t t);
BASEFUNCTION_API     void TTimeToCommonTime(const TTime* pmjd, COMMONTIME* pct);
BASEFUNCTION_API     gtime_t utc2gpst(gtime_t t);
BASEFUNCTION_API     gtime_t gpst2utc(gtime_t t);
BASEFUNCTION_API     double utc2gmst(gtime_t t, double ut1_utc);
BASEFUNCTION_API     int   str2time(const char *s, int i, int n, gtime_t *t);

//coordinate
BASEFUNCTION_API     void ecef2pos(const double *r, double *pos);
BASEFUNCTION_API     void xyz2enu(const double *pos, double *E);
BASEFUNCTION_API     void ecef2enu(const double *pos, const double *r, double *e);

BASEFUNCTION_API     Coordinate operator %(Coordinate coor1, Coordinate coor2);
BASEFUNCTION_API     Coordinate unit(Coordinate coor);
BASEFUNCTION_API     double Dot(int n, const double v1[], const double v2[]);
BASEFUNCTION_API     void sunxyz(TTime& mjd, double rs[]);

BASEFUNCTION_API     double *zeros(int n, int m);
BASEFUNCTION_API	 int *imat(int n, int m);
BASEFUNCTION_API     double *mat(int n, int m);
BASEFUNCTION_API     double *eye(int n);
BASEFUNCTION_API     double dot(const double *a, const double *b, int n);
BASEFUNCTION_API     double norm(const double *a, int n);
BASEFUNCTION_API     int normv3(const double *a, double *b);
BASEFUNCTION_API     void cross3(const double *a, const double *b, double *c);
BASEFUNCTION_API     void sunmoonpos(gtime_t tutc, const double *erpv, double *rsun,
									 double *rmoon, double *gmst);

BASEFUNCTION_API     int matinv(double *A, int n);
BASEFUNCTION_API     void matcpy(double *A, const double *B, int n, int m);
BASEFUNCTION_API     void matmul(const char *tr, int n, int k, int m, double alpha,
			                     const double *A, const double *B, double beta, double *C);
BASEFUNCTION_API	 int lsq(const double *A, const double *y, int n, int m, double *x,
							 double *Q);
BASEFUNCTION_API     int filter(double *x, double *P, const double *H, const double *v,
								const double *R, int n, int m);


BASEFUNCTION_API     int  cmpWLInfo(const void *p1, const void *p2);
BASEFUNCTION_API     int  cmpobs(const void *p1, const void *p2);
BASEFUNCTION_API     void OutPutObs(FILE *fp_Obs, GNSSDATA ObsData, u2 NAVSYS);
BASEFUNCTION_API     void GetHighRateProducts(FILE *ft, gtime_t GNowTime, IGS_Orbit_t *Orb, IGS_Clk_t *Clk);
BASEFUNCTION_API     int  GetHighRateTime(FILE *ft, gtime_t &gtime, int type, char line1[MAXLENGTH]);

BASEFUNCTION_API    void AddWLInfo(AMB_Info &AmbInfo, int Prn, int SSod, int ESod);

BASEFUNCTION_API    double GetBetaAngle(double XSun[3], double SatXYZ[3], double VSat[3]);
BASEFUNCTION_API    double GetIFCBCorr(char chSys, int prn, int MJDN, double dBeta, double mu);
BASEFUNCTION_API    double GetIFCBCorr(double dCos[4], double dSin[4], double mu);
BASEFUNCTION_API    void RePutUPD();
#endif