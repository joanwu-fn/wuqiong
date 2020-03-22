#ifndef _GNSSERROR_HEADER_20150126_
#define _GNSSERROR_HEADER_20150126_
#include "../Common//BaseFunction.h"
#include "../Common/Common.h"

BASEFUNCTION_API     double SD_var(double var, double el);
BASEFUNCTION_API     double var_LC(double i, double j, double k, double sig,int sys);
BASEFUNCTION_API     double varerr(double el, int type, u2 sys);
BASEFUNCTION_API     void   tidedisp(gtime_t tutc, const double *rr, int opt, const erp_t *erp,const double *odisp, double *dr);
BASEFUNCTION_API     double prectrop(gtime_t time, double pos[], double azel[],double x[], double dtdx[], double &var);
BASEFUNCTION_API     void   GetTropDelay(gtime_t time,double pos[], double azel[],double &Trop, double &Var);
BASEFUNCTION_API     int iontec(gtime_t time, const double *pos, tec_t *tec, int nt,
		                         const double *azel, int opt, double *delay, double *var);
BASEFUNCTION_API     void   GetBDSBias(double Satel, u2 Prn, u2 Fre, double &Bias);
BASEFUNCTION_API     void   GetBDSCodeBiasVariations(int prn, double Elevation, u2 Fre, double &bias);
BASEFUNCTION_API     double satazel(const double *pos, const double *e, double *azel);
BASEFUNCTION_API     double geodist(const double *rs, const double *rr, double *e);
BASEFUNCTION_API	 int    valsol(const double *azel, const double *v, int nv, int nx);
BASEFUNCTION_API	 double  MWObs(const GNSSDATA Obs, u2 Numb, double el, int i, int j);
BASEFUNCTION_API     double  GetCombination(const GNSSDATA Obs, u2 Numb, double el, double Pse[3], double Carr[3]);
BASEFUNCTION_API     bool DetectRobust(double *Val, int number, double std, double &Ave, double &stdval);
BASEFUNCTION_API     bool GetMWObs(const GNSSDATA Obs, u2 Numb, double &MW, double Coeff[5], 
									int Type, double Lamda[3]);
#endif