#ifndef _SPPPROCESS_HEADER_20150126_
#define _SPPPROCESS_HEADER_20150126_
#include "../Common/SVEphemeris.h"
#include "../Common/GNSSERROR.h"


#ifdef GNSSPROCESS_EXPORTS
#define GNSSPROCESS_API _declspec(dllexport)
#else
#define GNSSPROCESS_API _declspec(dllimport)
#endif

struct SPPResult
{
	double x[7];
	double Q[49];
};

class SPP{
public:
	 GNSSPROCESS_API	          SPP();
     GNSSPROCESS_API     bool     SPPProcess(SPPResult &Result, const GNSSDATA Obs,const double SatXYZ[], const double SatClk[], double SatAzEl[], double ephRMS[]);
						 bool     EstPos(SPPResult &Result,const GNSSDATA Obs, const double SatXYZ[], const double SatClk[], double Azel[], double EphRMS[]);
						 u2		  ResCode(u2 Iter, double x[], const GNSSDATA Obs, const double SatXYZ[], const double SatClk[], 
										double Azel[], double EphRMS[], double H[], double V[],double Var[], u2 &Nx, u2 RecCLK[]);
						 double   Prange(const GNSSDATA Obs, u2 Numb, double Satel, u2 SysNum);
public:
	tec_tt      tec;
	SPPTYPE     m_Type;
};


#endif