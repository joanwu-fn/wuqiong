#ifndef _SVEPPHEMERIS_HEADER_20150125
#define _SVEPPHEMERIS_HEADER_20150125

#include "../Common/Common.h"

#ifdef EMPHERIS_EXPORTS
#define EMPHERIS_API _declspec(dllexport)
#else
#define EMPHERIS_API _declspec(dllimport)
#endif

#define NMAX    8
class  Orb_Clk
{
public:
	EMPHERIS_API Orb_Clk();
	EMPHERIS_API ~Orb_Clk();
	EMPHERIS_API u2   GetOrb_Clk(GNSSDATA &Obs, double RecClk, double SatXYZ[], double SatClk[], double SatAzEl[], double ephRMS[],double rb[3], double pos[3]);
	EMPHERIS_API u2   GetSatellitePos(GNSSDATA &Obs, int prn, double SatXYZ[], double SatClk[], double SatAzEl[], double ephRMS[], double *rb, double *pos);
	EMPHERIS_API bool GetOrb_Clk_Brdc(gtime_t time, u2 prn, double SatXYZ[], double &SatClk, double SatAzEl[], double &ephRMS);
	EMPHERIS_API bool GetOrb_Clk_Igs(gtime_t time, u2 prn, double SatXYZ[], double &SatClk, double SatAzEl[], double &ephRMS);
	EMPHERIS_API void GetBetaAngle(gtime_t time, double Beta[]);
	EMPHERIS_API void Initialize(int MJDN);
	             bool Get_Orb_Igs(gtime_t time, u2 prn, double SatXYZ[],double &ephRMS);
				 bool Get_Clk_Igs(gtime_t time, u2 prn, double &SatClk);
				 void GetSat_Pos(gtime_t time, u2 prn, double SatXYZ[], double &ephRMS, double &SatClk);
				 void satantoff(gtime_t time, const double *rs, Pcv_t *pcv,double *dant);


public:
	EPHN_t      m_eph[MAXSAT];
	u2          m_nclk;
	u2          m_nsp3;
	u2          m_SampleRate[2];
	IGS_Clk_t   m_clk[7];
	IGS_Orbit_t m_sp3[21];
	Pcv_t       m_Pcv[MAXSAT];
	u2          m_type;
	/*---------------------*/
	bool        m_HighRate;
	gtime_t     m_Nowtime;
	double      _dSin[MAXSAT][4];
	double      _dCos[MAXSAT][4];
};
#endif