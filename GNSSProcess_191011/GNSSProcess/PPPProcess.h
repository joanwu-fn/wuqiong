#ifndef _PPPPROCESS_HEADER_20150202_
#define _PPPPROCESS_HEADER_20150202_
#include "../GNSSProcess/SPPProcess.h"
#include "../PreProcess/PreProcess.h"

struct UDSatStaInfo 
{
	u2     SYS;
	u2     SysNum;
	double Geometry;
	double Azel[4];
	double e[3];
	double Dtrop;
	double Vtrop;
	double TropMFuncB[3];
	double TropMFuncR[3];
	double Diono;
	double Viono;
	double Observation[NEFREQ*2];
	double OMC[NEFREQ*2];
	double Weight[NEFREQ*2];
};
class PPP
{
public:
	GNSSPROCESS_API       		 PPP();
	GNSSPROCESS_API     bool     PPPProcess(GNSSInfo &PInfo, double X[], const GNSSDATA Obs,const double SatXYZ[], const double SatClk[],
												double SatAzEl[], double ephRMS[]);
	GNSSPROCESS_API     void     ProcessAmb(GNSSInfo &PInfo, GNSSDATA Obs, double *Azel, AMB_Info &AmbInfo);
						void     ComputeAR(GNSSDATA Obs, int Sample, AMB_Info &AmbInfo);
	GNSSPROCESS_API     int      ResPPP(u2 Iter, GNSSInfo &PInfo, const GNSSDATA Obs, const double SatXYZ[], const double SatClk[], 
											double Azel[], double EphRMS[], double H[], double V[], double R[]);
// 						double   GetMeasureMent(GNSSDATA Obs, u2 Numb, double Satel, u2 SysNum, u2 Type, u2 Fre);



public:
	double m_C1,m_C2,m_C3,m_C4;
	double _DCB13[14],_DCB12[14];
};
				  int  Get_SingleDiff_Res(u2 Iter, GNSSInfo &PInfo, const GNSSDATA Obs, const double SatXYZ[], const double SatClk[], 
											double Azel[], double EphRMS[], double H[], double V[], double R[]);
				  int Get_UnDiff_Res(u2 Iter, GNSSInfo &PInfo, const GNSSDATA Obs, const double SatXYZ[], const double SatClk[], 
											double Azel[], double EphRMS[], double H[], double V[], double R[]);
GNSSPROCESS_API  bool BaseLineProcess(GNSSInfo &BaseLineInfo, const GNSSDATA Obs,const double SatXYZ[], const double SatClk[], double SatAzEl[], double ephRMS[]);
GNSSPROCESS_API  double GetMeasureMent(GNSSDATA Obs, u2 Numb, u2 SysNum, u2 Type, u2 Fre);

#endif