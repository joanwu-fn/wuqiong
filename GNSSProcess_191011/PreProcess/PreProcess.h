#ifndef _PREPROCESS_HEADER_20150203
#define _PREPROCESS_HEADER_20150203
#include "../Common/BaseFunction.h"

#define  MAX_TEQC_EPOCH   (300)
#define FCB_DETAIAL

#ifdef PREPROCESS_EXPORTS
#define PREPROCESS_API _declspec(dllexport)
#else
#define PREPROCESS_API _declspec(dllimport)
#endif
     
const int MINAMB = 3;

enum OBS_LEVEL						{ OBS_REFER, OBS_LOWLEVEL, OBS_TOCHECK, OBS_NOCHECK, OBS_CONSTRAIN };
// X Y Z dtGPS dtBDS dtGAL dtGLO trop1 trop2 Ambi ...  ...
#define IT(r)     (11+r)
#define IB(sat, Fre)     (IT(1) + sat  + Fre * MAXSAT)

PREPROCESS_API     void  PreProcess(const GNSSDATA Obs, GNSSInfo &PInfo, double X[]);
PREPROCESS_API     void  TEQC(const GNSSDATA Obs, GNSSInfo &PInfo);
PREPROCESS_API	   void  EstDCB(const GNSSDATA &Obs, double *Azel, GNSSInfo &PInfo);
PREPROCESS_API     void  EstFCB(const GNSSDATA &Obs, double *Azel, GNSSInfo &PInfo, FILE* ftMW);
PREPROCESS_API     void  AddGFIF(const GNSSDATA &Obs, double *Azel, GNSSInfo &PInfo, FILE* ftMW);
PREPROCESS_API     void  AddOneFcb(GNSSInfo &PInfo, int Prn);
PREPROCESS_API     void  ComputeAverage(UPDInfo *UPDIf, GNSSInfo *PInfo, int StaNum, int Type, FILE *ft, int nSat[2]);
PREPROCESS_API     void  GetAverageDCB(GNSSInfo &PInfo);
PREPROCESS_API     void  AddLastTime(GNSSInfo &PInfo);
                  double MPValue(const ObsData_t ObsData, enum MPTYPE type);
                   void  GetMP(const GNSSDATA Obs, GNSSInfo &PInfo);
				   void  IniAmbiguity(const GNSSDATA Obs, GNSSInfo &PInfo);
PREPROCESS_API	   void  DetectCS_GF2(const GNSSDATA Obs, SSat_t  SSat[MAXSAT]);
PREPROCESS_API     void  DetectCS_MW12(const GNSSDATA Obs, SSat_t SSat[MAXSAT]);
PREPROCESS_API	   void  DetectCS_MW23(const GNSSDATA Obs, SSat_t SSat[MAXSAT]);
PREPROCESS_API     void  DetectCS_GF2_13(const GNSSDATA Obs, SSat_t  SSat[MAXSAT]);

/*---------------------------------gross error detection with medial value----------------------------------*/
class MedialDetection
{
public:
	 MedialDetection(int maxobs, float *stdCutoff = NULL);
	 ~MedialDetection();

	 void	Initialize();
	 void	UpdateList(float val, OBS_LEVEL obsLevel = OBS_REFER);
	 bool	ReWeight(int& nNormal, int minObs = 5, int enlarge = 1);
	 float	GetOMCWei(int ind);
	 float	GetValue(int ind);
	 int	GetNumObs();
	 float	GetAverage();
	 float	GetSTD();
	 float	GetPercent(int percent);
public:
	int			_nOmc;
	int			_nSort;
	int           _nObsCheck;
	int				_maxObs;
	float			_average;
	float			_std;
	float			_stdCutoff[2];
	float*			_omcList;
	float*			_omcSort;
	float			_meanBak;
	float			_varBak;
	float*			_omcWeight;
	OBS_LEVEL*		_obsLevel;
};


class UPDModel
{
public:
	void       InitialUPD(UPDInfo* UPDIf, GNSSInfo *PInfo, int StaNum, int Type, FILE *ft, int nSat[2]);
	void       OutPutUPD(char chData[255]);
	void       FindStation(int &iRefSite, bool bFlag);
    void       GetIniUPD();
	void       GetNextUPD();
	void       GetAverageUPD();
	bool       ExamineConsistence(MedialDetection* pMedia1, MedialDetection* pMedia2, int iStaNo, int& iMediaNum,int prnlist[MAXSAT]);
    void       ComputeStatistical(UPDInfo* UPDIf, GNSSInfo *PInfo, int StaNum, int Type, FILE *ft, int nSat[2]);
public:
	UPDInfo*   _pUPD;
	int        _nSite;
	GNSSInfo*  _PInfo;
	int        _iPrnNum[2];  //Start Prn / End Prn
	int        _iType;
	//==================================
	bool       _bIniSat[MAXSAT];
	double     _dSatUPD[MAXSAT];
	bool       _bIniSite[MAXSTATION];
	bool       _bSiteSlip[MAXSTATION];
    double     _dRcvUPD[MAXSTATION];
	int        _iIteration;
	FILE*      _ftOut;
};
#endif