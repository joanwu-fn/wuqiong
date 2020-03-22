/**/
#ifndef _FILELOAD_HEADER_20150118
#define _FILELOAD_HERDER_20150118
//Í·ÎÄ¼þ
#include "../Common/Common.h"
#include "../Common//BaseFunction.h"

#include <stdlib.h>


#ifdef FILELOAD_EXPORTS
#define FILELOAD_API _declspec(dllexport)
#else
#define FILELOAD_API _declspec(dllimport)
#endif

FILELOAD_API     int readngspcv(const char *file, Pcv_t *pcvs);
FILELOAD_API     int readantex(const char *file, Pcv_t *pcvs);
struct OBSTYPE
{
	OBSTYPE()
	{
		ntype = 0;
		memset(TypeStrList, 0x00, sizeof(char) * MAXENTRY_NUMBER * 4);
		memset(Frequncy,    0x00, sizeof(u2)   * MAXENTRY_NUMBER);
		memset(Type,        0x00, sizeof(u2)   * MAXENTRY_NUMBER);
		memset(valid,       0x00, sizeof(u2)   * MAXENTRY_NUMBER);
		memset(priority,    0x00, sizeof(u2)   * MAXENTRY_NUMBER);
		memset(bP1C1,       0x00, sizeof(bool) * MAXENTRY_NUMBER);
		memset(bP2C2,       0x00, sizeof(bool) * MAXENTRY_NUMBER);
	}
	u2   ntype;
	u2   Frequncy[MAXENTRY_NUMBER];  //123  : 012
	u2   Type[MAXENTRY_NUMBER];      //CLSD : 0123
	u2   valid[MAXENTRY_NUMBER];     //0 1
	u2   priority[MAXENTRY_NUMBER];  
	char TypeStrList[MAXENTRY_NUMBER][4];
	bool bP1C1[MAXENTRY_NUMBER];
	bool bP2C2[MAXENTRY_NUMBER];
};


class RinexObsLoad
{
public:
FILELOAD_API	RinexObsLoad();
FILELOAD_API	~RinexObsLoad();
FILELOAD_API	bool Initial(char chFile[], u2 NavSys, FILE *ft_out);
FILELOAD_API    void InitRinexload(GNSSCodeBias bias, int RecType);
FILELOAD_API	bool GetObsData(gtime_t time, GNSSDATA &ObsData);
private:
	bool GetObsData_3(gtime_t time, GNSSDATA &ObsData);
	bool GetObsData_2(gtime_t time, GNSSDATA &ObsData);
	bool GetObsType_2(char buff[]);
	bool GetObsType_3(char buff[]);
	unsigned char obs2code(const char *obs, int sys, u2 &freq);
	bool ReadObsHeader(FILE *ft_out);

public:

	u2        m_Version;               //Rinex Version
	u2        m_SampRate;
	FILE      *m_File;                 //FILE *

	u2        m_NavSys;                //Navigation System to use:   1: GPS    2: BDS     4: Galileo    8: GLONASS
    GNSSDATA  m_NowObsData;            //Obs data of now epoch
	OBSTYPE   m_ObsType[MAXSYS];
	u2        m_StaNo;
	GNSSCodeBias _bias;
	int       _RecType;
};


class RinexNavLoad{
public:
	FILELOAD_API	RinexNavLoad();
	FILELOAD_API	~RinexNavLoad();
	FILELOAD_API	bool Initial(char chFile[], u2 NavSys = SYS_GPS);
	                bool ReadNavHeader();
	FILELOAD_API    bool ReadNavData(GNSSEPH_t &brdc);
					bool ReadNavData_3(GNSSEPH_t &brdc);
					bool ReadNavData_2(GNSSEPH_t &brdc);
private:

public:
	u2        m_Version;
	u2        m_SYS;
	FILE      *m_File;                 //FILE *
	double    m_Alpha[4],m_Beta[4];
};

class IGSClkLoad
{
public:
	FILELOAD_API	IGSClkLoad();
	FILELOAD_API	~IGSClkLoad();
	FILELOAD_API	bool Initial(char chFile[]);
	                bool ReadClkHeader();
	FILELOAD_API    bool ReadClkData(GNSSIGS_t &IGS);
private:

public:
	FILE      *m_File; 
};

class IGSOrbitLoad
{
public:
	FILELOAD_API	IGSOrbitLoad();
	FILELOAD_API	~IGSOrbitLoad();
	FILELOAD_API	bool Initial(char chFile[]);
					bool ReadOrbitHeader();
	FILELOAD_API    bool ReadOrbitData(GNSSIGS_t &IGS);
private:

public:
	FILE      *m_File;
};

FILELOAD_API  bool ReadTecFile(char *file, tec_tt * tec);
//FILELOAD_API  bool ReadDCB(char *dcbfile, double DCBValue[MAXSAT]);
#endif