#ifndef _PNETCONTROL_HEADER_20150129_
#define _PNETCONTROL_HEADER_20150129_
#include "../Common/Common.h"
#include "../Common/BaseFunction.h"
#include "../Common/debug.h"
#include "../Common/SVEphemeris.h"
#include "../Common/FileLoad.h"
#include "../GNSSProcess/SPPProcess.h"
#include "../GNSSProcess/PPPProcess.h"
#include "../GNSSProcess/PosDifferential.h"

//#define         HIGHRATEPRODUCTS
enum GNSSProcessType
{
	Process_SPP  = 1,
	Process_PPP  = 2,
	Process_RTK  = 3,
	Process_Net  = 4,
	Process_TEQC = 5,
	Process_FCB  = 6
};

#ifdef PNETCONTROL_EXPORTS
#define PNETCONTROL_API _declspec(dllexport)
#else
#define PNETCONTROL_API _declspec(dllimport)
#endif

//把SPP 和 PPP 写成不需要管理内存的dll，只提供一个接口， 参数由外部传入，输出并更新结果，借鉴四川CORS项目的设计

class Netcontrol
{
public:
	PNETCONTROL_API          Netcontrol();
    PNETCONTROL_API          ~Netcontrol();
	PNETCONTROL_API	   void  ReadCtl(char *ctrl_file);
	PNETCONTROL_API    void  IniNet();
	PNETCONTROL_API    void  LoadNav(int MJDN);
	PNETCONTROL_API    void  InitailFile(int MJDN);
	PNETCONTROL_API    void  Process(int MJDN);
	PNETCONTROL_API    void  ProcessCompareObs(int MJDN);
	PNETCONTROL_API    void  JudgeVisible(const GNSSDATA Obs, GNSSInfo &PInfo, double rb[3], double pos[3]);
	PNETCONTROL_API    void  Release();
	PNETCONTROL_API    void  OutPutUPDStatis();
					   void  UpDataEph(gtime_t time, TTime TNowTime);
					   void  UpdateHighRate(FILE *ft, gtime_t time);
					   void  OutPutMP();
					   void  OutPutDCB();
					   void  OutPutWLRinex(FILE *ft, int MJDN, int StaNo, AMB_Info &AInfo);

private:

public:
	/*  Function  */
	Orb_Clk       *m_OrbClk;            // compute satellite position 
	SPP             m_SPP;              // get the SPP result
	PPP             m_PPP;              // get the PPP result

    /* FILE  LOAD */
	RinexNavLoad    m_NavLoad;
	IGSOrbitLoad    m_OrbLoad;
	IGSClkLoad      m_ClkLoad;
	RinexObsLoad    *m_ObsLoad;

    /* Option and Flag */	
	ProcessOption   m_Option;
	bool            m_Initial;
	enum  GNSSProcessType m_ProcessType;
	/*Obs and eph memory*/
	GNSSDATA        m_ObsData;
	GNSSEPH_t       m_eph;
	GNSSInfo        *m_GNSSInfo;

	/*     OutPut    */
	FILE  *m_ftstatis;
	int   m_TTnumTEQC;

	UPDInfo         _UPDInfo;
	GNSSCodeBias    _bias;
};





#endif