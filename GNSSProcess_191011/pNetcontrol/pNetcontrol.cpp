#include "pNetcontrol.h"
#include <Windows.h>
#include <direct.h>

bool ReadDCB(char *dcbfile, double DCBValue[MAXSAT])
{
	memset(DCBValue, 0x00, sizeof(DCBValue));
	FILE* ft = NULL;
	ft = fopen(dcbfile, "r");
	if(ft == NULL)
	{
		return false;
	}
	char buff[MAXLENGTH];
	while (fgets(buff,MAXLENGTH, ft) != NULL)
	{
		if(buff[0] != '*')
		{
			continue;
		}
		while (fgets(buff,MAXLENGTH, ft) != NULL)
		{
			char chSat[4];
			double DCB, DCBStd;
			u2 Sys_Sys;
			sscanf(buff, "%s %lf %lf",chSat, &DCB, &DCBStd);
			if((JudgeSys(buff[0], SYS_ALL, &Sys_Sys)) == -1)
			{
				continue;
			}
			int prn = (int)str2num(buff,1,2);
			prn     = SatNo(prn, Sys_Sys);
			if(prn < 1 || prn > MAXSAT)
			{
				continue;
			}
			DCBValue[prn - 1] = DCB * 0.2998;
		}
	}
	fclose(ft);
	return true;
}
void ReadBDSCodeBias(char chBias[], GNSSCodeBias* bias)
{
	memset(bias->BDSRecBias, 0x00, sizeof(bias->BDSRecBias));
	bias->nRecType = 0;
	FILE* fCtrl = NULL;
	if ((fCtrl = fopen(chBias, "r")) == NULL)
	{
		return;
	}
	char chLineTemp[1000];
	while (fgets(chLineTemp, sizeof(chLineTemp), fCtrl))
	{
		if (chLineTemp[0] == '+')
		{
			int iFre = 0;
			while (fgets(chLineTemp, sizeof(chLineTemp), fCtrl))
			{
				if (chLineTemp[0] == '-')
				{
					break;
				}
				if(iFre > 2)
				{
					printf("!!!!!Code Bias Wrong!!!!\n");
					continue;
				}
				for (int prn_1 = 0; prn_1 < 14; prn_1++)
				{
					bias->BDSRecBias[bias->nRecType][NSATGPS + prn_1][iFre] = str2num(chLineTemp, 3 + prn_1 * 6, 6) * 0.2998;
				}
				iFre++;
			}
			bias->nRecType++;
		}
	}
	fclose(fCtrl);
}
/*******************************************************************/
/*   函数名：                 FormCommonObs                        */
/* 功能简介：     选择两个站（按卫星号排列）的共同观测值           */
/* 输入变量：    GNSSDATA &RoveObs,     GNSSDATA &BaseObs          */
/*     单位：            无                     无                 */
/* 输出变量：    GNSSDATA &RoveObs,     GNSSDATA &BaseObs          */
/*     单位：            无                     无                 */
/*   返回值：                      无                              */
/*     含义：                      无                              */
/*******************************************************************/
void FormCommonObs(GNSSDATA &RoveObs, GNSSDATA &BaseObs)
{
	int i,j,Common=0;
	for(i = 0; i < BaseObs.NumSats; i++)
	{
		for(j = Common; j < RoveObs.NumSats; j++)
		{
			if(BaseObs.ObsData[i].prn == RoveObs.ObsData[j].prn)
			{
				memcpy(&BaseObs.ObsData[Common], &BaseObs.ObsData[i], sizeof(ObsData_t));
				memcpy(&RoveObs.ObsData[Common], &RoveObs.ObsData[j], sizeof(ObsData_t));
				Common++;
			}
		}
	}
	BaseObs.NumSats = Common;
	RoveObs.NumSats = Common;
}

void ComputeDCBStatistical(GNSSInfo *PInfo, int StaNum, int type)
{
	int StaNo,prn_1,j;

	/*+++++++++++++++++++++++++++++++++++++++++++++++*/
	for(StaNo = 1; StaNo < StaNum; StaNo++)
	{
		int RefPrn_1 = -1;
		for(prn_1 = NSATGPS; prn_1 < NSATGPS + NSATBDS; prn_1++)
		{
			if(PInfo[StaNo].UDFcbFlag[prn_1][type] && PInfo[0].UDFcbFlag[prn_1][type])
			{
				RefPrn_1 = prn_1;
				break;
			}
		}
		if(RefPrn_1 != -1)
		{
			double RefDif = PInfo[StaNo].UDFcb[RefPrn_1][type] - PInfo[0].UDFcb[RefPrn_1][type];
			for(prn_1 = NSATGPS; prn_1 < NSATGPS + NSATBDS; prn_1++)
			{
				if(PInfo[StaNo].UDFcbFlag[prn_1][type])
				{
					PInfo[StaNo].UDFcb[prn_1][type] = PInfo[StaNo].UDFcb[prn_1][type] - RefDif;
				}
			}

			double deta[2] = {0};
			for(prn_1 = NSATGPS; prn_1 < NSATGPS + NSATBDS; prn_1++)
			{
				if(PInfo[StaNo].UDFcbFlag[prn_1][type])
				{
					if(PInfo[0].UDFcbFlag[prn_1][type])
					{
						double Fcbtemp = (PInfo[StaNo].UDFcb[prn_1][type] - PInfo[0].UDFcb[prn_1][type]);				
						deta[0] += Fcbtemp; 
						deta[1] += 1;
					}
				}
			}
			/*-------------------减去整体基准偏差---------------------*/
			if(deta[1] > 0)
			{
				deta[0] = deta[0] / deta[1];
				for(prn_1 = NSATGPS; prn_1 < NSATGPS + NSATBDS; prn_1++)
				{
					if(PInfo[StaNo].UDFcbFlag[prn_1][type])
					{
						PInfo[StaNo].UDFcb[prn_1][type] = PInfo[StaNo].UDFcb[prn_1][type] - deta[0];
					}
				}
			}
			/*-------------------减去整体基准偏差---------------------*/
		}		
	}
	/*+++++++++++++++++++++++++++++++++++++++++++++++*/
	return;
}
void FormSDObs(GNSSDATA &RoveObs, const GNSSDATA BaseObs)
{
	int j;
	for(int i = 0; i < RoveObs.NumSats; i++)
	{
		for(j = 0; j < 3; j++)
		{
			if(RoveObs.ObsData[i].P[j] > PRECISION && BaseObs.ObsData[i].P[j] > PRECISION && 
				fabs(RoveObs.ObsData[i].L[j]) > PRECISION && fabs(BaseObs.ObsData[i].L[j]) > PRECISION)
			{
				RoveObs.ObsData[i].P[j] -= BaseObs.ObsData[i].P[j];
				RoveObs.ObsData[i].L[j] -= BaseObs.ObsData[i].L[j];
			}
			else
			{
				RoveObs.ObsData[i].P[j] = 0.0;
				RoveObs.ObsData[i].L[j] = 0.0;
			}
		}	
	}
}
void InitialGNSSIfo(GNSSInfo & PInfo, enum GNSSProcessType ProcessType, int TTnumb)
{
	PInfo.Nx = 3 + 4 +1 +3 * MAXSAT;
	
	PInfo.Mode = PMODE_PPP_KINEMA;
	if(PInfo.m_X == NULL)
	{
		PInfo.m_X = zeros(PInfo.Nx, 1);
	}
	if(PInfo.m_P == NULL)
	{
		PInfo.m_P = zeros(PInfo.Nx, PInfo.Nx);
	}
	PInfo.Numb = 0;
	PInfo.Nfcb[0] = PInfo.Nfcb[1] = 0;
	PInfo.Fcb = NULL;


	memset(PInfo.bUPDExclude, 0x00, sizeof(bool) * MAXSAT * 2);
	memset(PInfo.UDFcbFlag,   0x00, sizeof(bool) * MAXSAT * 2);
	memset(PInfo.UDFcb,       0x00, sizeof(double) * MAXSAT * 2);

	memset(&PInfo.RMS, 0x00, sizeof(double) * 3);
	memset(PInfo.SSat, 0x00, sizeof(SSat_t)*MAXSAT);
	memset(&PInfo.LastObsTime, 0x00, sizeof(gtime_t));

	memset(PInfo.dPIFGF, 0x00, sizeof(double) * MAXSAT);
	memset(PInfo.iPIFGF, 0x00, sizeof(int) * MAXSAT);

	PInfo.UseFulNum = 0;
	PInfo.EpoNum = 0;
	for(int i = 0; i < MAXSAT; i++)
	{
		PInfo.MP_S[i].LastEpoch[0] = PInfo.MP_S[i].LastEpoch[1] = 0;
		memset(PInfo.MP_S[i].Std, 0x00, sizeof(double) * 8);
	}

	for(int i = 0; i < 3; i++)
	{
		PInfo.RMS[i] = PInfo.Ave[i] = 0.0;
	}
	memset(&PInfo.Visible,  0x00, sizeof(int) * MAXSAT * 2);
	return;
}

Netcontrol::Netcontrol()
{
	m_ProcessType = Process_SPP;
	m_ObsLoad  = NULL;
	m_OrbClk   = NULL;
	m_GNSSInfo  = NULL;
	memset(_UPDInfo._EWLRefFCB, 0x00, sizeof(_UPDInfo._EWLRefFCB));
	memset(_UPDInfo._WLRefFCB,  0x00, sizeof(_UPDInfo._WLRefFCB));
	memset(_UPDInfo._bEWLRef,   0x00, sizeof(_UPDInfo._bEWLRef));
	memset(_UPDInfo._bWLRef,    0x00, sizeof(_UPDInfo._bWLRef));

	_UPDInfo._dEWLSeri = new double[MAXSAT * 365];
	_UPDInfo._dWLSeri  = new double[MAXSAT * 365];
	_UPDInfo._bEWLSeri = new bool[MAXSAT * 365];
	_UPDInfo._bWLSeri  = new bool[MAXSAT * 365];
    _UPDInfo._ftRes    = fopen("UPDRes.txt", "w");
	memset(_UPDInfo._dEWLSeri, 0x00, sizeof(double) * MAXSAT * 365);
	memset(_UPDInfo._dWLSeri,  0x00, sizeof(double) * MAXSAT * 365);
	memset(_UPDInfo._bEWLSeri, 0x00, sizeof(bool)   * MAXSAT * 365);
	memset(_UPDInfo._bWLSeri,  0x00, sizeof(bool)   * MAXSAT * 365);
	_UPDInfo._nDay = 0;

	memset(&_bias, 0x00, sizeof(GNSSCodeBias));
}
void Netcontrol::Release()
{
	int i,j;
	for(i = 0; i < m_Option.StaNum; i++)
	{
		if(m_ObsLoad[i].m_File != NULL)
		{
			fclose(m_ObsLoad[i].m_File);
		}
	}

	for(i = 0; i < MAXSYS; i++)
	{
		if(m_eph.brdc[i].max > 0)
		{
			free(m_eph.brdc[i].eph);
			m_eph.brdc[i].eph = NULL;
			m_eph.brdc[i].n = m_eph.brdc[i].max = 0;
		}
	}
	if(m_eph.IGS.maxsp3 > 0)
	{
		free(m_eph.IGS.sp3);
		m_eph.IGS.sp3 = NULL;
		m_eph.IGS.maxsp3 = m_eph.IGS.nsp3 = 0;
	}
	if(m_eph.IGS.maxclk > 0)
	{
		free(m_eph.IGS.clk);
		m_eph.IGS.clk = NULL;
		m_eph.IGS.maxclk = m_eph.IGS.nclk = 0;
	}
	free(m_OrbClk);
	m_OrbClk = NULL;

	delete []m_ObsLoad;
	m_ObsLoad = NULL;
}

Netcontrol::~Netcontrol()
{
	for(int i = 0; i < m_Option.StaNum; i++)
	{
		free(m_GNSSInfo[i].m_X);
		m_GNSSInfo[i].m_X = NULL;
		free(m_GNSSInfo[i].m_P);
		m_GNSSInfo[i].m_P = NULL;
	}
	delete []m_GNSSInfo;
	delete []_UPDInfo._dEWLSeri; 
	delete []_UPDInfo._dWLSeri;
	delete []_UPDInfo._bEWLSeri;
	delete []_UPDInfo._bWLSeri;
	fclose(_UPDInfo._ftRes);
}

void Netcontrol::IniNet()
{
	if(m_Option.StaNum > 0)
	{
		m_GNSSInfo = new GNSSInfo[m_Option.StaNum];
		for(int i = 0; i < m_Option.StaNum; i++)
		{
			m_GNSSInfo[i].m_X = NULL;
			m_GNSSInfo[i].m_P = NULL;
			m_GNSSInfo[i].Station = NULL;
			InitialGNSSIfo(m_GNSSInfo[i], m_ProcessType, m_TTnumTEQC);
			m_GNSSInfo[i].Station = m_Option.Station;
		}
	}
	return;
}
void Netcontrol::InitailFile(int MJDN)
{
	int i,j;
	char ObsFile[MAXLENGTH],OutPut[MAXLENGTH];
	ASSERT(m_Option.StaNum < MAXSTATION);

	int year,week,Day;
	double Doy,SoW;
	TTime Time;
	COMMONTIME ct;

	Time.MJDN = MJDN;
	Time.SoD = 0;

	TTimeToCommonTime(&Time, &ct);
	t2doy(&Time,&year,&Doy);
	SoW = Ttime2GPS(Time, &week);
	Day = (int)(SoW / 86400);

	m_Option.ProcessDoY = (int)Doy;
	m_Option.ProcessYear = year;

	if(m_Option.StaNum > 0)
	{	
		/********************************************************************************/
		char DCBFile[MAXLENGTH];
		sprintf(DCBFile,"%s/P1C1%02d%02d.DCB",m_Option.Directory.DCBDir,year%100,ct.month);
		if(!ReadDCB(DCBFile, _bias.P1C1))
		{
			printf("Lack of %s file\n", DCBFile);
			Sleep(1000);
		}
		sprintf(DCBFile,"%s/P2C2%02d%02d.DCB",m_Option.Directory.DCBDir,year%100,ct.month);
		if(!ReadDCB(DCBFile, _bias.P2C2))
		{
			printf("Lack of %s file\n", DCBFile);
			Sleep(1000);
		}
		sprintf(DCBFile,"%s/P3C3%02d%02d.DCB",m_Option.Directory.DCBDir,year%100,ct.month);
		if(!ReadDCB(DCBFile, _bias.P3C3))
		{
			printf("Lack of %s file\n", DCBFile);
			Sleep(1000);
		}
		memset(DCBFile, 0x00, sizeof(DCBFile));
		sprintf(DCBFile,"%s/CodeBias.txt",m_Option.Directory.TableDir);
		ReadBDSCodeBias(DCBFile, &_bias);

		m_ObsLoad  = new RinexObsLoad[m_Option.StaNum];
		for(i = 0; i < m_Option.StaNum; i++)
		{
			m_ObsLoad[i].InitRinexload(_bias, m_Option.Station[i].RecType);
			memset(&ObsFile, 0x00, sizeof(char) * MAXLENGTH);
			sprintf(ObsFile,"%s\\%04d\\%03d\\%s%03d0.%02do",m_Option.Directory.ObsDir,year,(int)Doy,m_Option.Station[i].Name,(int)Doy,year%100);
			m_ObsLoad[i].Initial(ObsFile, m_Option.NavSys, NULL);
			m_ObsLoad[i].m_StaNo = i;
		}
	
		/********************************************************************************/
		memset(&OutPut, 0x00, sizeof(char)*MAXLENGTH);
		sprintf(OutPut,"%s\\endoutput\\%03d",m_Option.Directory.ProDir,(int)Doy);
		mkdir(OutPut);
		memset(&OutPut, 0x00, sizeof(char)*MAXLENGTH);
		sprintf(OutPut,"%s\\log\\%03d",m_Option.Directory.ProDir,(int)Doy);
		mkdir(OutPut);
		/*                        LOAD  FILE                   */
		for(i = 0; i < m_Option.StaNum; i++)
		{
			ecef2pos(m_Option.Station[i].XYZ, m_Option.Station[i].Pos);

			memset(&OutPut, 0x00, sizeof(char)*MAXLENGTH);
			sprintf(OutPut,"%s\\endoutput\\%03d\\%s%03d0.%02dmw",m_Option.Directory.ProDir,(int)Doy,m_Option.Station[i].Name,(int)Doy,year%100);
			if(m_Option.Station[i].ft_mw != NULL)
			{
				fclose(m_Option.Station[i].ft_mw);
				m_Option.Station[i].ft_mw = NULL;
			}
			m_Option.Station[i].ft_mw = fopen(OutPut, "w");
			ASSERT(m_Option.Station[i].ft_mw != NULL);		
		}

		m_OrbClk = (Orb_Clk *)realloc(m_OrbClk, sizeof(Orb_Clk));
		memset(m_OrbClk->m_eph, 0x00, sizeof(EPHN_t)*MAXSAT);
		memset(m_OrbClk->m_clk, 0x00, sizeof(IGS_Clk_t)*10);
		memset(m_OrbClk->m_sp3, 0x00, sizeof(IGS_Orbit_t)*10);
		memset(m_OrbClk->m_Pcv, 0x00, sizeof(Pcv_t)*MAXSAT);
		m_OrbClk->m_nsp3 = m_OrbClk->m_nclk = 0;
		memset(&m_OrbClk->m_Nowtime, 0x00, sizeof(gtime_t));

		if(!m_Option.HighRate)
		{
			memset(&ObsFile, 0x00, sizeof(char)*MAXLENGTH);
			sprintf(ObsFile,"%s\\%04d\\brdm%03d0.%02dp",m_Option.Directory.NavDir,year,(int)Doy,year%100);
// 			sprintf(ObsFile,"%s\\%04d\\brdc%03d0.%02dn",m_Option.Directory.NavDir,year,(int)Doy,year%100);
			m_NavLoad.Initial(ObsFile, SYS_NONE);

			memset(ObsFile, 0x00, sizeof(char)*MAXLENGTH);
			sprintf(ObsFile,"%s\\%04d\\gbm%04d%d.sp3",m_Option.Directory.Sp3Dir,year,week,Day);
			m_OrbLoad.Initial(ObsFile);

			memset(ObsFile, 0x00, sizeof(char)*MAXLENGTH);
			sprintf(ObsFile,"%s\\%04d\\gbm%04d%d.clk",m_Option.Directory.Sp3Dir,year,week,Day);
			m_ClkLoad.Initial(ObsFile);
		}

		memset(ObsFile, 0x00, sizeof(char)*MAXLENGTH);
		sprintf(ObsFile,"%s\\antenna.pcv",m_Option.Directory.TableDir);
		readantex(ObsFile, m_OrbClk->m_Pcv);
		
		m_Initial = true;
		return;
	}
	else
	{
		return;
	}
}
void Netcontrol::ReadCtl(char *ctrl_file)
{
	m_Initial = false;
	FILE *fCtrl;
	if((fCtrl = fopen(ctrl_file, "r")) == NULL)
	{
		return;
	}
	M_Station *StationTemp;
	char chLineTemp[MAXLENGTH],chRemark[MAXLENGTH],chData[MAXLENGTH];
	while(fgets(chLineTemp, MAXLENGTH, fCtrl) != NULL)
	{
		if('#' == chLineTemp[0])
		{
			continue;
		}
		xstrStretch(chLineTemp, 120);
		xstrtrim(chLineTemp, chLineTemp);
		if (!strcmp(chLineTemp, "+Process mode section"))
		{
			/*-----------------Process mode section-----------------*/
			while(fgets(chLineTemp, MAXLENGTH, fCtrl) != NULL)
			{
				if('#' == chLineTemp[0])
				{
					continue;
				}
				xstrmid(chRemark, chLineTemp, 0, 20);
				xstrtrim(chRemark, chRemark);
				xstrmid(chData, chLineTemp, 20, 60);
				xstrtrim(chData, chData);
				xstrtrim(chLineTemp, chLineTemp);
				if (!strcmp(chRemark, "Mode"))
				{
					m_Option.Mode = atoi(chData);
					if(m_Option.Mode < 0)
					{
						m_Option.Mode = 0;
					}
					else if(m_Option.Mode > 1)
					{
						m_Option.Mode =1;
					}
				}
				else if (!strcmp(chLineTemp, "-Process mode section"))
				{
					break;
				}
			}	
		   /*-----------------Process mode section-----------------*/
		}
		else if (!strcmp(chLineTemp, "+Session time section"))
		{
			if (1 == m_Option.Mode)
			{
				/*-----------------Session time section-----------------*/
				while (fgets(chLineTemp, MAXLENGTH, fCtrl))
				{
					if ('#' == chLineTemp[0])
					{
						continue;
					}
					xstrmid(chRemark, chLineTemp, 0, 20);
					xstrtrim(chRemark, chRemark);
					xstrmid(chData, chLineTemp, 20, 60);
					xstrtrim(chData, chData);
					xstrtrim(chLineTemp, chLineTemp);
					int ep[6];
					if (!strcmp(chRemark, "From"))
					{
						sscanf(chData, "%d %d %d", &ep[0], &ep[1], &ep[2]);
						m_Option.StartTime = cal2t(ep[0],ep[1],ep[2], 0, 0, 0.0);
					}
					else if (!strcmp(chRemark, "To"))
					{
						sscanf(chData, "%d %d %d",  &ep[0], &ep[1], &ep[2]);
						m_Option.EndTime = cal2t(ep[0],ep[1],ep[2], 0, 0, 0.0);
					}
					else if (!strcmp(chLineTemp, "-Session time section"))
					{
						break;
					}
				}
			}
			/*-----------------Session time section-----------------*/
		}
		else if (!strcmp(chLineTemp, "+Project directory/server section"))
		{
			/*-----------Project directory/server section-----------*/
			while (fgets(chLineTemp, MAXLENGTH, fCtrl))
			{
				if ('#' == chLineTemp[0])
				{
					continue;
				}
				xstrmid(chRemark, chLineTemp, 0, 20);
				xstrtrim(chRemark, chRemark);
				xstrmid(chData, chLineTemp, 20, 60);
				xstrtrim(chData, chData);
				xstrtrim(chLineTemp, chLineTemp);
				if (!strcmp(chRemark, "Project"))
				{
					sprintf(m_Option.Directory.ProDir, chData);
				}
				else if (!strcmp(chRemark, "Mass product"))
				{
					;
				}
				else if (!strcmp(chRemark, "Rinex"))
				{
					sprintf(m_Option.Directory.ObsDir, "%s",chData);
				}
				else if (!strcmp(chRemark, "RinexRove"))
				{
					sprintf(m_Option.Directory.ObsRefDir, "%s",chData);
				}
				else if (!strcmp(chRemark, "Table"))
				{
					sprintf(m_Option.Directory.TableDir, "%s",chData);
				}
				else if (!strcmp(chRemark, "Mass"))
				{
					sprintf(m_Option.Directory.MassDir, "%s",chData);
				}
				else if(!strcmp(chRemark, "Iono"))
				{
					sprintf(m_Option.Directory.TECDir, "%s", chData);
				}
				else if (!strcmp(chRemark, "Dcb"))
				{
					sprintf(m_Option.Directory.DCBDir, "%s",chData);
				}
				else if (!strcmp(chRemark, "Navigation"))
				{
					sprintf(m_Option.Directory.NavDir, chData);
				}
				else if (!strcmp(chRemark, "Igs orbit"))
				{
					sprintf(m_Option.Directory.Sp3Dir, chData);
				}
				else if (!strcmp(chLineTemp, "-Project directory/server section"))
				{
					break;
				}
			}
			/*-----------Project directory/server section-----------*/
		}
		else if (!strcmp(chLineTemp, "+Solution strategy section"))
		{
			/*-----------Project directory/server section-----------*/
			while (fgets(chLineTemp, MAXLENGTH, fCtrl))
			{
				if ('#' == chLineTemp[0])
				{
					continue;
				}
				xstrmid(chRemark, chLineTemp, 0, 20);
				xstrtrim(chRemark, chRemark);
				xstrmid(chData, chLineTemp, 20, 60);
				xstrtrim(chData, chData);
				xstrtrim(chLineTemp, chLineTemp);
				if (!strcmp(chRemark, "process interval"))
				{
					m_Option.Interval = atof(chData);
					m_TTnumTEQC = 86400 / (int)m_Option.Interval;
				}
				else if(!strcmp(chRemark, "Using High Rtae"))
				{
					if(!strcmp(chData, "YES"))
					{
						m_Option.HighRate = true;
					}
					else
					{
						m_Option.HighRate = false;
					}
				}
				else if(!strcmp(chRemark, "process type"))
				{
					if(!strcmp(chData, "SPP"))
					{
						m_ProcessType = Process_SPP;
					}
					else if(!strcmp(chData, "PPP"))
					{
						m_ProcessType = Process_PPP;
					}
					else if(!strcmp(chData, "RTK"))
					{
						m_ProcessType = Process_RTK;
					}
					else if(!strcmp(chData, "NET"))
					{
						m_ProcessType = Process_Net;
					}
					else if(!strcmp(chData, "TEQC"))
					{
						m_ProcessType = Process_TEQC;
					}
					else if(!strcmp(chData, "FCB"))
					{
						m_ProcessType = Process_FCB;
					}
				}
				else if(!strcmp(chRemark, "nav sys to process"))
				{
					if(chData[0] == 'G')
					{
						m_Option.NavSys = SYS_GPS;
					}
					else if(chData[0] == 'C')
					{
						m_Option.NavSys = SYS_BDS;
					}
					else if(chData[0] == 'E')
					{
						m_Option.NavSys = SYS_GAL;
					}
					else if(chData[0] == 'M')
					{
						m_Option.NavSys = SYS_GPS | SYS_BDS | SYS_GAL;
					}
					else
					{
						m_Option.NavSys = SYS_NONE;
					}
				}
				else if (!strcmp(chLineTemp, "-Solution strategy section"))
				{
					break;
				}
			}
			/*-----------Project directory/server section-----------*/
		}
		else if (!strcmp(chLineTemp, "+Site section"))
		{
			int nSite = 0;
			char dyn = 's';
			int itemp;
			double dtemp;
			/*-----------+Site section-----------*/
			while (fgets(chLineTemp, MAXLENGTH, fCtrl))
			{
				xstrtrim(chData, chLineTemp);
				if ('#' == chLineTemp[0])
				{
					continue;
				}
				else if(!strcmp(chData, "-Site section"))
				{
					break;
				}
				if(nSite >= m_Option.StaNum)
				{
					m_Option.StaNum = nSite + 10;
					if ((StationTemp = (M_Station *)realloc(m_Option.Station, sizeof(M_Station)*m_Option.StaNum)) == NULL)
					{
						ASSERT(10 < 0);
					}
					m_Option.Station = StationTemp;
					for(int i = nSite; i < m_Option.StaNum; i++)
					{
						m_Option.Station[i].ft_ewl    = NULL;
						m_Option.Station[i].ft_mw     = NULL;
						m_Option.Station[i].ft_Result = NULL;
					}
				}
				sscanf(chData,"%s %c %d %d %lf %lf %lf %lf %lf %d %d %lf %lf %lf %d",
					&m_Option.Station[nSite].Name, &dyn, &itemp, &itemp, 
					&dtemp, &dtemp, &dtemp, &dtemp, &dtemp, &itemp, &itemp,
					&m_Option.Station[nSite].XYZ[0], &m_Option.Station[nSite].XYZ[1], &m_Option.Station[nSite].XYZ[2],
				    &m_Option.Station[nSite].RecType);
				if(dyn == 's' || dyn == 'S')
				{
					m_Option.Station[nSite].Pro_Mode = PMODE_PPP_STATIC;
				}
				else if(dyn == 'k' || dyn == 'K')
				{
					m_Option.Station[nSite].Pro_Mode = PMODE_PPP_KINEMA;
				}
				nSite++;
			}
			m_Option.StaNum = nSite;
		  /*-----------+Site section-----------*/
		}
	}
	fclose(fCtrl);
}
void Netcontrol::JudgeVisible(const GNSSDATA Obs, GNSSInfo &PInfo, double rb[3], double pos[3])
{
	int prn_1,j,PrnNo = 0;
	bool Flag;
	double SatXYZ[4],SatAzEl[2],ephRMS[2],e[3];
	for(prn_1 = 0; prn_1 < MAXSAT; prn_1++)
	{
		if(0 == m_OrbClk->m_type)
		{
			Flag = m_OrbClk->GetOrb_Clk_Brdc(Obs.ObsTime, prn_1+1, &SatXYZ[0], SatXYZ[3], &SatAzEl[0], ephRMS[0]);
		}
		else
		{
			Flag = m_OrbClk->GetOrb_Clk_Igs(Obs.ObsTime, prn_1+1, &SatXYZ[0], SatXYZ[3], &SatAzEl[0], ephRMS[0]);
		}
		if(false == Flag)
		{
			continue;
		}
		geodist(&SatXYZ[0], rb, e);
		satazel(pos, e, &SatAzEl[0]);
		if(SatAzEl[1] > 15 * D2R)
		{
			PInfo.Visible[prn_1][0]++;
			for(j = PrnNo; j < Obs.NumSats; j++)
			{
				if(prn_1 == (Obs.ObsData[j].prn - 1))
				{
					PrnNo = j + 1;
					PInfo.Visible[prn_1][1]++;
					break;
				}
			}
		}
		
	}
}

void Netcontrol::UpdateHighRate(FILE *ft, gtime_t time)
{
	char buff[MAXLENGTH] = "\0";
	if(m_OrbClk->m_Nowtime.time < 1000)  //first time
	{
		GetHighRateTime(ft, m_OrbClk->m_Nowtime, 1, buff);
	}
	else 
	{
		double diff = timediff(time, m_OrbClk->m_Nowtime);
		while(diff > 0)
		{
			if(diff < -1.5)
			{
				GetHighRateProducts(ft, m_OrbClk->m_Nowtime, &m_OrbClk->m_sp3[0], &m_OrbClk->m_clk[0]);
				return;
			}
			else
			{
				int ret = GetHighRateTime(ft, m_OrbClk->m_Nowtime, 1, buff);
				if(ret == -1)
				{
					return;
				}
				diff = timediff(time, m_OrbClk->m_Nowtime);
			}
		}		
	}
	
}
void Netcontrol::UpDataEph(gtime_t time, TTime TNowTime)
{
	int i,j;
	u2 Prn,SYS,SysNum;
	char StrSys;
	double diff;

	if(((int)TNowTime.SoD) % 3600 <= 40)  //check whether need to update brdc
	{
// 		printf("%10d  \n",(int)TNowTime.SoD);
		double BrdcTime[MAXSYS] = {14401.0, 14401.0, 14401.0, 14401.0};
		double BrdcTime1[MAXSYS] = {7201.0, 3601.0, 3601.0, 3601.0};
		for(Prn = 0; Prn < MAXSAT; Prn++)
		{
			if(Prn2Sys(Prn+1, SYS, SysNum, StrSys) <= 0)
			{
				continue;
			}
			for(i = 0; i < m_eph.brdc[SysNum].n; i++)
			{
				if( m_eph.brdc[SysNum].eph[i].prn == (Prn+1))
				{
					diff = timediff(time, m_eph.brdc[SysNum].eph[i].TOE);
					if(fabs(diff) <= BrdcTime[SysNum])
					{
						memcpy(&m_OrbClk->m_eph[Prn], &m_eph.brdc[SysNum].eph[i], sizeof(EPHN_t));
						if(fabs(diff) <= BrdcTime1[SysNum])
						{
							break;
						}
					}
				}
			}
		}
	}
	if(m_OrbClk->m_type == 1)
	{
		/*-------------------------------------------orbit--------------------------------------------- */
		u2 Numb = 0;
		if(((int)TNowTime.SoD) % m_OrbClk->m_SampleRate[0] <= 40)  //check whether need to update sp3
		{
			for(i = 0; i < m_eph.IGS.nsp3; i++)
			{
				diff = timediff(time, m_eph.IGS.sp3[i].time);
				if(fabs(diff) < (m_OrbClk->m_SampleRate[0] *10 +1))
				{
					memcpy(&m_OrbClk->m_sp3[Numb], &m_eph.IGS.sp3[i], sizeof(IGS_Orbit_t));
					Numb++;
				}	
				else if(diff < (m_OrbClk->m_SampleRate[0] *10 +1))
				{
					break;
				}
			}
			ASSERT(Numb < 22);
			m_OrbClk->m_nsp3 = Numb;
		}
		/*--------------------------------------------Clock---------------------------------------------*/
		if(((int)TNowTime.SoD) % m_OrbClk->m_SampleRate[1] <= 40)  //check whether need to update sp3
		{
			Numb = 0;			
			for(i = 0; i < m_eph.IGS.nclk; i++)
			{
				diff = timediff(time, m_eph.IGS.clk[i].time);
				if(fabs(diff) < (m_OrbClk->m_SampleRate[1] * 3 +1))
				{
					memcpy(&m_OrbClk->m_clk[Numb], &m_eph.IGS.clk[i], sizeof(IGS_Clk_t));
					Numb++;
				}
				else if(diff < (m_OrbClk->m_SampleRate[1] * 3 +1))
				{
					break;
				}
			}	
			ASSERT(Numb < 8);
			m_OrbClk->m_nclk = Numb;
		}	
	}
	return;
}

void Netcontrol::ProcessCompareObs(int MJDN)
{
	int i,j,StaNo;
	gtime_t NowGTime;
	TTime   NowTTime;
	bool Flag;
	GNSSDATA Obs[2] = {0};

	for(StaNo = 0; StaNo < m_Option.StaNum; StaNo++)
	{
		NowTTime.MJDN = MJDN;
		NowTTime.SoD  = 0;
		int year;
		double doy;
		t2doy(&NowTTime, &year, &doy);

		NowGTime = gltime2rtktime(NowTTime);

		char chOutFile[MAXLENGTH];
		sprintf(chOutFile, "%s/midoutput/%s%03d0.%02ddiff", m_Option.Directory.ProDir,m_Option.Station[StaNo].Name, (int)doy, year%100);
		FILE* ftOut = fopen(chOutFile, "w");

		fprintf(ftOut,"DifSat	");
		for(i = 0; i < MAXSAT; i++)
		{
			if(i < MAXPRNGPS)
			{
				fprintf(ftOut,"G%02d	G%02d	G%02d	G%02d	",i+1, i+1, i+1, i+1);
			}
			else if(i < MAXPRNGPS + MAXPRNBDS)
			{
				fprintf(ftOut,"C%02d	C%02d	C%02d	C%02d	",i-31, i-31, i-31, i-31);
			}
		}
		fprintf(ftOut,"\n");

		m_ObsLoad  = new RinexObsLoad[2];
        m_ObsLoad[0].InitRinexload(_bias, m_Option.Station[StaNo].RecType);
        m_ObsLoad[1].InitRinexload(_bias, m_Option.Station[StaNo].RecType);
		
		char chObsFile[MAXLENGTH];
		memset(&chObsFile, 0x00, sizeof(char) * MAXLENGTH);
		sprintf(chObsFile,"%s\\%04d\\%03d\\%s%03d0.%02do",
			m_Option.Directory.ObsDir, year, (int)doy,m_Option.Station[StaNo].Name, (int)doy, year%100);
		m_ObsLoad[0].Initial(chObsFile, m_Option.NavSys, NULL);
		m_ObsLoad[0].m_StaNo = StaNo;

		memset(&chObsFile, 0x00, sizeof(char) * MAXLENGTH);
		sprintf(chObsFile,"%s\\%04d\\%03d\\%s%03d0.%02do",
			m_Option.Directory.ObsRefDir, year,(int)doy ,m_Option.Station[StaNo].Name, (int)doy, year%100);
		m_ObsLoad[1].Initial(chObsFile, m_Option.NavSys, NULL);
		m_ObsLoad[1].m_StaNo = StaNo;

		int ObsType[6] = {0};
		while((int)NowTTime.SoD < 86400)
		{				
			Flag = m_ObsLoad[0].GetObsData(NowGTime, Obs[0]);
			if(Flag == false || Obs[0].NumSats < 1)
			{
				NowTTime.SoD  += m_Option.Interval;
				NowGTime.time += m_Option.Interval;
				continue;
			}	
			ASSERT(Obs[0].NumSats < (MAXCHN));

			Flag = m_ObsLoad[1].GetObsData(NowGTime, Obs[1]);
			if(Flag == false || Obs[1].NumSats < 1)
			{
				NowTTime.SoD  += m_Option.Interval;
				NowGTime.time += m_Option.Interval;
				continue;
			}	
			ASSERT(Obs[1].NumSats < (MAXCHN));
			int DifPrn = Obs[0].NumSats - Obs[1].NumSats;
			int RovNum = 0;
			double PObsDif[MAXSAT][2] = {0}, CObsDif[MAXSAT][2] = {0}, LObsDif[MAXSAT][2] = {0};
			for (i = 0; i < Obs[0].NumSats; i++)
			{
				for(j = RovNum; j < Obs[1].NumSats; j++)
				{
					if(Obs[0].ObsData[i].prn == Obs[1].ObsData[j].prn)
					{
						if(Obs[0].ObsData[i].P[0] > 0.0 && Obs[1].ObsData[j].P[0] > 0.0)
						{
							PObsDif[Obs[0].ObsData[i].prn-1][0] = Obs[0].ObsData[i].P[0] - Obs[1].ObsData[j].P[0];
							if(fabs(PObsDif[Obs[0].ObsData[i].prn-1][0]) > 0.1)
							{
								ObsType[0]++;
							}
						}
						if(Obs[0].ObsData[i].P[1] > 0.0 && Obs[1].ObsData[j].P[1] > 0.0)
						{
							PObsDif[Obs[0].ObsData[i].prn-1][1] = Obs[0].ObsData[i].P[1] - Obs[1].ObsData[j].P[1];
							if(fabs(PObsDif[Obs[0].ObsData[i].prn-1][1]) > 0.1)
							{
								ObsType[1]++;
							}
						}
						if(Obs[0].ObsData[i].C[0] > 0.0 && Obs[1].ObsData[j].C[0] > 0.0)
						{
							CObsDif[Obs[0].ObsData[i].prn-1][0] = Obs[0].ObsData[i].C[0] - Obs[1].ObsData[j].C[0];
							if(fabs(CObsDif[Obs[0].ObsData[i].prn-1][0]) > 0.1)
							{
								ObsType[2]++;
							}
						}
						if(Obs[0].ObsData[i].C[1] > 0.0 && Obs[1].ObsData[j].C[1] > 0.0)
						{
							CObsDif[Obs[0].ObsData[i].prn-1][1] = Obs[0].ObsData[i].C[1] - Obs[1].ObsData[j].C[1];
							if(fabs(CObsDif[Obs[0].ObsData[i].prn-1][1]) > 0.1)
							{
								ObsType[3]++;
							}
						}
						if(fabs(Obs[0].ObsData[i].L[0]) > 0.0 && fabs(Obs[1].ObsData[j].L[0]) > 0.0)
						{
							LObsDif[Obs[0].ObsData[i].prn-1][0] = Obs[0].ObsData[i].L[0] - Obs[1].ObsData[j].L[0];
						}
						if(fabs(Obs[0].ObsData[i].L[1]) > 0.0 && fabs(Obs[1].ObsData[j].L[1]) > 0.0)
						{
							LObsDif[Obs[0].ObsData[i].prn-1][1] = Obs[0].ObsData[i].L[1] - Obs[1].ObsData[j].L[1];
						}
						RovNum++;
						break;
					}
				}
			}
			fprintf(ftOut,"%3d	",DifPrn);
			for(i = 0; i < MAXSAT; i++)
			{
				fprintf(ftOut,"%4.2f	%4.2f	%4.2f	%4.2f	",LObsDif[i][0], LObsDif[i][1], CObsDif[i][0], CObsDif[i][1]);
			}
			fprintf(ftOut,"\n");
			NowTTime.SoD  += m_Option.Interval;
			NowGTime.time += m_Option.Interval;
		}
		printf("%s %d %d %d %d\n", m_Option.Station[StaNo].Name, ObsType[0], ObsType[1], ObsType[2], ObsType[3]);
		fclose(ftOut);
		delete[] m_ObsLoad;
	}
	return;
}
void InitialAmbInfo(char *FCBFile,char *DcbFile, AMB_Info &AInfo, int interval, int Type)
{
	int Number = 86400 / interval;
	
	for(int i = 0; i < MAXSAT; i++)
	{
		if(Type == 0)
		{
			AInfo.EWL[i] = (double *)calloc(sizeof(double), Number);
			AInfo.WL[i]  = (double *)calloc(sizeof(double), Number);
			AInfo.Var[i] = (double *)calloc(sizeof(double), Number);
		}
		memset(AInfo.EWL[i], 0, sizeof(double) * Number);
		memset(AInfo.WL[i],  0, sizeof(double) * Number);
		memset(AInfo.Var[i], 0, sizeof(double) * Number);
	}
	
	memset(AInfo.Numb, 0, sizeof(int) * MAXSAT);
	memset(AInfo.FixNumber,  0, sizeof(int) * 4 * 3);
	memset(AInfo.RcvFcb,  0, sizeof(double) * MAXSYS * 2);
	memset(AInfo.SatFloAmb, 0, sizeof(double) * MAXSAT * 2);
	memset(AInfo.DCBValue, 0, sizeof(double) * MAXSAT * 2);
	memset(AInfo.DCBNumb, 0, sizeof(int) * MAXSAT);
	
	
	if(Type == 0)
	{
		AInfo.NWLInfo[0] = 0;
		AInfo.NWLInfo[1] = 0;
		AInfo.WLInfo = NULL;
		memset(AInfo.FCB,    0, sizeof(int) * MAXSAT * 2);
		memset(AInfo.P_DCB,  0, sizeof(int) * MAXSAT * 2);
		FILE *fp = NULL;
		fp = fopen(FCBFile,"r");
		if(fp == NULL)
		{
			return;
		}
		else
		{
			char chLineTemp[MAXLENGTH];
			while(fgets(chLineTemp, MAXLENGTH, fp) != NULL)
			{
				if(chLineTemp[0] == '*')
				{
					int Type;
					if(chLineTemp[1] == 'E')
					{
						Type = 0;
					}
					else
					{
						Type = 1;
					}
					char Temp[10];
					int FCBNum,Prn;
					sscanf(chLineTemp, "%s %d",Temp, &FCBNum);
					for(int i = 0; i < FCBNum; i++)
					{
						if(fgets(chLineTemp, MAXLENGTH, fp) != NULL)
						{
							Prn = str2num(chLineTemp, 1, 2);
							if(chLineTemp[0] == 'C')
							{
								Prn += NSATGPS;
							}
							AInfo.FCB[Prn-1][Type] = str2num(chLineTemp, 3, 7);
						}
						else
						{
							break;
						}
					}

				}
			}
			fclose(fp);
			fp = NULL;
		}	
		fp = fopen(DcbFile, "r");
		if(fp == NULL)
		{
			return;
		}
		else
		{
			char chLineTemp[MAXLENGTH];
			while(fgets(chLineTemp, MAXLENGTH, fp) != NULL)
			{
				if(chLineTemp[0] == '*')
				{
					int Type;
					if(chLineTemp[2] == '1')
					{
						Type = 0;
					}
					else
					{
						Type = 1;
					}
					char Temp[10];
					int FCBNum,Prn;
					sscanf(chLineTemp, "%s %d",Temp, &FCBNum);
					for(int i = 0; i < FCBNum; i++)
					{
						if(fgets(chLineTemp, MAXLENGTH, fp) != NULL)
						{
							Prn = str2num(chLineTemp, 1, 2);
							if(chLineTemp[0] == 'C')
							{
								Prn += NSATGPS;
							}
							AInfo.P_DCB[Prn-1][Type] = str2num(chLineTemp, 3, 7);
						}
						else
						{
							break;
						}
					}

				}
			}
			fclose(fp);
			fp = NULL;
		}
	}
	else
	{
		AInfo.NWLInfo[0] = 0;
		if(AInfo.NWLInfo[1] > 0)
		{
			memset(AInfo.WLInfo, 0x00, sizeof(AInfo.WLInfo));
		}
	}
}
void OutPutIFB(FILE *ft, GNSSDATA Obs)
{
	int i;
	double f[4];
	double IFB[MAXSAT] = {0};
	f[0] = FREB1 * FREB1 / (FREB1 * FREB1 - FREB2 * FREB2);
	f[1] = 1 - f[0];
	f[2] = FREB1 * FREB1 / (FREB1 * FREB1 - FREB5 * FREB5);
	f[3] = 1 - f[2];
	f[0] = f[0] - f[2];
	for(i = 0; i < Obs.NumSats; i++)
	{
		if(fabs(Obs.ObsData[i].L[0]) > 0 && fabs(Obs.ObsData[i].L[1]) > 0 && fabs(Obs.ObsData[i].L[2]) > 0)
		{
			int prn_1 = Obs.ObsData[i].prn - 1;
			IFB[prn_1] = Obs.ObsData[i].L[0] * G_Lam[1][0] * f[0] + Obs.ObsData[i].L[1] * G_Lam[1][1] * f[1]
						- Obs.ObsData[i].L[2] * G_Lam[1][2] * f[3];
		}
	}
	for(i = MAXPRNGPS; i < MAXPRNGPS + MAXPRNBDS; i++)
	{
		fprintf(ft, "%10.5f	",IFB[i]);
	}
	fprintf(ft,"\n");
}
void OutputLog(FILE *ft_log, AMB_Info Ambinfo)
{
	int i;
	for(i = 0; i < Ambinfo.NWLInfo[0]; i++)
	{
		fprintf(ft_log, "C%02d %6d %6d %8d %8d\n",
			Ambinfo.WLInfo[i].Prn,Ambinfo.WLInfo[i].SSod, Ambinfo.WLInfo[i].ESod,
			Ambinfo.WLInfo[i].EWL, Ambinfo.WLInfo[i].WL);
	}
}
void UpdateWLInfo(AMB_Info AInfo, WL_Info *SatInfo, int Sod, int &ANum)
{
	for(int i = ANum; i < AInfo.NWLInfo[0]; i++)
	{
		if(AInfo.WLInfo[i].SSod == Sod)
		{
			memcpy(&SatInfo[AInfo.WLInfo[i].Prn-1], &AInfo.WLInfo[i], sizeof(WL_Info));
			ANum++;
		}
		else
		{
			return;
		}
	}
	return;
}
void ReSetPeusedorange(ObsData_t &Obs, int EWL, int WL, double FCB[2], double PDCB[2],double *deta)
{
	double P_EWL,P_WL,Iono;
	if(Obs.prn > MAXPRNGPS)
	{
		P_EWL = (Obs.L[2] - Obs.L[1] + EWL + FCB[0]) * G_EWL[1];
		P_WL  = (Obs.L[0] - Obs.L[1] +  WL + FCB[1]) * G_WL[1];
		Iono  = P_EWL - P_WL;
		P_EWL = P_EWL + Iono * -1.983053644 + PDCB[0];
		P_WL  = P_WL  + Iono *  1.271304966 + PDCB[1];
		deta[0] = Obs.P[0] - P_EWL;
		deta[1] = Obs.P[1] - P_WL;
		Obs.P[0] = P_EWL;
		Obs.P[1] = P_WL;
	}
	else
	{
		;
	}
	return;
}

void Netcontrol::OutPutWLRinex(FILE *ft, int MJDN, int StaNo, AMB_Info &AInfo)
{
	gtime_t NowGTime;
	TTime   NowTTime;
	bool    Flag;
	WL_Info SatInfo[MAXSAT] = {0};
	int     i,ANum = 0;

	double *SatXYZ,*Clk,*Azel,*EphRMS;

	SatXYZ = zeros(3, MAXCHN);
	Clk    = zeros(1, MAXCHN);
	Azel   = zeros(2, MAXCHN);
	EphRMS = zeros(1, MAXCHN);

	NowTTime.MJDN = MJDN;
	NowTTime.SoD = 0;
	NowGTime = gltime2rtktime(NowTTime);

	int year;
	double doy;
	t2doy(&NowTTime, &year, &doy);

	char ObsFile[200];
	sprintf(ObsFile,"%s\\%04d\\%03d\\%s%03d0.%02do",
		m_Option.Directory.ObsDir,year,(int)doy,m_Option.Station[StaNo].Name,(int)doy,year%100);
	m_ObsLoad[StaNo].Initial(ObsFile, m_Option.NavSys, ft);
	m_ObsLoad[StaNo].m_StaNo = StaNo;

	memset(&m_ObsData, 0x00, sizeof(GNSSDATA));

	while((int)NowTTime.SoD < 86400)
	{
		if(!m_Option.HighRate)
		{
			UpDataEph(NowGTime, NowTTime);
		}
		
		UpdateWLInfo(AInfo, SatInfo, (int)NowTTime.SoD, ANum);
		Flag = m_ObsLoad[StaNo].GetObsData(NowGTime, m_ObsData);

		NowTTime.SoD  += m_Option.Interval;
		NowGTime.time += m_Option.Interval;

		if(Flag == false || m_ObsData.NumSats < 1)
		{
			continue;
		}	

		ASSERT(m_ObsData.NumSats < (MAXCHN));

		if(m_OrbClk->GetOrb_Clk(m_ObsData, 0.0, SatXYZ, Clk, Azel, EphRMS,m_Option.Station[StaNo].XYZ, m_Option.Station[StaNo].Pos) < 1)
		{
			printf("%10d lack of enough empheris\n",(int)(NowTTime.SoD - m_Option.Interval));
			continue;
		}

		int ObsNum = 0;
		double Deta[MAXSAT][2] = {0};
		for(i = 0; i < m_ObsData.NumSats; i++)
		{
			int prn_1 = m_ObsData.ObsData[i].prn - 1;
			if(NowTTime.SoD >= SatInfo[prn_1].SSod && NowTTime.SoD < SatInfo[prn_1].ESod)
			{
				if(fabs(m_ObsData.ObsData[i].L[0]) > 0 && fabs(m_ObsData.ObsData[i].L[1]) > 0 
					&& fabs(m_ObsData.ObsData[i].L[2]) > 0)
				{
					ReSetPeusedorange(m_ObsData.ObsData[i], SatInfo[prn_1].EWL, SatInfo[prn_1].WL,
						AInfo.FCB[prn_1], AInfo.P_DCB[prn_1], &Deta[prn_1][0]);
		
					AInfo.DCBValue[prn_1][0] = (AInfo.DCBValue[prn_1][0] * AInfo.DCBNumb[prn_1] + Deta[prn_1][0])
						/(AInfo.DCBNumb[prn_1]  + 1);
					AInfo.DCBValue[prn_1][1] = (AInfo.DCBValue[prn_1][1] * AInfo.DCBNumb[prn_1]  + Deta[prn_1][1])
						/(AInfo.DCBNumb[prn_1]  + 1);
					AInfo.DCBNumb[prn_1]++;

					memcpy(&m_ObsData.ObsData[ObsNum], &m_ObsData.ObsData[i], sizeof(ObsData_t));
					ObsNum++;
				}
			}
		}
		m_ObsData.NumSats = ObsNum;
// 		for(i = MAXPRNGPS; i < MAXPRNGPS + MAXPRNBDS; i++)
// 		{
// 			fprintf(ft,"%7.3f	%7.3f	",Deta[i][0],Deta[i][1]);
// 		}
// 		fprintf(ft,"\n");

		if(ObsNum > 0)
		{
			double ep[6];
			time2epoch(m_ObsData.ObsTime, ep);
			fprintf(ft,"> %04d %02d %02d %02d %02d %10.7f  0 %02d\n",
				(int)ep[0], (int)ep[1], (int)ep[2], (int)ep[3], (int)ep[4], ep[5],ObsNum);
			for(i = 0; i < ObsNum; i++)
			{
				fprintf(ft,"C%02d%14.3f  %14.3f  %14.3f  %14.3f  %14.3f  %14.3f\n",m_ObsData.ObsData[i].prn-32,
					m_ObsData.ObsData[i].P[0], m_ObsData.ObsData[i].P[1], 0.0,
					m_ObsData.ObsData[i].L[0], m_ObsData.ObsData[i].L[1], 0.0);
			}
		}
	}
	free(SatXYZ); free(Clk);  free(Azel);   free(EphRMS);
}
void Netcontrol::LoadNav(int MJDN)
{
	TTime   NowTTime;
	NowTTime.MJDN = MJDN;
	NowTTime.SoD  = 0;
	int year;
	double doy;
	t2doy(&NowTTime, &year, &doy);

	m_OrbClk->m_type = 0;
	m_NavLoad.ReadNavData(m_eph);
	if(m_OrbLoad.ReadOrbitData(m_eph.IGS))
	{
		m_OrbClk->m_type = 1;
		m_OrbClk->m_SampleRate[0] = timediff(m_eph.IGS.sp3[1].time, m_eph.IGS.sp3[0].time);
	}
	m_ClkLoad.ReadClkData(m_eph.IGS);	
	m_OrbClk->m_SampleRate[1] = m_eph.IGS.SamplRate[1];
	return;
}

void Netcontrol::Process(int MJDN)
{
	int i,j,StaNo;
	gtime_t NowGTime;
	TTime   NowTTime;
	bool Flag;
	double *SatXYZ,*Clk,*Azel,*EphRMS;
	
	SatXYZ = zeros(3, MAXCHN);
	Clk    = zeros(1, MAXCHN);
	Azel   = zeros(2, MAXCHN);
	EphRMS = zeros(1, MAXCHN);

	LoadNav(MJDN);

	m_OrbClk->Initialize(MJDN);

	NowTTime.MJDN = MJDN;
	NowTTime.SoD  = 0;

	for(StaNo = 0; StaNo < m_Option.StaNum; StaNo++)
	{
		m_GNSSInfo[StaNo].N_Rove = StaNo;
		m_GNSSInfo[StaNo].lastepoch = false;
	}

	for(StaNo = 0; StaNo < m_Option.StaNum; StaNo++)
	{
		printf("%s\n",m_Option.Station[StaNo].Name);
		NowTTime.SoD = 0;
		NowGTime = gltime2rtktime(NowTTime);
		while((int)NowTTime.SoD < 86400)
		{	
			UpDataEph(NowGTime, NowTTime);	
			m_GNSSInfo[StaNo].LogSoD = (int)NowTTime.SoD;

			Flag = m_ObsLoad[StaNo].GetObsData(NowGTime, m_ObsData);
			
			NowTTime.SoD  += m_Option.Interval;
			NowGTime.time += m_Option.Interval;

			if(Flag == false || m_ObsData.NumSats < 1)
			{
				continue;
			}	
			ASSERT(m_ObsData.NumSats < (MAXCHN));

			m_GNSSInfo[StaNo].UseFulNum++;
			if(m_OrbClk->GetOrb_Clk(m_ObsData, 0.0, SatXYZ, Clk, Azel, EphRMS,m_Option.Station[StaNo].XYZ, m_Option.Station[StaNo].Pos) < 1)
			{
				continue;
			}
			AddGFIF(m_ObsData, Azel, m_GNSSInfo[StaNo] ,m_Option.Station[StaNo].ft_mw);
			EstFCB(m_ObsData, Azel, m_GNSSInfo[StaNo], m_Option.Station[StaNo].ft_mw);
		}	//while
// 		AddLastTime(m_GNSSInfo[StaNo]);
	}
	free(SatXYZ); free(Clk);  free(Azel);   free(EphRMS);
	return;
}

void Netcontrol::OutPutUPDStatis()
{
	FILE* fp = fopen("WLStatis.txt", "w");
	//======================================OutPut
	for(int StaNo = 0; StaNo < m_Option.StaNum; StaNo++)
	{
		AddLastTime(m_GNSSInfo[StaNo]);
	}
	int nSat[2] = {0, NSATGPS + NSATBDS};
	if(m_Option.NavSys == SYS_GPS)
	{
		nSat[0] = 0;
		nSat[1] = NSATGPS;
	}
	else if(m_Option.NavSys == SYS_BDS)
	{
		nSat[0] = NSATGPS;
		nSat[1] = NSATGPS + NSATBDS;
	}
	ComputeAverage(&_UPDInfo, m_GNSSInfo, m_Option.StaNum, 0, fp, nSat);
// 	ComputeAverage(&_UPDInfo, m_GNSSInfo, m_Option.StaNum, 1, fp, nSat);
	fclose(fp);
	//======================================OutPut
	int RefSite = -1, MAXObs = 0;
	int MinArc = 100;
	for(int StaNo = 0; StaNo < m_Option.StaNum; StaNo++)
	{
		int nTemp = 0;
		for(int prn_1 = nSat[0]; prn_1 < nSat[1]; prn_1++)
		{
			if(m_GNSSInfo[StaNo].iPIFGF[prn_1] > MinArc)
			{
				nTemp++;
			}
		}
		if(nTemp > MAXObs)
		{
			RefSite = StaNo;
			MAXObs = nTemp;
		}
	}
   FILE* ftPIFGF;
   char chFile[255];
   sprintf(chFile,"PIFGF.txt");
   ftPIFGF = fopen(chFile, "w");
   for(int StaNo = 0; StaNo < m_Option.StaNum; StaNo++)
   {
	   fprintf(ftPIFGF,"%s, ", m_GNSSInfo[StaNo].Station[StaNo].Name);
	   double dTmp = 0.0;
	   int nTmp = 0;
	   for(int prn_1 = nSat[0]; prn_1 < nSat[1]; prn_1++)
	   {
		   if(m_GNSSInfo[StaNo].iPIFGF[prn_1] > MinArc && m_GNSSInfo[RefSite].iPIFGF[prn_1] > MinArc)
		   {
			   dTmp += m_GNSSInfo[StaNo].dPIFGF[prn_1] - m_GNSSInfo[RefSite].dPIFGF[prn_1];
			   nTmp++;
		   }
	   }
	   if(nTmp > 0)
	   {
		   dTmp = dTmp / nTmp;
	   }
	   for(int prn_1 = nSat[0]; prn_1 < nSat[1]; prn_1++)
	   {
		   if(m_GNSSInfo[StaNo].iPIFGF[prn_1] > MinArc && m_GNSSInfo[RefSite].iPIFGF[prn_1] > MinArc)
		   {
			   fprintf(ftPIFGF,"%6.2f, ", m_GNSSInfo[StaNo].dPIFGF[prn_1] - m_GNSSInfo[RefSite].dPIFGF[prn_1] - dTmp);
		   }
		   else
		   {
			   fprintf(ftPIFGF,"   -  , ");
		   }
	   }
	   fprintf(ftPIFGF,"\n");
   }
   fclose(ftPIFGF);

	FILE *ftWL = NULL, *ftEWL;

	ftWL  = fopen("WLStatistic.txt",  "w");
	ftEWL = fopen("EWLStatistic.txt", "w");
	for(int i = 0; i < _UPDInfo._nDay; i++)
	{
		fprintf(ftWL, "Day_%05d ", m_Option.StartTime.MJDN + i);
		fprintf(ftEWL,"Day_%05d ", m_Option.StartTime.MJDN + i);
		for(int prn_1 = 0; prn_1 < MAXSAT; prn_1++)
		{
			fprintf(ftWL,  "%7.4f ", _UPDInfo._dWLSeri[i * MAXSAT + prn_1]);
			fprintf(ftEWL, "%7.4f ", _UPDInfo._dEWLSeri[i * MAXSAT + prn_1]);
		}
		fprintf(ftWL, "\n");
		fprintf(ftEWL, "\n");
	}
	fclose(ftWL);
	fclose(ftEWL);
}