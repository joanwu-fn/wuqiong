#include "../Common/FileLoad.h"
#include "../Common/BaseFunction.h"
#include "../Common/debug.h"
IGSClkLoad::~IGSClkLoad()
{
	;
}
IGSClkLoad::IGSClkLoad()
{
	m_File == NULL;
}
bool IGSClkLoad::Initial(char chFile[])
{
	m_File = fopen(chFile, "r");
	if(m_File == NULL)
		return false;
	return ReadClkHeader();
}
bool IGSClkLoad::ReadClkHeader()
{
	char buff[MAXLENGTH],*label=buff+60;

	while (fgets(buff,MAXLENGTH,m_File) != NULL) {
		if (strstr(label,"END OF HEADER") != NULL) 
		{
			return true;
		}
		else if(strstr(label,"INTERVAL") != NULL)
		{
			;
		}
	}
	return false;
}
bool IGSClkLoad::ReadClkData(GNSSIGS_t &IGS)
{
	char buff[MAXLENGTH],StrSys;
	double ep[6]={0},clk,detaT;
	u2  temp,prn,NAVSYS = SYS_GPS;
	int sys;
	gtime_t NowTime= {0},LastTime ={0};
	IGS_Clk_t *IGSCLK;

	if(m_File == NULL)
		return false;
	while (fgets(buff,MAXLENGTH,m_File) != NULL)
	{
		if(buff[0] != 'A' || buff[1] != 'S')
			continue;
		StrSys =   buff[3];
		if((sys = JudgeSys(StrSys, SYS_GPS|SYS_BDS|SYS_GAL, &NAVSYS)) == -1)
			continue;
		xstrReplace(buff,'D','E',4);
		temp    =   (int)str2num(buff, 4, 2);
		prn = SatNo(temp, NAVSYS);
		ep[0]  =    str2num(buff,8,4);         ep[1]  =    str2num(buff,13,2);
		ep[2]  =    str2num(buff,16,2);        ep[3]  =    str2num(buff,19,2);
		ep[4]  =    str2num(buff,22,2);        ep[5]  =    str2num(buff,25,9);
		NowTime = epoch2time(ep);
		clk = str2num(buff, 40, 19);
		detaT = NowTime.time - LastTime.time + NowTime.sec -NowTime.sec;
		if(fabs(detaT) > 1e-5)
		{
			if(IGS.maxclk <= IGS.nclk)
			{
				IGS.maxclk += 300;
				if ((IGSCLK=(IGS_Clk_t *)realloc(IGS.clk, sizeof(IGS_Clk_t)*IGS.maxclk)) == NULL)
				{
					ASSERT(2 < 0);
				}
				IGS.clk = IGSCLK;
			}
			if(IGS.nclk == 1)
			{
				IGS.SamplRate[1] = (u2)timediff(NowTime, LastTime);
			}
			memset(&IGS.clk[IGS.nclk], 0x00, sizeof(IGS_Clk_t));
			memcpy(&LastTime, &NowTime, sizeof(gtime_t));
			memcpy(&IGS.clk[IGS.nclk].time, &LastTime, sizeof(gtime_t));		
			IGS.nclk++;		
		}
		IGS.clk[IGS.nclk - 1].prn[prn-1] = prn;
		IGS.clk[IGS.nclk - 1].SatClk[prn-1] = clk;
	}
	fclose(m_File);
	m_File = NULL;
	return true;
}