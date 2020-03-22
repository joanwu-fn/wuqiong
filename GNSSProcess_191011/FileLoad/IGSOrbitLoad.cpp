#include "../Common/FileLoad.h"
#include "../Common/BaseFunction.h"
#include "../Common/debug.h"
IGSOrbitLoad::IGSOrbitLoad()
{
	m_File = NULL;
}
IGSOrbitLoad::~IGSOrbitLoad()
{
	;
}
bool IGSOrbitLoad::Initial(char chFile[])
{
	m_File = fopen(chFile, "r");
	if(m_File == NULL)
		return false;
	return true;
}

bool IGSOrbitLoad::ReadOrbitHeader()
{
	return true;
}
bool IGSOrbitLoad::ReadOrbitData(GNSSIGS_t &IGS)
{
	char buff[MAXLENGTH],StrSys;
	double ep[6]={0},pos[3];
	u2 nPrn=0, temp,prn,NAVSYS = SYS_GPS;
	int sys;
	gtime_t NowTime= {0};
	IGS_Orbit_t *IGSOrbit;
	bool InitialFlag = false;
	
	if(m_File == NULL)
		return false;

	while (fgets(buff,MAXLENGTH,m_File) != NULL)
	{
		if(buff[0] != '*' && buff[0] != 'P')
			continue;
		if(IGS.maxsp3 <= IGS.nsp3)
		{
			IGS.maxsp3 += 100;
			if ((IGSOrbit=(IGS_Orbit_t *)realloc(IGS.sp3, sizeof(IGS_Orbit_t)*IGS.maxsp3)) == NULL)
			{
				ASSERT(3 < 0);
			}
			IGS.sp3 = IGSOrbit;
		}
		if(buff[0] == '*')
		{
			InitialFlag = true;
			ep[0] = str2num(buff, 3, 4);         ep[1] = str2num(buff, 8, 2);
			ep[2] = str2num(buff, 11, 2);        ep[3] = str2num(buff, 14, 2);
			ep[4] = str2num(buff, 17, 2);        ep[5] = str2num(buff, 20, 11);
			NowTime = epoch2time(ep);
			memset(&IGS.sp3[IGS.nsp3], 0x00, sizeof(IGS_Orbit_t));
			memcpy(&IGS.sp3[IGS.nsp3].time, &NowTime, sizeof(gtime_t));			
			IGS.nsp3++;
		}
		else if(buff[0] == 'P' && InitialFlag == true)
		{
			StrSys =   buff[1];
			if((sys = JudgeSys(StrSys, SYS_GPS| SYS_BDS | SYS_GAL, &NAVSYS)) == -1)
				continue;
			temp = (int)str2num(buff, 2, 2);
			prn = SatNo(temp, NAVSYS);
			pos[0] = str2num(buff,  5, 13)*1000.0;
			pos[1] = str2num(buff, 19, 13)*1000.0;
			pos[2] = str2num(buff, 33, 13)*1000.0;
			
			IGS.sp3[IGS.nsp3-1].prn[prn-1] = prn;
			IGS.sp3[IGS.nsp3-1].Orbit[prn-1][0] = pos[0];
			IGS.sp3[IGS.nsp3-1].Orbit[prn-1][1] = pos[1];
			IGS.sp3[IGS.nsp3-1].Orbit[prn-1][2] = pos[2];
		}
	}  // while
	fclose(m_File);
	m_File = NULL;
	return true;
}