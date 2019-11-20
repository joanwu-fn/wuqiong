#include "../Common/BaseFunction.h"

/*****************************************************************************
* Name        : trim
* Description : Trim a string (remove preceding and succeding spaces)
* Parameters  :
* Name                           |Da|Unit|Description
* char  *line                     IO N/A  String to trim
* Returned value (char*)          O  N/A  Pointer to *line
*****************************************************************************/
char *trim(char *line) {
	int i;
	int pos=0;
	int lastnotspace=0;
	int state=0;
	// 0 -> Trimming first spaces
	// 1 -> Waitting to end string

	for (i=0;i<strlen(line);i++){
		if (state==0) {
			if (line[i]!=' ') {
				state = 1;
				line[pos] = line[i];
				pos++;
				lastnotspace = pos;
			}
		} else if (state==1) {
			line[pos] = line[i];
			pos++;
			if (line[i]!=' ' && line[i]!='\n') {
				lastnotspace = pos;
			}
		}
	}
	line[lastnotspace]='\0';
	return &line[0];
}
/*****************************************************************************
* Name        : getstr
* Description : Get a substring from a string and trim it
* Parameters  :
* Name                           |Da|Unit|Description
* char  *line                     I  N/A  Input string
* char  *out                      O  N/A  Substring
* int  ini                        I  N/A  Initial character position
* int  length                     I  N/A  Length of the substring
* Returned value (char*)          O  N/A  Pointer to *out
*****************************************************************************/
char *getstr(char *out,char *line,int ini, int length) {
	int end;

	end = strlen(line);
	if (ini>=end) {
		out[0]='\0';
		return &out[0];
	}
	strncpy(out,&line[ini],length);
	out[length]='\0';
	trim(out);
	return &out[0];
}

int JudgeSys(char strSys, u2 NAVSYS, u2 *SYS)
{
	if(strSys == 'G' || strSys == ' '|| strSys == 'M')  //'M'针对几何法定轨的地面站
	{
		if( SYS != NULL)
		{
			*SYS = SYS_GPS;
		}
		if(NAVSYS == -1)
			return 0;
		else
		{
			if(NAVSYS & 1)
				return 0;
			else 
				return -1;
		}
	}
	else if(strSys == 'C')
	{
		if( SYS != NULL)
		{
			*SYS = SYS_BDS;
		}
		if(NAVSYS == -1)
			return 1;
		else
		{
			if(NAVSYS & 2)
				return 1;
			else 
				return -1;
		}
	}
	else if(strSys == 'E')
	{
		if( SYS != NULL)
		{
			*SYS = SYS_GAL;
		}
		if(NAVSYS == -1)
			return 2;
		else
		{
			if(NAVSYS & 4)
				return 2;
			else 
				return -1;
		}
	}
	else if(strSys == 'R')
	{
		if( SYS != NULL)
		{
			*SYS = SYS_GLO;
		}
		if(NAVSYS == -1)
			return 3;
		else
		{
			if(NAVSYS & 8)
				return 3;
			else 
				return -1;
		}
	}
	else 
	{
		if( SYS != NULL)
		{
			*SYS = SYS_NONE;
		}
		return -1;
	}
}

u2 SatNo(u2 Sat, u2 NAVSYS)
{
	if(Sat <= 0)
		return 0;
	switch(NAVSYS)
	{
	   case SYS_GPS:
		   if (Sat < MINPRNGPS || MAXPRNGPS < Sat) return 0;
		   return Sat-MINPRNGPS+1;
	   case SYS_BDS:
		   if (Sat < MINPRNBDS||MAXPRNBDS < Sat) return 0;
		   return NSATGPS+Sat-MINPRNBDS+1;
	   case SYS_GAL:
		   if(Sat < MINPRNGAL||MAXPRNGAL < Sat) return 0;
		   return NSATGPS+NSATBDS+Sat-MINPRNGAL+1;
	   case SYS_GLO:
		   if(Sat < MINPRNGLO||MAXPRNGLO < Sat) return 0;
		   return NSATGPS+NSATBDS+NSATGAL+Sat-MINPRNGLO;
	   default:
		   return 0;
	}
	return 0;
}

u2 Prn2Sys(u2 Prn, u2 &SYS, u2 &SysNum, char &StrSys)
{
	if(Prn > MAXSAT || Prn <= 0)     //NONE
	{
		SYS = SYS_NONE;
		StrSys = 'N';
		SysNum = 0;
		return 0;
	}
	else if(Prn <= NSATGPS)                               //GPS
	{
		SysNum = 0;
		SYS = SYS_GPS;
		StrSys = 'G';
		return Prn;
	}
	else if(Prn <= (NSATGPS + NSATBDS))                   //BDS
	{
		SYS = SYS_BDS;
		StrSys = 'C';
		SysNum = 1;
		return Prn - NSATGPS;
	}
	else if(Prn <= (NSATGPS + NSATBDS + NSATGAL))         //GALILEO
	{
		SYS = SYS_GAL;
		StrSys = 'E';
		SysNum = 2;
		return Prn - NSATGPS - NSATBDS;
	}
	else if(Prn <= (NSATGPS + NSATBDS +NSATGAL + NSATGLO))   //GLONASS
	{
		SYS = SYS_GLO;
		StrSys = 'R';
		SysNum = 3;
		return Prn - NSATGPS - NSATBDS -NSATGAL;
	}
	else
	{
		SYS = SYS_NONE;
		StrSys = 'N';
		SysNum = 0;
		return 0;
	}		
}
/* satellite id to satellite number --------------------------------------------
* convert satellite id to satellite number
* args   : char   *id       I   satellite id (nn,Gnn,Rnn,Enn,Jnn,Cnn or Snn)
* return : satellite number (0: error)
* notes  : 120-138 and 193-195 are also recognized as sbas and qzss
*-----------------------------------------------------------------------------*/
int Satid2No(const char *id)
{
	int sys,prn;
	char code;

	if (sscanf(id,"%d",&prn)==1) {
		if      (MINPRNGPS<=prn&&prn<=MAXPRNGPS) sys=SYS_GPS;
		else return 0;
		return SatNo(prn, sys);
	}
	if (sscanf(id,"%c%d",&code,&prn)<2) return 0;

	switch (code) {
		case 'G': sys=SYS_GPS;  break;
		case 'C': sys=SYS_BDS;  break;
		case 'R': sys=SYS_GLO;  break;
		case 'E': sys=SYS_GAL;  break;	
		default: return 0;
	}
	return SatNo(prn, sys);
}
void OutPutObs(FILE *fp_Obs, GNSSDATA ObsData, u2 NAVSYS)
{
	if(fp_Obs == NULL)
		return;
	u2 i,prn,SYS,SysNum;
	double ep[6];
	char StrSys='G';
	time2epoch(ObsData.ObsTime,ep);
	fprintf(fp_Obs,"> %04d %02d %02d %02d %02d %10.7f  0 %2d\n",(int)ep[0],(int)ep[1],(int)ep[2],(int)ep[3],(int)ep[4],ep[5],ObsData.NumSats);
	for(i =0; i < ObsData.NumSats; i++)
	{
		prn = Prn2Sys(ObsData.ObsData[i].prn, SYS, SysNum, StrSys);
		fprintf(fp_Obs,"%c%02d%14.3f  %14.3f  %14.3f  %14.3f  %14.3f  %14.3f  \n",
			StrSys,prn,ObsData.ObsData[i].P[0],ObsData.ObsData[i].P[1],ObsData.ObsData[i].P[2],
			ObsData.ObsData[i].L[0],ObsData.ObsData[i].L[1],ObsData.ObsData[i].L[2]);
	}
}
void AddWLInfo(AMB_Info &AmbInfo, int Prn, int SSod, int ESod)
{
	if(AmbInfo.Flag[Prn-1][0] && AmbInfo.Flag[Prn-1][1])
	{
		WL_Info *WLTemp;
		if(AmbInfo.NWLInfo[0] >= AmbInfo.NWLInfo[1])
		{
			AmbInfo.NWLInfo[1] += 200;
			if ((WLTemp=(WL_Info *)realloc(AmbInfo.WLInfo, sizeof(WL_Info) * AmbInfo.NWLInfo[1])) == NULL)
			{
				exit(0);
			}
			AmbInfo.WLInfo = WLTemp;
		}
		AmbInfo.WLInfo[AmbInfo.NWLInfo[0]].Prn  = Prn;
		AmbInfo.WLInfo[AmbInfo.NWLInfo[0]].ESod = ESod;
		AmbInfo.WLInfo[AmbInfo.NWLInfo[0]].SSod = SSod;
		AmbInfo.WLInfo[AmbInfo.NWLInfo[0]].EWL  = AmbInfo.Amb[Prn-1][0];
		AmbInfo.WLInfo[AmbInfo.NWLInfo[0]].WL   = AmbInfo.Amb[Prn-1][1];
		AmbInfo.NWLInfo[0]++;
	}	
}