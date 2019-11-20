#include "../Common/FileLoad.h"
#include "../Common/BaseFunction.h"
#include "../Common/debug.h"

const char codes[]="CLDS";    /* obs type codes */

char *obscodes[]=
{       /* observation code strings */
	""  ,"1C","1P","1W","1Y", "1Q","1N","1S","1L","1E",
	"1A","1B","1X","1Z","2C", "2D","2S","2L","2X","2P",
	"2W","2Y","5C","5W","5I", "5Q","5X","7I","7Q","7X",
	"6A","6B","6C","6X","6Z", "6S","6L","8L","8Q","8X",
	"2I","2Q","6I","6Q","3P", "3L","3C","1I",""
};
unsigned char obsfreqs[]={ /* 1:L1,2:L2,3:L5,4:L6,5:L7,6:L8 */
	0, 1, 1, 1, 1,      1, 1, 1, 1, 1, 
	1, 1, 1, 1, 2,      2, 2, 2, 2, 2,
	2, 2, 3, 3, 3,      3, 3, 2, 5, 5, 
	4, 4, 4, 4, 4,      4, 4, 6, 6, 6,
	1, 2, 3, 4, 3,      3, 3, 1, 0
};
char codepris[4][16]={  /* code priority table {L1,L2,L5,L6,L7,L8} */
	{"CPWXYMNSLQ"}, /* GPS */
	{"CI"},                /* BDS */	
	{"QCABXZ"},   /* GAL */
	{"PC"}                   /* GLO */
};
unsigned char RinexObsLoad::obs2code(const char *obs, int sys, u2 &freq)
{
	freq = 0;
	if(0 == sys)  //GPS
	{
		if(   strcmp("1C",obs) == 0 || strcmp("1P",obs) == 0 || strcmp("1W",obs) == 0)
		{
			freq = 1;
		}
		else if(strcmp("2C",obs) == 0 || strcmp("2P",obs) == 0 || strcmp("2W",obs) == 0 
			 || strcmp("2X",obs) == 0)
		{
			freq = 2;
		}
		else if(strcmp("5C",obs) == 0 || strcmp("5W",obs) == 0 || strcmp("5Q",obs) == 0
			|| strcmp("5I",obs) == 0  || strcmp("5X",obs) == 0 || strcmp("5S",obs) == 0)
		{
			freq = 3;
		}
		else
		{
			freq = 0;
		}
	} 
	else if(1 == sys) //BDS
	{
// 		if(    strcmp("1C",obs) == 0 || strcmp("1I",obs) == 0 || strcmp("1Q",obs) == 0
// 			|| strcmp("2I",obs) == 0 || strcmp("1P",obs) == 0 || strcmp("2Q",obs) == 0
// 			|| strcmp("1W",obs) == 0)
		if(strcmp("2I",obs) == 0 || strcmp("1I",obs) == 0)
		{
			freq = 1;
		}
// 		else if(strcmp("2C",obs) == 0 || strcmp("2P",obs) == 0 || strcmp("7I",obs) == 0
// 			 || strcmp("7Q",obs) == 0 || strcmp("2W",obs) == 0)
		else if(strcmp("6I",obs) == 0)
		{
			freq = 3;
		}
// 		else if(strcmp("3C",obs) == 0 || strcmp("3P",obs) == 0 || strcmp("6I",obs) == 0
// 			|| strcmp("6Q",obs) == 0)
		else if(strcmp("5I",obs) == 0 || strcmp("5Q",obs) == 0 || strcmp("7I",obs) == 0)
		{
			freq = 2;
		}
		else
		{
			freq = 0;
		}
	}
	else if(2 == sys)  //GAL
	{
		if(   strcmp("1C",obs) == 0 || strcmp("1P",obs) == 0 || strcmp("1W",obs) == 0 || strcmp("1Q",obs) == 0
			||strcmp("1I",obs) == 0 || strcmp("1X",obs) == 0)
		{
			freq = 1;
		}
		else if(strcmp("2C",obs) == 0 || strcmp("2P",obs) == 0 || strcmp("2W",obs) == 0 || strcmp("5Q",obs) == 0
			 || strcmp("5I",obs) == 0 || strcmp("5X",obs) == 0 || strcmp("5S",obs) == 0)
		{
			freq = 2;
		}
		else if(strcmp("3C",obs) == 0 || strcmp("3P",obs) == 0 || strcmp("7Q",obs) == 0 || strcmp("7X",obs) == 0)
		{
			freq = 3;
		}
		else
		{
			freq = 0;
		}
	} 
	return 0;
// 	int i;
// 	freq=0;
// 	for (i=1;*obscodes[i];i++) 
// 	{
// 		if (strcmp(obscodes[i],obs) != 0)
// 		{
// 			continue;
// 		}
// 		 freq=(u2)obsfreqs[i];
// 		return (unsigned char)i;
// 	}
// 	return 0;
}
RinexObsLoad::RinexObsLoad()
{
	memset(&_bias, 0x00, sizeof(_bias));
	_RecType = 0;
}
void RinexObsLoad::InitRinexload(GNSSCodeBias bias, int RecType)
{
	m_Version      =  0;
	m_File         =  NULL;
	m_NavSys       =  0;
	_RecType       =  RecType;
	memcpy(&_bias,  &bias, sizeof(GNSSCodeBias));
	memset(&m_NowObsData, 0x00, sizeof(GNSSDATA));
	return;
}
RinexObsLoad::~RinexObsLoad()
{
	;
}

//NAME: Initial()
//Initial the file name and system flag , reading the file header
//true :successful    false: failed
bool RinexObsLoad::Initial(char chFile[], u2 NavSys, FILE *ft_out)
{
	if(m_File != NULL)
	{
		fclose(m_File);
		m_File = NULL;
	}
	m_File   = fopen(chFile, "r");
	if(m_File == NULL)
	{
		return false;
	}

	m_NavSys = NavSys;
	memset(&m_NowObsData, 0x00, sizeof(GNSSDATA));
	return ReadObsHeader(ft_out);
}
bool RinexObsLoad::GetObsType_2(char buff[])
{
	u2 number,j,k;
	char aux[4];
	const char *p;

	number = (int)str2num(buff,0,6);
	m_ObsType[0].ntype = number;
	for (j=0,k=10; j<number; j++,k+=6) 
	{
		if(k > 58)
		{
			if (!fgets(buff,MAXLENGTH,m_File)) 
			{
				return false;
			}
			k = 10;
		}
		getstr(aux, buff, k, 3);
		if(aux[0] == 'C' || aux[0] == 'c')
		{
			aux[2] = 'W';
		}
		else		
		{
			aux[2] = 'C';
		}

		obs2code(&aux[1], 0, m_ObsType[0].Frequncy[j]);
		m_ObsType[0].Type[j]     = (p=strchr(codes,aux[0]))?(int)(p-codes):0;
		m_ObsType[0].priority[j] = 16 - ((p=strchr(codepris[0],aux[2]))?(int)(p-codepris[0]):16);

		if((aux[0] == 'C' || aux[0] == 'c') && (aux[1] == '1')) //P1C1
		{
			m_ObsType[0].bP1C1[j] = true;
		}
		else
		{
			m_ObsType[0].bP1C1[j] = false;
		}

		if((aux[0] == 'C' || aux[0] == 'c') && (aux[1] == '2'))
		{
			m_ObsType[0].bP2C2[j] = true;
		}
		else
		{
			m_ObsType[0].bP2C2[j] = false;
		}
		
		if(m_ObsType[0].Type[j] < 2 && m_ObsType[0].Frequncy[j] > 0) //Puesdo-range  and  Carrier phase
		{
			m_ObsType[0].valid[j] = 1;
		}
	}
	memcpy(&m_ObsType[1], &m_ObsType[0], sizeof(OBSTYPE));
	memcpy(&m_ObsType[2], &m_ObsType[0], sizeof(OBSTYPE));
	memcpy(&m_ObsType[3], &m_ObsType[0], sizeof(OBSTYPE));
	return true;
}
bool RinexObsLoad::GetObsType_3(char buff[])
{
	u2 number;
	u2 j,k;
	int sys;

	sys =JudgeSys(buff[0], m_NavSys, NULL);
	if(sys == -1)
	{
		number = (int)str2num(buff,3,3);
		if(number <= 13)
		{
			return false;
		}
		fgets(buff,MAXLENGTH,m_File);
		return false;
	}
	number = (int)str2num(buff,3,3);
	m_ObsType[sys].ntype = number;

	ASSERT(number <= MAXENTRY_NUMBER);
	char aux[4];
	const char *p;
	bool bC2I = false, bC1I = false;
	int iPos = 0;
	for (j=0,k=7; j<number; j++,k+=4) 
	{
		if (k>58)
		{
			if (!fgets(buff,MAXLENGTH,m_File)) 
			{
				return false;
			}
			k=7;
		}
		getstr(aux, buff, k, 3);
		if(aux[2] == '\0')
		{
			aux[2] = 'C';
		}

		obs2code(&aux[1], sys, m_ObsType[sys].Frequncy[j]);

		if(aux[1] == '1' && aux[2] == 'I')
		{
			iPos = j;
			bC1I = true;
		}
		if(aux[1] == '2' && aux[2] == 'I')
		{
			bC2I = true;
		}

		m_ObsType[sys].Type[j]     = (p=strchr(codes,aux[0]))?(int)(p-codes):0;
		m_ObsType[sys].priority[j] = 16 - ((p=strchr(codepris[sys],aux[2]))?(int)(p-codepris[sys]):16);
		if(aux[0] == 'C' && (aux[2] == 'C' || aux[2] == 'c') && sys == 0 && aux[1] == '1')   //C1C
		{
			m_ObsType[sys].bP1C1[j] = true;
		}
		else
		{
			m_ObsType[sys].bP1C1[j] = false;
		}
		if(aux[0] == 'C' && (aux[2] == 'X' || aux[2] == 'C') && sys == 0 && aux[1] == '2')  //C2X C2C
		{
			m_ObsType[sys].bP2C2[j] = true;
		}
		else
		{
			m_ObsType[sys].bP2C2[j] = false;
		}

		if(m_ObsType[sys].Type[j] < 2 && m_ObsType[sys].Frequncy[j] > 0) //Puesdo-range  and  Carrier phase
		{
			m_ObsType[sys].valid[j] = 1;
		}
	}
	if(buff[0] == 'C' && bC1I && bC2I)
	{
		m_ObsType[sys].Frequncy[iPos] = -1;
		m_ObsType[sys].priority[iPos] = -1;
		m_ObsType[sys].Type[iPos]     = -1;
	}
	return true;
}
//NAME: ReadObsHeader()
//reading the file header
//true :successful    false: failed
bool RinexObsLoad::ReadObsHeader(FILE *ft_out)
{
	char buff[MAXLENGTH],*label=buff+60;
	u2 sys = 0;

	while (fgets(buff,MAXLENGTH, m_File) != NULL) {

		if (strlen(buff)<=60) 
		{
			continue;
		}
		else if (strstr(label,"RINEX VERSION / TYPE") != NULL) 
		{
			m_Version = (int)str2num(buff,0,9);
			if(ft_out != NULL)
			{
				fprintf(ft_out,"     3.02           OBSERVATION DATA    Mixed(MIXED)        RINEX VERSION / TYPE\n");
			}
		}
		else if (strstr(label,"SYS / # / OBS TYPES") != NULL)   //Version 3.0
		{
			GetObsType_3(buff);
			if(ft_out != NULL && buff[0] == 'C')
			{
				fprintf(ft_out,"C    6 C2I C7I C6I L2I L7I L6I                              SYS / # / OBS TYPES\n");
			}
		}
		else if(strstr(label,"# / TYPES OF OBSERV") != NULL)
		{
			GetObsType_2(buff);
			if(ft_out != NULL)
			{
				fprintf(ft_out,"C    6 C2I C7I C6I L2I L7I L6I                              SYS / # / OBS TYPES\n");
			}
		}
		else if (strstr(label,"END OF HEADER") != NULL) 
		{
			if(ft_out != NULL)
			{
				fprintf(ft_out,"%s",buff);
			}
			return true;
		}
		else if(strstr(label,"PRN / # OF OBS") != NULL)
		{
			continue;
		}
		else
		{
			if(ft_out != NULL)
			{
				fprintf(ft_out,"%s",buff);
			}
		}
	}
	return false;
}
bool RinexObsLoad::GetObsData_3(gtime_t time, GNSSDATA &ObsData)
{
	char buff[MAXLENGTH];
	u2 i,j,k,SatNum;
	double ep[6];
	GNSSDATA DataTemp = {0};
	double val[MAXENTRY_NUMBER];
	
	while(fgets(buff,MAXLENGTH,m_File) != NULL)
	{
		if(buff[0] == '>')
		{
			ep[0] = str2num(buff, 2, 4); //year
			ep[1] = str2num(buff, 7, 2); //month
			ep[2] = str2num(buff, 10,2); //day
			ep[3] = str2num(buff, 13,2); //hour
			ep[4] = str2num(buff, 16,2); //minute
			ep[5] = str2num(buff, 19,10);
			DataTemp.ObsTime = epoch2time(ep);

			SatNum = (int)str2num(buff,33,2);
			
			int sys,Sat=0;
			u2 priority[2][3]={0},NAVSYS=SYS_NONE;
			for(i = 0; i < SatNum; i++)
			{
				memset(priority, 0x00, sizeof(u2)*6);
				if(fgets(buff,MAXLENGTH,m_File) != NULL)
				{
					if((sys = JudgeSys(buff[0], m_NavSys, &NAVSYS)) == -1)
					{
						continue;
					}
					else
					{
						DataTemp.ObsData[Sat].prn = (int)str2num(buff,1,2);
						u2 PrnNumb = SatNo(DataTemp.ObsData[Sat].prn, NAVSYS);
						if(PrnNumb == 0)
						{
							continue;
						}
						DataTemp.ObsData[Sat].prn  = PrnNumb;
						int prn_1 = PrnNumb - 1;
						for(j=0,k=3; j<m_ObsType[sys].ntype; j++,k+=16)
						{	
							val[j] = str2num(buff,k,14);
							if(0 == m_ObsType[sys].valid[j] || fabs(val[j]) < 1e-5)
							{
								continue;
							}
							int freq = m_ObsType[sys].Frequncy[j] -1;
							
							if(freq > 2 || freq < 0)
							{
								continue;
							}

							if( m_ObsType[sys].Type[j] == 0) //code
							{	
// 								if(m_ObsType[sys].priority[j] >= priority[0][freq])
								{
									if (m_ObsType[sys].bP1C1[j] || m_ObsType[sys].bP2C2[j])
									{
										DataTemp.ObsData[Sat].C[freq] = val[j] + _bias.BDSRecBias[_RecType][prn_1][freq];
										if(freq == 0)
										{
											DataTemp.ObsData[Sat].C[freq] += _bias.P1C1[prn_1];
										}
										else if(freq == 1)
										{
											DataTemp.ObsData[Sat].C[freq] += _bias.P2C2[prn_1];
										}
										else
										{
											DataTemp.ObsData[Sat].C[freq] += _bias.P3C3[prn_1];
										}
									}
									else 
									{
										DataTemp.ObsData[Sat].P[freq] = val[j] + _bias.BDSRecBias[_RecType][prn_1][freq];
									}			
									priority[0][freq] = m_ObsType[sys].priority[j];
								}					
							}
							else if(m_ObsType[sys].Type[j] == 1)  //phase
							{
								if(m_ObsType[sys].priority[j] >= priority[1][freq])
								{
									DataTemp.ObsData[Sat].L[freq] = val[j];
									priority[1][freq] = m_ObsType[sys].priority[j];
								}
							}					
						}
						Sat++;
					}
				}
				else 
				{
					return false;
				}
			}//SatNum
			ASSERT(Sat < MAXCHN);
			
			DataTemp.NumSats = Sat;
			double detaT;
			detaT = timediff(time, DataTemp.ObsTime);
			if(detaT < -PRECISION)
			{
				qsort(DataTemp.ObsData, DataTemp.NumSats, sizeof(ObsData_t), cmpobs);      //sort by sat
				memcpy(&m_NowObsData, &DataTemp, sizeof(GNSSDATA));
				break;
			}
			else if(fabs(detaT) < PRECISION)
			{
				qsort(DataTemp.ObsData, DataTemp.NumSats, sizeof(ObsData_t), cmpobs);      //sort by sat
				memcpy(&m_NowObsData, &DataTemp, sizeof(GNSSDATA));
				memcpy(&ObsData, &DataTemp, sizeof(GNSSDATA));
				return true;
			}			
		}
		else 
		{
			continue;
		}
	}
	return false;
}
bool RinexObsLoad::GetObsData_2(gtime_t time, GNSSDATA &ObsData)
{
	char buff[MAXLENGTH];
	u2 i,j,k,SatNum;
	double ep[6];
	
	while(fgets(buff,MAXLENGTH,m_File) != NULL)
	{
		GNSSDATA DataTemp = {0};
		ep[0] = str2num(buff, 1, 2) + 2000; //year
		ep[1] = str2num(buff, 4, 2); //month
		ep[2] = str2num(buff, 7,2); //day
		ep[3] = str2num(buff, 10,2); //hour
		ep[4] = str2num(buff, 13,2); //minute
		ep[5] = str2num(buff, 16,10);
		DataTemp.ObsTime = epoch2time(ep);
		SatNum = (int)str2num(buff,30,2);
		int Flag[MAXCHN];
		memset(Flag, 0x00, sizeof(int) * MAXCHN);
		u2 NAVSYS=SYS_NONE;
		char aux[4];
		int Sat=0;
		for (i=0,j=32; i<SatNum; i++,j+=3) 
		{
			if (j>=68) 
			{
				if(fgets(buff,MAXLENGTH,m_File) == NULL) 
				{
					break;
				}
				j=32;
			}
			getstr(aux,buff, j, 3);
			if((JudgeSys(aux[0], m_NavSys, &NAVSYS)) == -1)
			{
				Flag[i] = 1;
				continue;
			}
			int prn = Satid2No(aux);
			if(prn == 0)
			{
				Flag[i] = 1;
				continue;
			}
			Flag[i] = 0;
			DataTemp.ObsData[Sat].prn = prn;
			Sat++;
		}

		u2 SatN, priority[2][3]={0};
		Sat = 0;
		for(SatN = 0; SatN < SatNum; SatN++)
		{			
			if(fgets(buff,MAXLENGTH,m_File) == NULL) 
			{
				break;
			}
			if(1 == Flag[SatN])
			{
				for (i=0,j=0; i< m_ObsType[0].ntype; i++,j+=16) 
				{
					if (j>=80) 
					{
						if(fgets(buff,MAXLENGTH,m_File) == NULL) 
							break;
						j=0;
					}
				}
				continue;
			}
			else
			{
				int prn_1 = DataTemp.ObsData[Sat].prn - 1;
				memset(priority, 0x00, sizeof(u2)*6);
				double val[MAXCHN] = {0.0};

				for (i=0,j=0; i< m_ObsType[0].ntype; i++,j+=16) 
				{
					if (j>=80) 
					{
						if(fgets(buff,MAXLENGTH,m_File) == NULL) 
							break;
						j=0;
					}
					val[i]=str2num(buff,j,14);
					if(0 == m_ObsType[0].valid[i] || fabs(val[i]) < 1e-5)
					{
						continue;
					}
					int freq = m_ObsType[0].Frequncy[i] -1;

					if(freq > 2 || freq < 0)
					{
						continue;
					}
					if( m_ObsType[0].Type[i] == 0)
					{	
// 						if(m_ObsType[0].priority[i] >= priority[0][freq])
						{
							if(m_ObsType[0].bP1C1[i] || m_ObsType[0].bP2C2[i])
							{
								DataTemp.ObsData[Sat].C[freq] = val[i] + _bias.BDSRecBias[_RecType][prn_1][freq];
								if(freq == 0)
								{
									DataTemp.ObsData[Sat].C[freq] += _bias.P1C1[prn_1];
								}
								else if(freq == 1)
								{
									DataTemp.ObsData[Sat].C[freq] += _bias.P2C2[prn_1];
								}
								else
								{
									DataTemp.ObsData[Sat].C[freq] += _bias.P3C3[prn_1];
								}
							}
							else
							{
								DataTemp.ObsData[Sat].P[freq] = val[i] + _bias.BDSRecBias[_RecType][prn_1][freq];
							}
							priority[0][freq] = m_ObsType[0].priority[i];
						}					
					}
					else if(m_ObsType[0].Type[i] == 1)
					{
						if(m_ObsType[0].priority[i] >= priority[1][freq])
						{
							DataTemp.ObsData[Sat].L[freq] = val[i];
							priority[1][freq] = m_ObsType[0].priority[i];
						}
					}					
				}  //ntype
				Sat++;
			}	
		} //SatNum
		ASSERT(Sat < MAXCHN);

		DataTemp.NumSats = Sat;
		double detaT;
		detaT = timediff(time, DataTemp.ObsTime);
		if(detaT < -PRECISION)
		{
			qsort(DataTemp.ObsData, DataTemp.NumSats, sizeof(ObsData_t), cmpobs);      //sort by sat
			memcpy(&m_NowObsData, &DataTemp, sizeof(GNSSDATA));
			break;
		}
		else if(fabs(detaT) < PRECISION)
		{
			qsort(DataTemp.ObsData, DataTemp.NumSats, sizeof(ObsData_t), cmpobs);      //sort by sat
			memcpy(&m_NowObsData, &DataTemp, sizeof(GNSSDATA));
			memcpy(&ObsData, &DataTemp, sizeof(GNSSDATA));
			return true;
		}
	}
	return false;
}
bool RinexObsLoad::GetObsData(gtime_t time, GNSSDATA &ObsData)
{
	double detaT;
	if(m_File == NULL)
	{
		return false;
	}
	detaT = timediff(time, m_NowObsData.ObsTime);
	if(m_NowObsData.NumSats == 0) //the first epoch
	{
		if(3 == m_Version)
		{
			return GetObsData_3(time, ObsData);
		}
		else if(2 == m_Version)
		{
			return GetObsData_2(time, ObsData);
		}
		else
		{
			return false;
		}
	}
	else
	{
		if(detaT < -PRECISION) // small than 0
		{
			return false;
		}
		else if(detaT < PRECISION) //the correct time
		{
			memcpy(&ObsData, &m_NowObsData, sizeof(GNSSDATA));
			return true;
		}
		else
		{
			if(3 == m_Version)
			{
			     return GetObsData_3(time, ObsData);
			}
			else if(2 == m_Version)
			{
				return GetObsData_2(time, ObsData);
			}
			else
			{
				return false;
			}

		}

	}	
	return true;
}