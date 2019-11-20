#include "../Common/FileLoad.h"
#include "../Common/BaseFunction.h"
#include "../Common/debug.h"

RinexNavLoad::RinexNavLoad()
{
	m_File    = NULL;
	m_Version = 0;
}
RinexNavLoad::~RinexNavLoad()
{
	;
}

bool RinexNavLoad::Initial(char chFile[], u2 NavSys)
{
	if(m_File != NULL)
	{
		fclose(m_File);
		m_File = NULL;
	}
	m_File = fopen(chFile,"r");
	m_SYS  = NavSys;
	if(m_File == NULL)
		return false;
	return  ReadNavHeader();
}
bool RinexNavLoad::ReadNavHeader()
{
	char buff[MAXLENGTH],*label=buff+60;

	while (fgets(buff,MAXLENGTH,m_File) != NULL) {

		if (strlen(buff)<=60) 
		{
			continue;
		}
		else if (strstr(label,"RINEX VERSION / TYPE") != NULL) 
		{
			m_Version = (int)str2num(buff,0,9);
		}
		else if (strstr(label,"ION ALPHA") != NULL) 
		{
			xstrReplace(buff,'D','E',4);
			sscanf(buff,"%lf %lf %lf %lf",&m_Alpha[0], &m_Alpha[1], &m_Alpha[2], &m_Alpha[3]);
		}
		else if (strstr(label,"ION BETA") != NULL) 
		{
			xstrReplace(buff,'D','E',4);
			sscanf(buff,"%lf %lf %lf %lf",&m_Beta[0], &m_Beta[1], &m_Beta[2], &m_Beta[3]);
		}
		else if (strstr(label,"END OF HEADER") != NULL) 
		{
			return true;
		}
	}
	return false;
}
bool RinexNavLoad::ReadNavData(GNSSEPH_t &brdc)
{
	char buff[MAXLENGTH],*label=buff+60;
	bool Flag;
	if(m_File == NULL)
	{
		return false;
	}
	if(3 == m_Version)
	{
		Flag = ReadNavData_3(brdc);
		if(m_File != NULL)
		{
			fclose(m_File);
			m_File = NULL;
		}
		 return Flag;
	}
	else if(2 == m_Version)
	{
		Flag = ReadNavData_2(brdc);
		if(m_File != NULL)
		{
			fclose(m_File);
			m_File = NULL;
		}
		return Flag;
	}
	return true;
}

bool RinexNavLoad::ReadNavData_3(GNSSEPH_t &brdc)
{
	char buff[MAXLENGTH],StrSys,strtemp[10];
	u2 i,row=0,Sys_Sys;
	EPHN_t eph = {0},*nav_eph;
	int ep[6],sys;
	double IODE,IODC,TOE=0,week,health,URAindex,temp;

	while (fgets(buff,MAXLENGTH,m_File) != NULL) {
		xstrReplace(buff,'D','E',4);
		switch(row)
		{
			case 0:
				if(buff[0] == 'S' || buff[0] == 'R')
				{
					for(i=0; i<3; i++)   //QZSS
					{
						if(fgets(buff,MAXLENGTH,m_File) == NULL)
						{
							return false;
						}
					}
					break;
				}
				if((sys = JudgeSys(buff[0],SYS_ALL, &Sys_Sys)) == -1)
				{
					for(i=0; i<7; i++)   //QZSS
					{
						if(fgets(buff,MAXLENGTH,m_File) == NULL)
						{
							return false;
						}
					}
					break;
				}
				
				StrSys = buff[0];
				eph.prn = (int)str2num(buff,1,2);
				eph.prn = SatNo(eph.prn, Sys_Sys);
				sscanf(buff,"%s%d%d%d%d%d%d%lf%lf%lf",&strtemp,&ep[0],&ep[1],&ep[2],&ep[3],&ep[4],&ep[5],
						&eph.clock_bias, &eph.clock_drift, &eph.clock_driftrate);
				if(buff[0] == 'E' && ep[4] > 0.0)  //GAL
				{
					for(i=0; i<7; i++)   //QZSS
					{
						if(fgets(buff,MAXLENGTH,m_File) == NULL)
						{
							return false;
						}
					}
					break;
				}
				row++;
				break;
			case 1:
				sscanf(buff,"%lf%lf%lf%lf",&IODE, &eph.Crs, &eph.Delta_n, &eph.M0);
				eph.IODE = (int)IODE;
				row++;
				break;
			case 2:
				sscanf(buff,"%lf%lf%lf%lf",&eph.Cuc, &eph.e, &eph.Cus, &eph.sqrt_A);
				eph.sqrt_A = eph.sqrt_A*eph.sqrt_A;
				row++;
				break;
			case 3:
				sscanf(buff,"%lf%lf%lf%lf",&TOE, &eph.Cic, &eph.OMEGA0, &eph.Cis);
				eph.toes = TOE;
				row++;
				break;
			case 4:
				sscanf(buff,"%lf%lf%lf%lf",&eph.i0, &eph.Crc, &eph.omega, &eph.OMEGADOT);
				row++;
				break;
			case 5:
				sscanf(buff,"%lf%lf%lf%lf",&eph.IDOT, &temp, &week, &temp);
				if(((int)week) < 1000 && StrSys == 'C')
				{
					week += 1356.0;
					eph.TOE = gpst2time((int)week,TOE);
					eph.TOC = eph.TOE;
					eph.TOE.time +=  14;
					eph.TOC.time +=  14;
				}
				else
				{
					eph.TOE = gpst2time((int)week,TOE);
					eph.TOC = eph.TOE;
				}
				
				row++;
				break;
			case 6:
				if(StrSys == 'C')
				{
					sscanf(buff,"%lf%lf%lf%lf",&URAindex, &health, &eph.TGD[0],&eph.TGD[1]);
				}
				else
				{
					eph.TGD[1] = 0.0;
					sscanf(buff,"%lf%lf%lf%lf",&URAindex, &health, &eph.TGD[0], &IODC);
				}
				eph.URAindex = (int)URAindex;
				eph.SVhealth = (int)health;
				eph.IODC     = (int)IODC;
				row++;
				break;
			case 7:
				if(StrSys == 'C')
				{
					sscanf(buff,"%lf%lf",&temp, &IODC);
					eph.IODC     = (int)IODC;
				}
				row=0;
				if(sys< 0 || sys > 3)
				{
					ASSERT(sys>= 0 && sys <4);
				}
				if(brdc.brdc[sys].max <= brdc.brdc[sys].n)
				{
					brdc.brdc[sys].max += 800;
					if ((nav_eph=(EPHN_t *)realloc(brdc.brdc[sys].eph, sizeof(EPHN_t)*brdc.brdc[sys].max)) == NULL)
					{
						ASSERT(1 < 0);
					}
					brdc.brdc[sys].eph = nav_eph;
				}
				memcpy(&brdc.brdc[sys].eph[brdc.brdc[sys].n], &eph, sizeof(EPHN_t));
				brdc.brdc[sys].n++;
				memset(&eph, 0x00, sizeof(EPHN_t));
				break;
			default:
				break;
		}
	}
	return true;
}
bool RinexNavLoad::ReadNavData_2(GNSSEPH_t &brdc)
{
	char buff[MAXLENGTH],StrSys = 'G',strtemp[10];
	u2 i,row=0,Sys_Sys;
	EPHN_t eph = {0},*nav_eph;
	int ep[6],sys = 0;
	double sec,IODE,IODC,TOE=0,week,health,URAindex,temp;

	while (fgets(buff,MAXLENGTH,m_File) != NULL) {
		xstrReplace(buff,'D','E',4);
		switch(row)
		{
		case 0:
			sscanf(buff,"%d%d%d%d%d%d%lf%lf%lf%lf",&eph.prn,&ep[0],&ep[1],&ep[2],&ep[3],&ep[4],&sec,
				&eph.clock_bias, &eph.clock_drift, &eph.clock_driftrate);
			if(ep[0] < 50)
			{
				ep[0] += 2000;
			}
			else
			{
				ep[0] += 1900;
			}
			row++;
			break;
		case 1:
			sscanf(buff,"%lf%lf%lf%lf",&IODE, &eph.Crs, &eph.Delta_n, &eph.M0);
			eph.IODE = (int)IODE;
			row++;
			break;
		case 2:
			sscanf(buff,"%lf%lf%lf%lf",&eph.Cuc, &eph.e, &eph.Cus, &eph.sqrt_A);
			eph.sqrt_A = eph.sqrt_A*eph.sqrt_A;
			row++;
			break;
		case 3:
			sscanf(buff,"%lf%lf%lf%lf",&TOE, &eph.Cic, &eph.OMEGA0, &eph.Cis);
			eph.toes = TOE;
			row++;
			break;
		case 4:
			sscanf(buff,"%lf%lf%lf%lf",&eph.i0, &eph.Crc, &eph.omega, &eph.OMEGADOT);
			row++;
			break;
		case 5:
			sscanf(buff,"%lf%lf%lf%lf",&eph.IDOT, &temp, &week, &temp);
			if(((int)week) < 1000 && StrSys == 'C')
			{
				week += 1356.0;
				eph.TOE = gpst2time((int)week,TOE);
				eph.TOC = eph.TOE;
				eph.TOE.time +=  14;
				eph.TOC.time +=  14;
			}
			else
			{
				eph.TOE = gpst2time((int)week,TOE);
				eph.TOC = eph.TOE;
			}

			row++;
			break;
		case 6:
			if(StrSys == 'C')
			{
				sscanf(buff,"%lf%lf%lf%lf",&URAindex, &health, &eph.TGD[0],&eph.TGD[1]);
			}
			else
			{
				eph.TGD[1] = 0.0;
				sscanf(buff,"%lf%lf%lf%lf",&URAindex, &health, &eph.TGD[0], &IODC);
			}
			eph.URAindex = (int)URAindex;
			eph.SVhealth = (int)health;
			eph.IODC     = (int)IODC;
			row++;
			break;
		case 7:
			if(StrSys == 'C')
			{
				sscanf(buff,"%lf%lf",&temp, &IODC);
				eph.IODC     = (int)IODC;
			}
			row=0;
			if(brdc.brdc[sys].max <= brdc.brdc[sys].n)
			{
				brdc.brdc[sys].max += 800;
				if ((nav_eph=(EPHN_t *)realloc(brdc.brdc[sys].eph, sizeof(EPHN_t)*brdc.brdc[sys].max)) == NULL)
				{
					ASSERT(1 < 0);
				}
				brdc.brdc[sys].eph = nav_eph;
			}
			memcpy(&brdc.brdc[sys].eph[brdc.brdc[sys].n], &eph, sizeof(EPHN_t));
			brdc.brdc[sys].n++;
			memset(&eph, 0x00, sizeof(EPHN_t));
			break;
		default:
			break;
		}
	}
	return true;
}