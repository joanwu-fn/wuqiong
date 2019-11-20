#include "pNetcontrol.h"

void Netcontrol::OutPutMP()
{
	int i,j,StaNo;
	int StartS = MAXPRNGPS, EndS = MAXPRNGPS;

	/*                   OutPut Header                     */
	fprintf(m_Option.ft_Static[0],"Station	");
	fprintf(m_Option.ft_Static[1],"Station	");
	fprintf(m_Option.ft_Static[2],"Station	");
	if(m_Option.NavSys & SYS_GPS)
	{
		StartS = 0;
		for(i = 0; i < MAXPRNGPS; i++)
		{
			fprintf(m_Option.ft_Static[0]," G%02d	",i+1);
			fprintf(m_Option.ft_Static[1]," G%02d	",i+1);
			fprintf(m_Option.ft_Static[2],"          G%02d      	",i+1);
		}
	}
	if(m_Option.NavSys & SYS_BDS)
	{
		EndS = MAXPRNGPS + MAXPRNBDS;
		for(i = 0; i < MAXPRNBDS; i++)
		{
			fprintf(m_Option.ft_Static[0]," C%02d	",i+1);
			fprintf(m_Option.ft_Static[1]," C%02d	",i+1);
			fprintf(m_Option.ft_Static[2],"          G%02d      	",i+1);
		}
	}

	fprintf(m_Option.ft_Static[0],"\n");
	fprintf(m_Option.ft_Static[1],"\n");
	fprintf(m_Option.ft_Static[2],"\n");
	/*                   OutPut Header                     */

	for(StaNo = 0; StaNo < m_Option.StaNum; StaNo++)
	{
		for(i = StartS; i < EndS; i++)
		{
			for(j = m_GNSSInfo[StaNo].MP_S[i].LastEpoch[0] - 1; j < m_GNSSInfo[StaNo].MP_S[i].LastEpoch[0]+1; j++)
			{
				m_GNSSInfo[StaNo].MP_S[i].MP_1[j] = 0.0;
			}
			for(j = m_GNSSInfo[StaNo].MP_S[i].LastEpoch[1] - 1; j < m_GNSSInfo[StaNo].MP_S[i].LastEpoch[1]+1; j++)
			{
				m_GNSSInfo[StaNo].MP_S[i].MP_2[j] = 0.0;
			}
		}

		for(j = 0; j < m_TTnumTEQC; j++)
		{
			for(i = StartS; i < EndS; i++)
			{
				double Temp = m_GNSSInfo[StaNo].MP_S[i].MP_1[j];
				fprintf_s(m_Option.Station[StaNo].ft_ewl, "%10.2f	", Temp);
				Temp = m_GNSSInfo[StaNo].MP_S[i].MP_2[j];
				fprintf_s(m_Option.Station[StaNo].ft_mw, "%10.2f	", Temp);
			}
			fprintf(m_Option.Station[StaNo].ft_ewl,"\n");
			fprintf(m_Option.Station[StaNo].ft_mw,"\n");
		}

		for(j = 0; j < MAXSAT; j++)
		{
			free(m_GNSSInfo[StaNo].MP_S[j].MP_1);
			free(m_GNSSInfo[StaNo].MP_S[j].MP_2);
			m_GNSSInfo[StaNo].MP_S[j].MP_1 = NULL;
			m_GNSSInfo[StaNo].MP_S[j].MP_2 = NULL;
		}


		fprintf(m_Option.ft_Static[0],"%s	",m_GNSSInfo[StaNo].Station[StaNo].Name);
		fprintf(m_Option.ft_Static[1],"%s	",m_GNSSInfo[StaNo].Station[StaNo].Name);
		fprintf(m_Option.ft_Static[2],"%s	",m_GNSSInfo[StaNo].Station[StaNo].Name);
		for(i = StartS; i < EndS; i++)
		{
			if(m_GNSSInfo[StaNo].MP_S[i].Std[1] > 0.0)
			{
				fprintf(m_Option.ft_Static[0],"%4.2f	",m_GNSSInfo[StaNo].MP_S[i].Std[0] / m_GNSSInfo[StaNo].MP_S[i].Std[2]);
			}
			else
			{
				fprintf(m_Option.ft_Static[0],"%4.2f	",0.0);
			}
			if(m_GNSSInfo[StaNo].MP_S[i].Std[3] > 0.0)
			{
				fprintf(m_Option.ft_Static[1],"%4.2f	",m_GNSSInfo[StaNo].MP_S[i].Std[3] / m_GNSSInfo[StaNo].MP_S[i].Std[5]);
			}
			else
			{
				fprintf(m_Option.ft_Static[1],"%4.2f	",0.0);
			}
			fprintf(m_Option.ft_Static[2],"%5d %5d / %5d	",
				m_GNSSInfo[StaNo].SSat[i].CSNum[2],m_GNSSInfo[StaNo].SSat[i].CSNum[1], m_GNSSInfo[StaNo].SSat[i].CSNum[0]);
		}
		fprintf(m_Option.ft_Static[0],"\n");
		fprintf(m_Option.ft_Static[1],"\n");
		fprintf(m_Option.ft_Static[2],"\n");
	}

	fprintf(m_Option.ft_Static[3],"Station  Continuity(%%)     CycleSlip    MP1(m)  MP2(m)\n");
	for(StaNo = 0; StaNo < m_Option.StaNum; StaNo++)
	{
		double conti = (m_GNSSInfo[StaNo].UseFulNum * 1.0 / (m_TTnumTEQC * 1.0)) * 100.0;
		int cycleslip[3] = {0}, Visible[2] = {0};
		double MP[2] = {0};
		int MPNum[2] = {0};
		double partconti;
		for(i = StartS; i < EndS; i++)
		{
			cycleslip[0] += m_GNSSInfo[StaNo].SSat[i].CSNum[0];
			cycleslip[1] += m_GNSSInfo[StaNo].SSat[i].CSNum[1];	
			cycleslip[2] += m_GNSSInfo[StaNo].SSat[i].CSNum[2];	
			Visible[0]   += m_GNSSInfo[StaNo].Visible[i][0];
			Visible[1]   += m_GNSSInfo[StaNo].Visible[i][1];
			if(fabs(m_GNSSInfo[StaNo].MP_S[i].Std[0]) > PRECISION)
			{
				MP[0] += m_GNSSInfo[StaNo].MP_S[i].Std[0] / m_GNSSInfo[StaNo].MP_S[i].Std[2];
				MPNum[0]++;
			}
			if(fabs(m_GNSSInfo[StaNo].MP_S[i].Std[3]) > PRECISION)
			{
				MP[1] += m_GNSSInfo[StaNo].MP_S[i].Std[3] / m_GNSSInfo[StaNo].MP_S[i].Std[5];
				MPNum[1]++;
			}
		}
		partconti    = (Visible[1] * 1.0 / ((Visible[0] + 1) * 1.0)) * 100.0;
		if(MPNum[0] > 0)
		{
			MP[0] = MP[0] / MPNum[0];
		}
		else
		{
			MP[0] = 0.0;
		}
		if(MPNum[1] > 0)
		{
			MP[1] = MP[1] / MPNum[1];
		}
		else
		{
			MP[1] = 0.0;
		}
		fprintf(m_Option.ft_Static[3],"%s    %5.2f   %5.2f  %5d %5d / %7d  %5.2f   %5.2f\n",
			m_Option.Station[StaNo].Name, 100 - conti, 100 - partconti, cycleslip[2], cycleslip[1], cycleslip[0], MP[0], MP[1]);
	}

	/*---------------------------------------œÍœ∏–≈œ¢--------------------------------*/
	fprintf(m_Option.ft_Static[3],"\n\n\nSat	");
	for(StaNo = 0; StaNo < m_Option.StaNum; StaNo++)
	{
		fprintf(m_Option.ft_Static[3],"          %s	           ",m_Option.Station[StaNo].Name);
	}
	fprintf(m_Option.ft_Static[3],"\n     ");
	for(StaNo = 0; StaNo < m_Option.StaNum; StaNo++)
	{
		fprintf(m_Option.ft_Static[3]," Lose  Part  CSlip           ");
	}
	for(j = StartS; j < EndS; j++)
	{
		if(j < MAXPRNGPS)
		{
			fprintf(m_Option.ft_Static[3],"\nG%02d	",j+1);
		}
		else
		{
			fprintf(m_Option.ft_Static[3],"\nC%02d	",j-31);
		}
		for(StaNo = 0; StaNo < m_Option.StaNum; StaNo++)
		{
			double Continuity[3];
			if(m_GNSSInfo[StaNo].Visible[j][0] == 0)
			{
				m_GNSSInfo[StaNo].Visible[j][0] += 1;
			}
			if(m_GNSSInfo[StaNo].SSat[j].CSNum[0] == 0)
			{
				m_GNSSInfo[StaNo].SSat[j].CSNum[0] += 1;
			}
			Continuity[0] = m_GNSSInfo[StaNo].Visible[j][1] * 1.0 / (m_GNSSInfo[StaNo].Visible[j][0] * 1.0 );
			Continuity[1] = m_GNSSInfo[StaNo].SSat[j].CSNum[2] * 1.0 / (m_GNSSInfo[StaNo].SSat[j].CSNum[0] * 1.0);
			Continuity[2] = m_GNSSInfo[StaNo].SSat[j].CSNum[1] * 1.0 / (m_GNSSInfo[StaNo].SSat[j].CSNum[0] * 1.0);
			
			fprintf(m_Option.ft_Static[3],"%6.2f %6.2f %6.2f   ***   ",
				100.0 - Continuity[0]*100,Continuity[1]*100,Continuity[2]*100);
		}
	}
	
}



void Netcontrol::OutPutDCB()
{
	int StaNo,i,Ref,ToTal=0,RefNum=0;
	double DetaDCB[14],Ave=0;
	for(StaNo = 0; StaNo < m_Option.StaNum; StaNo++)
	{
		GetAverageDCB(m_GNSSInfo[StaNo]);
		ToTal = 0;
		for(i = 0; i < 14; i++)
		{
			if(m_GNSSInfo[StaNo].SSat[32+i].DCBNum > 300)
			{
				ToTal++;
			}
			else
			{
				m_GNSSInfo[StaNo].SSat[32+i].DetaDCB = 0.0;
			}
		}
		if(ToTal > RefNum)
		{
			Ref = StaNo;
			RefNum = ToTal;
		}
	}
	//////////////////////////////////////////////////////////////////////////
	RefNum = 0;
	for(i = 0; i < 14; i++)
	{
		DetaDCB[i] = m_GNSSInfo[Ref].SSat[i+32].DetaDCB;
		Ave += m_GNSSInfo[Ref].SSat[i+32].DetaDCB;
		RefNum++;
	}
	Ave = Ave / RefNum;
	//////////////////////////////////////////////////////////////////////////
	double DetaTemp = 0;
	char FileDir[300];
	sprintf(FileDir,"%s\\EstDcb%02d%03d.txt",m_Option.Directory.DCBDir, m_Option.ProcessYear-2000, m_Option.ProcessDoY);
	FILE *fp = fopen(FileDir,"w");
	for(StaNo = 0; StaNo < m_Option.StaNum; StaNo++)
	{
		ToTal = 0;
		RefNum = 0;
		for(i = 0; i < 14; i++)
		{
			if(fabs(m_GNSSInfo[StaNo].SSat[32+i].DetaDCB) > PRECISION && fabs(DetaDCB[i]) > PRECISION)
			{
				DetaTemp += (m_GNSSInfo[StaNo].SSat[32+i].DetaDCB - DetaDCB[i]);
				RefNum++;		
			}
		}
		if(RefNum > 0)
		{
			DetaTemp = DetaTemp / RefNum;
			fprintf(fp,"%s	",m_GNSSInfo[StaNo].Station[StaNo].Name);
			for(i = 0; i < 14; i++)
			{
				if(fabs(m_GNSSInfo[StaNo].SSat[32+i].DetaDCB) > PRECISION)
				{
					fprintf(fp,"%10.3f	",m_GNSSInfo[StaNo].SSat[32+i].DetaDCB - Ave - DetaTemp);		
				}
				else
				{
					fprintf(fp,"%10.3f	",0.0);
				}
			}
			fprintf(fp,"\n");
		}

		for(i = 0; i < MAXSAT; i++)
		{
			free(m_GNSSInfo[StaNo].MP_S[i].MP_1);
			free(m_GNSSInfo[StaNo].MP_S[i].MP_2);
			m_GNSSInfo[StaNo].MP_S[i].MP_1 = NULL;
			m_GNSSInfo[StaNo].MP_S[i].MP_2 = NULL;
		}
	}
	fclose(fp);
}