#include "../PreProcess/PreProcess.h"
double MPValue(const ObsData_t ObsData, enum MPTYPE type)
{
	u2 SYS,SysNum;
	char StrSys;
	double MP;

	if( Prn2Sys(ObsData.prn, SYS, SysNum, StrSys) <= 0)
		return 0.0;
	if(type == MP1)
	{
		if(ObsData.P[0] < PRECISION || fabs(ObsData.L[0]) < PRECISION || fabs(ObsData.L[1]) < PRECISION)
		{
			return 0.0;
		}
		double C1 = (G_Fre[SysNum][0] * G_Fre[SysNum][0] + G_Fre[SysNum][1] * G_Fre[SysNum][1]) / (G_Fre[SysNum][0] * G_Fre[SysNum][0] - G_Fre[SysNum][1] * G_Fre[SysNum][1]);
		double C2 = 1 - C1;
		MP = ObsData.P[0] - ( C1 * ObsData.L[0] * G_Lam[SysNum][0] + C2 * ObsData.L[1] * G_Lam[SysNum][1]);
	}
	else if(type == MP2)
	{
		if(ObsData.P[1] < PRECISION || fabs(ObsData.L[0]) < PRECISION || fabs(ObsData.L[1]) < PRECISION)
		{
			return 0.0;
		}
		double C1 = 2 * G_Fre[SysNum][0] * G_Fre[SysNum][0] / (G_Fre[SysNum][0] * G_Fre[SysNum][0] - G_Fre[SysNum][1] * G_Fre[SysNum][1]);
		double C2 = 1 - C1;
		MP = ObsData.P[1] - ( C1 * ObsData.L[0] * G_Lam[SysNum][0] + C2 * ObsData.L[1] * G_Lam[SysNum][1]);
	}
	else if(type == MP3)
	{
		if(ObsData.P[2] < PRECISION || fabs(ObsData.L[0]) < PRECISION || fabs(ObsData.L[2]) < PRECISION)
		{
			return 0.0;
		}
		double C1 = (G_Fre[SysNum][0] * G_Fre[SysNum][0] + G_Fre[SysNum][2] * G_Fre[SysNum][2]) / (G_Fre[SysNum][0] * G_Fre[SysNum][0] - G_Fre[SysNum][2] * G_Fre[SysNum][2]);
		double C2 = 1 - C1;
		MP = ObsData.P[2] - ( C1 * ObsData.L[0] * G_Lam[SysNum][0] + C2 * ObsData.L[2] * G_Lam[SysNum][2]);
	}
	else
	{
		return 0.0;
	}
	return MP;
}
void GetSTD(double *Data, int Numb, double *STD, int type)
{
	double Ave = 0;
	STD[0] = 0.0;
	STD[1] = 0.0;
	int i;
	if(Numb < 1)
	{
		return;
	}
	for(i = 0; i < Numb; i++)
	{
		Ave +=Data[i];
	}
	Ave /= Numb;
	for(i = 1; i < Numb; i++)
	{
		STD[0] += (Data[i] -Ave) * (Data[i] -Ave);
		STD[1] += (Data[i] -Data[i-1]) * (Data[i] - Data[i-1]);				
	}

	if(type == 0)
	{			
		for(i = 0; i < Numb; i++)
		{
			Data[i] = Data[i] - Ave;
		}
	}
	STD[0] = sqrt(STD[0] / Numb);
	STD[1] = sqrt(STD[1] / Numb);
	STD[1] = sqrt(fabs(STD[0] * STD[0] - STD[1] * STD[1] / 2.0));
}
void SetZero(double *Data, int Numb)
{
	for(int i = 0; i < Numb; i++)
	{
		Data[i] = 0.0;
	}
}
void GetMP(const GNSSDATA Obs, GNSSInfo &PInfo)
{
	int i,j,prn_1,NowEpo,Last;
	double Std[2];
	NowEpo  = PInfo.EpoNum;
	for(i = 0; i < Obs.NumSats; i++)
	{
		prn_1 = Obs.ObsData[i].prn - 1;
		
		PInfo.MP_S[prn_1].MP_1[NowEpo] = MPValue(Obs.ObsData[i], MP1);
		PInfo.MP_S[prn_1].MP_2[NowEpo] = MPValue(Obs.ObsData[i], MP2);
	
		if(PInfo.SSat[prn_1].CSlip[0] == false && PInfo.SSat[prn_1].CSlip[1] == false)      // there is no cycle slip
		{
			if(fabs(PInfo.MP_S[prn_1].MP_1[NowEpo]) < PRECISION)
			{
				Last = PInfo.MP_S[prn_1].LastEpoch[0];
				if( (NowEpo - Last) > MAX_TEQC_EPOCH)
				{
					GetSTD(&PInfo.MP_S[prn_1].MP_1[Last], NowEpo - Last,Std , 0);
					//printf("%3d %10.3f\n",prn_1+1,Std);

					PInfo.MP_S[prn_1].Std[0] += Std[0] * (NowEpo - Last);
					PInfo.MP_S[prn_1].Std[1] += Std[1] * (NowEpo - Last);
					PInfo.MP_S[prn_1].Std[2] += (NowEpo - Last);
				}
				else
				{
					SetZero(&PInfo.MP_S[prn_1].MP_1[Last], NowEpo - Last);
				}
				PInfo.MP_S[prn_1].LastEpoch[0] = NowEpo + 1;
			}
			if(fabs(PInfo.MP_S[prn_1].MP_2[NowEpo]) < PRECISION)
			{
				Last = PInfo.MP_S[prn_1].LastEpoch[1];
				if( (NowEpo - Last) > MAX_TEQC_EPOCH)
				{
					GetSTD(&PInfo.MP_S[prn_1].MP_2[Last], NowEpo - Last, Std, 1);
					//printf("%3d %10.3f\n",prn_1+1,Std);

					PInfo.MP_S[prn_1].Std[3] += Std[0] * (NowEpo - Last);
					PInfo.MP_S[prn_1].Std[4] += Std[1] * (NowEpo - Last);
					PInfo.MP_S[prn_1].Std[5] += (NowEpo - Last);
				}
				else
				{
					//SetZero(&PInfo.MP_S[prn_1].MP_2[Last], NowEpo - Last);
				}
				PInfo.MP_S[prn_1].LastEpoch[1] = NowEpo + 1;
			}
			
		}
		else     // there is cycle slip occur
		{		
			Last = PInfo.MP_S[prn_1].LastEpoch[0];
			if( (NowEpo - PInfo.MP_S[prn_1].LastEpoch[0]) > MAX_TEQC_EPOCH)
			{				
				GetSTD(&PInfo.MP_S[prn_1].MP_1[Last], NowEpo - Last, Std, 0);
				//printf("%3d %10.3f\n",prn_1+1,Std);

				PInfo.MP_S[prn_1].Std[0] += Std[0] * (NowEpo - Last);
				PInfo.MP_S[prn_1].Std[1] += Std[1] * (NowEpo - Last);
				PInfo.MP_S[prn_1].Std[2] += (NowEpo - Last);
			}
			else
			{
				SetZero(&PInfo.MP_S[prn_1].MP_1[Last], NowEpo - Last);
			}
			if(fabs(PInfo.MP_S[prn_1].MP_1[NowEpo]) < PRECISION)
			{
				PInfo.MP_S[prn_1].LastEpoch[0] = NowEpo +1;
			}
			else
			{
				PInfo.MP_S[prn_1].LastEpoch[0] = NowEpo;
			}

			Last = PInfo.MP_S[prn_1].LastEpoch[1];
			if( (NowEpo - PInfo.MP_S[prn_1].LastEpoch[1]) > MAX_TEQC_EPOCH)
			{			
				GetSTD(&PInfo.MP_S[prn_1].MP_2[Last], NowEpo - Last, Std,  1);
				//printf("%3d %10.3f\n",prn_1+1,Std);

				PInfo.MP_S[prn_1].Std[3] += Std[0] * (NowEpo - Last);
				PInfo.MP_S[prn_1].Std[4] += Std[1] * (NowEpo - Last);
				PInfo.MP_S[prn_1].Std[5] += (NowEpo - Last);
			}
			else
			{
				//SetZero(&PInfo.MP_S[prn_1].MP_2[Last], NowEpo - Last);
			}

			if(fabs(PInfo.MP_S[prn_1].MP_2[NowEpo]) < PRECISION)
			{
				PInfo.MP_S[prn_1].LastEpoch[1] = NowEpo +1;
			}
			else
			{
				PInfo.MP_S[prn_1].LastEpoch[1] = NowEpo;
			}		
		}
	}
	j = 0;
	for(prn_1 = 0; prn_1 < MAXSAT; prn_1++)
	{
		if((prn_1 + 1) == Obs.ObsData[j].prn)
		{
			j++;
		}
		else
		{
			Last = PInfo.MP_S[prn_1].LastEpoch[0];
			if( (NowEpo - PInfo.MP_S[prn_1].LastEpoch[0]) > MAX_TEQC_EPOCH)
			{		
				GetSTD(&PInfo.MP_S[prn_1].MP_1[Last], NowEpo - Last, Std, 0);
				//printf("%3d %10.3f\n",prn_1+1,Std);

				PInfo.MP_S[prn_1].Std[0] += Std[0] * (NowEpo - Last);
				PInfo.MP_S[prn_1].Std[1] += Std[1] * (NowEpo - Last);
				PInfo.MP_S[prn_1].Std[2] += (NowEpo - Last);
			}
			else
			{
				SetZero(&PInfo.MP_S[prn_1].MP_1[Last], NowEpo - Last);
			}
			PInfo.MP_S[prn_1].LastEpoch[0] = NowEpo + 1;

			Last = PInfo.MP_S[prn_1].LastEpoch[1];
			if( (NowEpo - PInfo.MP_S[prn_1].LastEpoch[1]) > MAX_TEQC_EPOCH)
			{			
				GetSTD(&PInfo.MP_S[prn_1].MP_2[Last], NowEpo - Last, Std, 1);
				//printf("%3d %10.3f\n",prn_1+1,Std);

				PInfo.MP_S[prn_1].Std[3] += Std[0] * (NowEpo - Last);
				PInfo.MP_S[prn_1].Std[4] += Std[1] * (NowEpo - Last);
				PInfo.MP_S[prn_1].Std[5] += (NowEpo - Last);
			}
			else
			{
				//SetZero(&PInfo.MP_S[prn_1].MP_2[Last], NowEpo - Last);
			}
			PInfo.MP_S[prn_1].LastEpoch[1] = NowEpo + 1;
		}
	}
}
void TEQC(const GNSSDATA Obs, GNSSInfo &PInfo)
{
	DetectCS_GF2(Obs, PInfo.SSat);
	int i,prn_1;
	for(i = 0; i < Obs.NumSats; i++)
	{
		prn_1 = Obs.ObsData[i].prn - 1;
		PInfo.SSat[prn_1].CSNum[0] += 1;
		if(PInfo.SSat[prn_1].CSlip[0] == true || PInfo.SSat[prn_1].CSlip[0] == true)
		{
			PInfo.SSat[prn_1].CSNum[1] += 1;
		}
		if(fabs(Obs.ObsData[i].L[0]) < 1e-5 || fabs(Obs.ObsData[i].L[1]) < 1e-5)
		{
			PInfo.SSat[prn_1].CSNum[2] += 1;
		}
	}
	GetMP(Obs, PInfo);
	PInfo.EpoNum++;
	
}