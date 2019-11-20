#include "../PreProcess/PreProcess.h"
void EstDCB(const GNSSDATA &Obs, double *Azel, GNSSInfo &PInfo)
{
	int i;
	for(i = 0; i < Obs.NumSats; i++)
	{
		double C1 = G_Fre[1][0] * G_Fre[1][0] / (G_Fre[1][0] * G_Fre[1][0] - G_Fre[1][1] * G_Fre[1][1]);
		double C2 = 1 - C1;
		double C3 = G_Fre[1][0] * G_Fre[1][0] / (G_Fre[1][0] * G_Fre[1][0] - G_Fre[1][2] * G_Fre[1][2]);
		double C4 = 1 - C3;
		int prn_1 = Obs.ObsData[i].prn-1;
		if(Obs.ObsData[i].P[0] > PRECISION && Obs.ObsData[i].P[1] > PRECISION && Obs.ObsData[i].P[2] > PRECISION && Azel[i*2 +1] > 30 * D2R)
		{
			double Deta = (C1 - C3) * Obs.ObsData[i].P[0] + C2 * Obs.ObsData[i].P[1] - C4 * Obs.ObsData[i].P[2];
			PInfo.SSat[prn_1].DetaDCB = (PInfo.SSat[prn_1].DetaDCB * PInfo.SSat[prn_1].DCBNum + Deta)/ (PInfo.SSat[prn_1].DCBNum +1);
			PInfo.SSat[prn_1].DCBRMS = sqrt((PInfo.SSat[prn_1].DCBRMS * PInfo.SSat[prn_1].DCBRMS * PInfo.SSat[prn_1].DCBNum + Deta * Deta) / (PInfo.SSat[prn_1].DCBNum+1));
			PInfo.SSat[prn_1].DCBNum++;
			PInfo.MP_S[prn_1].MP_1[PInfo.MP_S[prn_1].LastEpoch[0]] = Deta;
			PInfo.MP_S[prn_1].MP_2[PInfo.MP_S[prn_1].LastEpoch[0]] = Azel[i*2 +1];
			PInfo.MP_S[prn_1].LastEpoch[0]++;
		}
	}
}

void GetAverageDCB(GNSSInfo &PInfo)
{
	int prn_1,i;
	for(prn_1 = 32; prn_1 < 32 + 14; prn_1++)
	{
		if(PInfo.MP_S[prn_1].LastEpoch[0] > 300)
		{
			double cigma = sqrt(PInfo.SSat[prn_1].DCBRMS * PInfo.SSat[prn_1].DCBRMS - PInfo.SSat[prn_1].DetaDCB * PInfo.SSat[prn_1].DetaDCB);
			double Ave = 0.0,Weight = 0.0,Temp;
			int Num = 0;
			for(i = 0; i < PInfo.MP_S[prn_1].LastEpoch[0]; i++)
			{
				if(fabs(PInfo.MP_S[prn_1].MP_1[i] - PInfo.SSat[prn_1].DetaDCB) < 2.5 * cigma)
				{
					if(PInfo.MP_S[prn_1].MP_2[i] < 40 * D2R)
					{
						Temp = 0.4;
					}
					else if(PInfo.MP_S[prn_1].MP_2[i] < 60 * D2R)
					{
						Temp = 0.6;
					}
					else
					{
						Temp = 1.0;
					}	
					Ave = (Ave * Weight + PInfo.MP_S[prn_1].MP_1[i] * Temp) /(Weight + Temp);
					Weight += Temp;
					Num++;
				}
			}
			PInfo.SSat[prn_1].DetaDCB = Ave;
			PInfo.SSat[prn_1].DCBNum  = Num;
		}
		else
		{
			PInfo.SSat[prn_1].DCBRMS  = 0.0;
			PInfo.SSat[prn_1].DetaDCB = 0.0;
			PInfo.SSat[prn_1].DCBNum  = 0;
		}
		
	}
}