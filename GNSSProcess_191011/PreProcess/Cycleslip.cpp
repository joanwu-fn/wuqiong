#include "../PreProcess/PreProcess.h"
#include "../Common/GNSSERROR.h"

double GFObs_L1L2(const GNSSDATA Obs, u2 Numb)
{
	double L1,L2;

	L1 = Obs.ObsData[Numb].L[0];
	L2 = Obs.ObsData[Numb].L[1];

	if(fabs(L1) < 1e-5 || fabs(L2) < 1e-5)
	{
		return 0.0;
	}
	else
	{
		u2 SYS,Sys_Num;
		char StrSys;
		u2 prn = Prn2Sys(Obs.ObsData[Numb].prn, SYS, Sys_Num, StrSys);
		if(Sys_Num == 1 && prn > 16)
		{
			double Lambda2 = CLIGHT / FREB5I;
			return L1 * G_Lam[Sys_Num][0] - L2 * Lambda2;
		}
		else
		{
			return L1 * G_Lam[Sys_Num][0] - L2 * G_Lam[Sys_Num][1];
		}
	}
}

double GFObs_L1L3(const GNSSDATA Obs, u2 Numb)
{
	double L1,L3;

	L1 = Obs.ObsData[Numb].L[0];
	L3 = Obs.ObsData[Numb].L[2];

	if(fabs(L1) < 1e-5 || fabs(L3) < 1e-5)
	{
		return 0.0;
	}
	else
	{
		u2 SYS,Sys_Num;
		char StrSys;
		Prn2Sys(Obs.ObsData[Numb].prn, SYS, Sys_Num, StrSys);
		return L1 * G_Lam[Sys_Num][0] - L3 * G_Lam[Sys_Num][2];
	}
}
void DetectCS_MW12(const GNSSDATA Obs, SSat_t SSat[MAXSAT])
{
	int i,prn_1;
	double MW;
	for(i = 0; i < Obs.NumSats; i++)
	{
		prn_1 = Obs.ObsData[i].prn-1;
		if(fabs(MW  = MWObs(Obs, i, 90, 1, 2)) < 1e-5)
		{
			continue;
		}
		if(SSat[prn_1].CSMW[2] > 0)
		{
			if(fabs(MW - SSat[prn_1].CSMW[0]) > 4 * SSat[prn_1].CSMW[1])
			{
				SSat[prn_1].CSlip[0] = SSat[prn_1].CSlip[1] = true;
				SSat[prn_1].CSMW[0] = MW;
				SSat[prn_1].CSMW[1] = 0.5;
				SSat[prn_1].CSMW[2] = 1;
			}
			else
			{
				double lmw = SSat[prn_1].CSMW[0];
				double var = SSat[prn_1].CSMW[1];
				double numb = SSat[prn_1].CSMW[2];
				SSat[prn_1].CSMW[0] = (lmw * numb + MW) / (numb + 1);
				SSat[prn_1].CSMW[1] = sqrt((var * var * numb + (MW - lmw) * (MW - lmw)) / (numb + 1));
				SSat[prn_1].CSMW[2] += 1;
			}
		}
		else
		{
			SSat[prn_1].CSMW[0] = MW;
			SSat[prn_1].CSMW[1] = 0.5;
			SSat[prn_1].CSMW[2] = 1;
		}	
	}
	return;
}

void DetectCS_MW23(const GNSSDATA Obs, SSat_t SSat[MAXSAT])
{
	int i,prn_1;
	double MW;
	for(i = 0; i < Obs.NumSats; i++)
	{
		prn_1 = Obs.ObsData[i].prn-1;
		if(prn_1 < NSATGPS)
		{
			if(fabs(MW  = MWObs(Obs, i, 90, 2, 3)) < 1e-5)
			{
				continue;
			}
		}
		else  //need modification if use more system
		{
			if(fabs(MW  = MWObs(Obs, i, 90, 3, 2)) < 1e-5)
			{
				continue;
			}
		}
		
		if(SSat[prn_1].CSMW[5] > 0)
		{
			if(fabs(MW - SSat[prn_1].CSMW[3]) > 4 * SSat[prn_1].CSMW[4])
			{
				SSat[prn_1].CSlip[0] = SSat[prn_1].CSlip[1] = true;
				SSat[prn_1].CSMW[3] = MW;
				SSat[prn_1].CSMW[4] = 0.5;
				SSat[prn_1].CSMW[5] = 1;
			}
			else
			{
				double lmw = SSat[prn_1].CSMW[3];
				double var = SSat[prn_1].CSMW[4];
				double numb = SSat[prn_1].CSMW[5];
				SSat[prn_1].CSMW[3] = (lmw * numb + MW) / (numb + 1);
				SSat[prn_1].CSMW[4] = sqrt((var * var * numb + (MW - lmw) * (MW - lmw)) / (numb + 1));
				SSat[prn_1].CSMW[5] += 1;
			}
		}
		else
		{
			SSat[prn_1].CSMW[3] = MW;
			SSat[prn_1].CSMW[4] = 0.5;
			SSat[prn_1].CSMW[5] = 1;
		}	
	}
	return;
}
void DetectCS_GF2(const GNSSDATA Obs, SSat_t SSat[MAXSAT])
{
	double GfNow,GfPre1,GfPre2,DetaGf;
	u2 i,prn_1;
	double Gf[MAXSAT] = {0};
	for(i=0; i < MAXSAT; i++)
	{
		SSat[i].Flag = false;
		SSat[i].CSlip[0] = SSat[i].CSlip[1] = SSat[i].CSlip[2] = false;
	}
	for(i = 0; i < Obs.NumSats; i++)
	{
		if(Obs.ObsData[i].prn <= 0 || Obs.ObsData[i].prn > MAXSAT)
		{
			continue;
		}
		prn_1 = Obs.ObsData[i].prn-1;
		if(fabs(GfNow  = GFObs_L1L2(Obs, i)) < 1e-5)
		{
			continue;
		}
		Gf[prn_1] = GfNow;
		GfPre1 = SSat[prn_1].Gf[0];
		GfPre2 = SSat[prn_1].Gf[1];
		SSat[prn_1].Gf[0] = GfPre2;
		SSat[prn_1].Gf[1] = GfNow;

		if(fabs(GfPre1) < 1e-5 && fabs(GfPre2) < 1e-5)  //The First Epoch
		{
			SSat[prn_1].TrackNum = 0;
			continue;
		}
		else if(fabs(GfPre1) < 1e-5)  //The Second Epoch
		{
			SSat[prn_1].TrackNum += 1;
			if (fabs(GfPre2 - GfNow) > 0.02) 
			{
				SSat[prn_1].CSlip[0] = true;
				SSat[prn_1].CSlip[1] = true;
			}
		}
		else
		{
			DetaGf = GfNow - 2*GfPre2 + GfPre1;
			if(fabs(DetaGf) < 0.05)  
			{
				if(fabs(GfPre2 - GfPre1) < 5.0)    //没有周跳
				{
					SSat[prn_1].TrackNum += 1;
					continue;
				}
				else
				{
					SSat[prn_1].CSlip[0] = true;
					SSat[prn_1].CSlip[1] = true;
				}
			}
			else
			{
				if(fabs(GfNow - GfPre2) < 0.025)    //前一历元周跳
				{   
					SSat[prn_1].TrackNum += 1;
					continue;
				}
				else                               //当前历元周跳
				{
					SSat[prn_1].CSlip[0] = true;
					SSat[prn_1].CSlip[1] = true;
				}
			}		
		}
		//for i : 0 to Number
	}	
}


void DetectCS_GF2_13(const GNSSDATA Obs, SSat_t SSat[MAXSAT])
{
	double GfNow,GfPre1,GfPre2,DetaGf;
	u2 i,prn_1;
	double Gf[MAXSAT] = {0};

	for(i = 0; i < Obs.NumSats; i++)
	{
		prn_1 = Obs.ObsData[i].prn-1;
		if(fabs(GfNow  = GFObs_L1L3(Obs, i)) < 1e-5)
		{
			continue;
		}
		Gf[prn_1] = GfNow;
		GfPre1 = SSat[prn_1].Gf[2];
		GfPre2 = SSat[prn_1].Gf[3];
		SSat[prn_1].Gf[2] = GfPre2;
		SSat[prn_1].Gf[3] = GfNow;

		if(fabs(GfPre1) < 1e-5 && fabs(GfPre2) < 1e-5)  //The First Epoch
		{
			continue;
		}
		else if(fabs(GfPre1) < 1e-5)  //The Second Epoch
		{
			if (fabs(GfPre2 - GfNow) > 0.02) {
				SSat[prn_1].CSlip[0] = true;
				SSat[prn_1].CSlip[2] = true;
			}
		}
		else
		{
			DetaGf = GfNow - 2*GfPre2 + GfPre1;
			if(fabs(DetaGf) < 0.05)  
			{
				if(fabs(GfPre2 - GfPre1) < 5.0)    //没有周跳
				{
					continue;
				}
				else
				{
					SSat[prn_1].CSlip[0] = true;
					SSat[prn_1].CSlip[2] = true;
				}
			}
			else
			{
				if(fabs(GfNow - GfPre2) < 0.025)    //前一历元周跳
				{   
					continue;
				}
				else                               //当前历元周跳
				{
					SSat[prn_1].CSlip[0] = true;
					SSat[prn_1].CSlip[2] = true;
				}
			}		
		}
		//for i : 0 to Number
	}	
}