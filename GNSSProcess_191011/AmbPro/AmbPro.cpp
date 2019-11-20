#include "../Common/AmbPro.h"
#include "../Common/GNSSERROR.h"
#define ROUND(x)        (int)floor((x)+0.5)
#define FIX_THRES       0.129       /* fix threshold (cycle): p0=0.9999 */
const double EWL_Fra[14] = {0.01, 0.40, 0.50, 0.39, 0.57, 0.99, 0.44, 0.26, 0.40, 0.13, 0.80, 0.08, 0.00, 0.55};
const double  WL_Fra[14] = {0.01, 0.72, 0.50, 0.98, 0.75, 0.26, 0.85, 0.70, 0.75, 0.02, 0.13, 0.98, 0.00, 0.67};


double Guass(double x)
{
	if(x < -5.0)
	{
		return 0.0;
	}
	else if(x > 5.0)
	{
		return 1.0;
	}
	else
	{
		long double w,y=0,z,sum,b=1.0;
		int n=1, a=1, m = 1;
		while(n < 101)
		{
			w = pow(x,n);
			z = a * w /(b*n);
			y+=z;
			a*=-1;
			b*=2*m;
			m++;
			n+=2;
		}
		sum = 0.5 + y/sqrt(2*PI);
		return sum;
	}	
}
double probability(double ave, double sig)
{
	double deta = ave,fai1,fai2;
	double gua1,gua2;
	while(fabs(deta) > 0.5)
	{
		if(deta < 0.0)
		{
			deta += 1.0;
		}
		else
		{
			deta -= 1.0;
		}
	}
	fai1 = ( 0.5 - deta) / sig;
	fai2 = (-0.5 - deta) / sig;
	gua1 = Guass(fai1);
	gua2 = Guass(fai2);
	return gua1 - gua2;
}
void GetAverageAmb(const GNSSDATA &Obs, SSat_t *Ssat, double Azel[],char *Station)
{
	int i,j;
	u2 SYS,SysNum;
	double Amb[2],Var[2];
	char StrSys;
	double Pse[3],Carr[3];
	double lam_WL[2] = {CLIGHT/(FREQ1 - FREQ2), CLIGHT/(FREB1 - FREB2)},lam_EWL[2] = {CLIGHT/(FREQ2 - FREQ5), CLIGHT/(FREB5 - FREB2)};

	for(i = 0; i < Obs.NumSats; i++)
	{
		if(fabs(Azel[1 + 2 * i]) < 1e-5)
		{
			continue;
		}

		Prn2Sys(Obs.ObsData[i].prn, SYS, SysNum, StrSys);
		/*                 EWL                      */
		if(SYS == SYS_BDS)
		{
			Amb[0] = MWObs(Obs, i, Azel[i*2+1], 3, 2); // 0, -1, 1
		}
		else if(SYS == SYS_GPS)
		{
			Amb[0] = MWObs(Obs, i, Azel[i*2+1], 2, 3); // 0, 1, -1
		}

		Var[0]=SD_var(var_LC(0, 1, 1, 0.3, SYS),Azel[1 + 2 * i]) / (lam_EWL[SysNum] * lam_EWL[SysNum]); 
		/*                 WL                      */
		if(SYS == SYS_BDS)
		{
			Pse[0] = 0.4889/lam_WL[SysNum];  Pse[1] = 0.1882/lam_WL[SysNum];   Pse[2] = 0.3229/lam_WL[SysNum];
			Carr[0] = 1.0;    Carr[1] = -1.0;    Carr[2] = 0.0;
			Amb[1] = GetCombination(Obs, i, Azel[i*2+1], Pse, Carr);
			if(fabs(Amb[0]) < 1e-5 || fabs(Amb[1]) < 1e-5)
			{
				Ssat[Obs.ObsData[i].prn - 1].Number[0] =  0;
				Ssat[Obs.ObsData[i].prn - 1].Number[0] =  0;
				continue;
			}
			Amb[1] = Amb[1] - 0.4631*Amb[0]*lam_EWL[SysNum] / lam_WL[SysNum];
			Var[1]=SD_var(var_LC(0, 1, 1, 0.3, SYS),Azel[1 + 2 * i]) / (lam_WL[SysNum] * lam_WL[SysNum] * 2.25);       //除以1.5是因为超宽巷辅助了宽巷的固定
		}
		else
		{
			Amb[1] = MWObs(Obs, i, Azel[i*2+1], 1, 2); // 1, -1, 0
			if(fabs(Amb[1]) < 1e-5)
			{
				continue;
			}
			Var[1]=SD_var(var_LC(0, 1, 1, 0.3, SYS),Azel[1 + 2 * i]) / (lam_WL[SysNum] * lam_WL[SysNum]);       //除以1.5是因为超宽巷辅助了宽巷的固定
		}

		if(Ssat[Obs.ObsData[i].prn - 1].CSlip[0] == true || Ssat[Obs.ObsData[i].prn - 1].CSlip[1] == true || Ssat[Obs.ObsData[i].prn - 1].CSlip[2] == true)
		{
			if(fabs(Amb[0]) > 1e-5)
			{
				Ssat[Obs.ObsData[i].prn - 1].Mw[0] =  Amb[0];
				Ssat[Obs.ObsData[i].prn - 1].Var[0] = Var[0];
				Ssat[Obs.ObsData[i].prn - 1].Number[0] =  1;
			}

			if(fabs(Amb[1]) > 1e-5)
			{
				Ssat[Obs.ObsData[i].prn - 1].Mw[1] =  Amb[1];
				Ssat[Obs.ObsData[i].prn - 1].Var[1] = Var[1];
				Ssat[Obs.ObsData[i].prn - 1].Number[1] =  1;
			}
		}
		else
		{
			if(fabs(Amb[0]) > 1e-5)
			{
				Ssat[Obs.ObsData[i].prn - 1].Mw[0] +=  (Amb[0] - Ssat[Obs.ObsData[i].prn - 1].Mw[0])/(Ssat[Obs.ObsData[i].prn - 1].Number[0] +1);
				Ssat[Obs.ObsData[i].prn - 1].Var[0] +=  (Var[0] - Ssat[Obs.ObsData[i].prn - 1].Var[0])/(Ssat[Obs.ObsData[i].prn - 1].Number[0] +1);
				Ssat[Obs.ObsData[i].prn - 1].Number[0] +=1;
			}

			if(fabs(Amb[1]) > 1e-5)
			{
				Ssat[Obs.ObsData[i].prn - 1].Mw[1] +=  (Amb[1] - Ssat[Obs.ObsData[i].prn - 1].Mw[1])/(Ssat[Obs.ObsData[i].prn - 1].Number[1] +1);
				Ssat[Obs.ObsData[i].prn - 1].Var[1] +=  (Var[1] - Ssat[Obs.ObsData[i].prn - 1].Var[1])/(Ssat[Obs.ObsData[i].prn - 1].Number[1] +1);
				Ssat[Obs.ObsData[i].prn - 1].Number[1] +=1;
			}
		}
	}
}
void Fix_EWL_Amb(const GNSSDATA &Obs, SSat_t *Ssat, double Azel[],char *Station)
{
	int i,j,m,Ref = -1,number[2] = {0};
	double lam_WL[2] = {CLIGHT/(FREQ1 - FREQ2), CLIGHT/(FREB1 - FREB2)},lam_EWL[2] = {CLIGHT/(FREQ2 - FREQ5), CLIGHT/(FREB5 - FREB2)};

	for(i = 0; i < MAXSAT; i++)
	{
		Ssat[i].FixWL[0]   = -999;
		Ssat[i].FixFlag[0] = 0;
		Ssat[i].FixFlag[1] = 0;
	}
	GetAverageAmb(Obs, Ssat, Azel,Station);

	u2 Sys,SysNum;
	char StrSys;
	int TT=0;
	for(m = 1; m < 2; m++)
	{
		for(i = 0; i < Obs.NumSats; i++)
		{
			Prn2Sys(Obs.ObsData[i].prn, Sys, SysNum, StrSys);
			if((0 == m && Sys != SYS_GPS) || (1 == m && Sys != SYS_BDS))
				continue;
			if(-1 == Ref)
			{
				Ref = i;
			}
			else
			{
				if(Ssat[Obs.ObsData[i].prn - 1].Number[0] > Ssat[Obs.ObsData[Ref].prn - 1].Number[0])
				{
					Ref = i;
				}
			}	
		}
		if(Ref < 0)
		{
			continue;
		}
		double EWL,iEWL,Var;
		int fix;
		for(i = 0; i < Obs.NumSats; i++)
		{
			Prn2Sys(Obs.ObsData[i].prn, Sys, SysNum, StrSys);
			if((0 == m && Sys != SYS_GPS) || (1 == m && Sys != SYS_BDS) || i == Ref)
				continue;

			TT++;
			EWL = Ssat[Obs.ObsData[i].prn - 1].Mw[0] - Ssat[Obs.ObsData[Ref].prn - 1].Mw[0] - (EWL_Fra[Obs.ObsData[i].prn - 33] - EWL_Fra[Obs.ObsData[Ref].prn - 33]);
			Var = Ssat[Obs.ObsData[i].prn - 1].Var[0] / Ssat[Obs.ObsData[i].prn - 1].Number[0] + Ssat[Obs.ObsData[Ref].prn - 1].Var[0] / Ssat[Obs.ObsData[Ref].prn - 1].Number[0];
			iEWL = ROUND(EWL);
			fix=fabs(EWL - iEWL) <= FIX_THRES*2 && sqrt(Var) <= FIX_THRES;
			if(fix)
			{
				Ssat[Obs.ObsData[i].prn - 1].FixFlag[0]   = 1;
				Ssat[Obs.ObsData[Ref].prn - 1].FixFlag[0] = 1;
				Ssat[Obs.ObsData[i].prn - 1].FixWL[0]     = iEWL;
				Ssat[Obs.ObsData[Ref].prn - 1].FixWL[0]   = 0;
				
				EWL = Ssat[Obs.ObsData[i].prn - 1].Mw[1] - Ssat[Obs.ObsData[Ref].prn - 1].Mw[1] - (WL_Fra[Obs.ObsData[i].prn - 33] - WL_Fra[Obs.ObsData[Ref].prn - 33]);
				Var = Ssat[Obs.ObsData[i].prn - 1].Var[1] / Ssat[Obs.ObsData[i].prn - 1].Number[1] + Ssat[Obs.ObsData[Ref].prn - 1].Var[1] / Ssat[Obs.ObsData[Ref].prn - 1].Number[1];
				iEWL = ROUND(EWL);
				fix=fabs(EWL - iEWL) <= FIX_THRES*1.5 && sqrt(Var) <= FIX_THRES;
				number[0]++;
				if(fix && Ssat[Obs.ObsData[i].prn - 1].FixCount > 10)
				{
					Ssat[Obs.ObsData[i].prn - 1].FixWL[1] = iEWL;
					Ssat[Obs.ObsData[Ref].prn - 1].FixWL[1] = 0;
					Ssat[Obs.ObsData[i].prn - 1].FixFlag[1] = 1;
					Ssat[Obs.ObsData[Ref].prn - 1].FixFlag[1] = 1;
					number[1]++;
				}
				else if(fix)
				{
					Ssat[Obs.ObsData[i].prn - 1].FixCount++;
				}
			}
		}
	}//m=0:2
}