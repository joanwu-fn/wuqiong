#include "../PreProcess/PreProcess.h"

void DetectError(const GNSSDATA Obs, GNSSInfo &PInfo)
{
	u2 i,prn_1;
	for(i = 0; i < Obs.NumSats; i++)
	{
		prn_1 = Obs.ObsData[i].prn-1;

		if(fabs(Obs.ObsData[i].P[0] - Obs.ObsData[i].P[1]) > 50.0)
		{
			PInfo.SSat[prn_1].Flag = true;
		}
		if(Obs.ObsData[i].P[2] > 1e-5)
		{
			if(fabs(Obs.ObsData[i].P[1] - Obs.ObsData[i].P[2]) > 50.0)
			{
				PInfo.SSat[prn_1].Flag = true;
			}
		}
	}
}
double GetAmbiguity(GNSSDATA Obs, u2 Numb, u2 Type)
{
	u2 i,SYS,Sys_Num;
	char StrSys;
	double Fre[3],P[3],L[3];
	
	
	Prn2Sys(Obs.ObsData[Numb].prn, SYS, Sys_Num, StrSys);

	for(i = 0; i < 3; i++)
	{
		Fre[i] = G_Fre[Sys_Num][i];
		P[i]   = Obs.ObsData[Numb].P[i];
		L[i]   = Obs.ObsData[Numb].L[i];
	}
	if(0 == Type) //N1   1,0,0
	{
		if(fabs(P[0]) < 1e-5 || fabs(L[0]) < 1e-5 )
		{
			return 0.0;
		}
		if(fabs(P[1]) < 1e-5)
		    return (P[0] / G_Lam[Sys_Num][0] - L[0]);
		else
		{
			double b = 2 * (G_Fre[Sys_Num][1]*G_Fre[Sys_Num][1]) / (G_Fre[Sys_Num][1]*G_Fre[Sys_Num][1] - G_Fre[Sys_Num][0]*G_Fre[Sys_Num][0]);
			double a = 1-b;
			return (a * P[0] + b*P[1])/ G_Lam[Sys_Num][0] - L[0];
		}
	}
	else if(1 == Type) //WL   1,-1,0
	{
		if(fabs(P[0]) < 1e-5 || fabs(P[1]) < 1e-5 || fabs(L[0]) < 1e-5 || fabs(L[1]) < 1e-5)
		{
			return 0.0;
		}
		return ((Fre[0] * P[0] + Fre[1] * P[1])*(Fre[0] - Fre[1])/((Fre[0]+Fre[1])*CLIGHT)-(L[0] - L[1]));
	}
	else if(2 == Type) //EWL  GPS:0,1,-1     BDS 0,-1,1
	{
		if(fabs(P[1]) < 1e-5 || fabs(P[2]) < 1e-5 || fabs(L[1]) < 1e-5 || fabs(L[2]) < 1e-5)
		{
			return 0.0;
		}
		if(SYS == SYS_BDS)
		{
			return ((Fre[1] * P[1] + Fre[2] * P[2])*(-Fre[1] + Fre[2])/((Fre[1]+Fre[2])*CLIGHT)-(-L[1] + L[2]));
		}
		else if(SYS == SYS_GPS)
		{
			return ((Fre[1] * P[1] + Fre[2] * P[2])*(Fre[1] - Fre[2])/((Fre[1] + Fre[2])*CLIGHT)-(L[1] - L[2]));
		}
		else
		{
			return 0.0;
		}	
	}
	else if(3 == Type)  // 1 2 Iono Free Ambiguity
	{
		double C1,C2;
		C1 =  Fre[0]*Fre[0]/(Fre[0]*Fre[0] - Fre[1]*Fre[1]);
		C2 = -Fre[1]*Fre[1]/(Fre[0]*Fre[0] - Fre[1]*Fre[1]);
		if(fabs(P[0]) < 1e-5 || fabs(P[1]) < 1e-5 || fabs(L[0]) < 1e-5 || fabs(L[1]) < 1e-5)
		{
			return 0.0;
		}
		return ((C1*P[0] +  C2*P[1]) -(C1 * L[0]* CLIGHT/Fre[0] + C2 * L[1] * CLIGHT/Fre[1]));
	}
	else if(4 == Type)  // 1 3 Iono Free Ambiguity
	{
		double C1,C2;
		C1 =  Fre[0]*Fre[0]/(Fre[0]*Fre[0] - Fre[2]*Fre[2]);
		C2 = -Fre[2]*Fre[2]/(Fre[0]*Fre[0] - Fre[2]*Fre[2]);
		if(fabs(P[0]) < 1e-5 || fabs(P[2]) < 1e-5 || fabs(L[0]) < 1e-5 || fabs(L[2]) < 1e-5)
		{
			return 0.0;
		}
		return ((C1*P[0] +  C2*P[2]) -(C1 * L[0]* CLIGHT/Fre[0] + C2 * L[2] * CLIGHT/Fre[2]));

	}
	return 0.0;
}
void InitialPara(GNSSInfo &PInfo, double Para, double Var, u2 Numb)
{
	u2 i;

	PInfo.m_X[Numb] = Para;

	for(i = 0; i < PInfo.Nx; i++)
	{
		 PInfo.m_P[Numb+i*PInfo.Nx] = PInfo.m_P[i+Numb*PInfo.Nx] = i == Numb? Var:0.0;
	}
	return;
}

void IniAmbiguity(const GNSSDATA Obs, GNSSInfo &PInfo)
{
	u2 i,j,prn_1,Fre=0,Numb;
    double Ambi[MAXCHN];
	double VarAm[3] = {900.0, 25.0, 1.0};
	if(0)
	{
		for(Fre = 0; Fre < NEFREQ; Fre++)
		{
			for(i = 0; i < Obs.NumSats; i++)
			{
				prn_1 = Obs.ObsData[i].prn - 1;

				if(PInfo.SSat[prn_1].CSlip[Fre] == true)
				{
					PInfo.m_X[IB(Obs.ObsData[i].prn, Fre)] = 0.0; //?
				}
				Ambi[i] = GetAmbiguity(Obs, i, Fre);
			}

			if(Fre == 2)
				continue;
			for(i = 0; i < Obs.NumSats; i++)
			{
				Numb = IB(Obs.ObsData[i].prn, Fre);
				if(fabs(PInfo.m_X[Numb]) < 1e-5 && fabs(Ambi[i]) > 1e-5)
				{
					InitialPara(PInfo, Ambi[i], VarAm[Fre], Numb);
				}
				else
				{
					continue;
				}
			}
		}
	}
	else
	{
		for(Fre = 0; Fre < 2; Fre++)
		{
			for(i = 0; i < Obs.NumSats; i++)
			{
				prn_1 = Obs.ObsData[i].prn - 1;

				if(0 == Fre)
				{
					if(PInfo.SSat[prn_1].CSlip[0] == true || PInfo.SSat[prn_1].CSlip[1] == true)
					{
						PInfo.m_X[IB(Obs.ObsData[i].prn, Fre)] = 0.0; //?
					}
				}
				else if(1 == Fre)
				{
					if(PInfo.SSat[prn_1].CSlip[0] == true || PInfo.SSat[prn_1].CSlip[2] == true)
					{
						PInfo.m_X[IB(Obs.ObsData[i].prn, Fre)] = 0.0; //?
					}
				}
				
				Ambi[i] = GetAmbiguity(Obs, i, Fre + 3);
			}

			for(i = 0; i < Obs.NumSats; i++)
			{
				Numb = IB(Obs.ObsData[i].prn, Fre);
				if(fabs(PInfo.m_X[Numb]) < 1e-5 && fabs(Ambi[i]) > 1e-5)
				{
					InitialPara(PInfo, Ambi[i], VarAm[Fre], Numb);
				}
				else
				{
					continue;
				}
			}
		}

	}
	
}
void UpdateTrop(GNSSInfo &PInfo, gtime_t NowObsTime)
{
	int m;
	m = PInfo.Mode > PMODE_PPP_FIXHOL? 2:1;
	for(int i = 0; i < m; i++)
	{
		u2 Numb = IT(i);
		double DetaT;

		if(fabs(PInfo.m_X[Numb]) < 1E-5)
		{
			InitialPara(PInfo, 0.15, 0.4, Numb);
		}
		else
		{
			DetaT = timediff(NowObsTime, PInfo.LastObsTime);
			PInfo.m_P[Numb + Numb * PInfo.Nx] += DetaT * 1e-3 * 1e-3;
		}
	}
	return;
}

void UpdateClk_Pos(GNSSInfo &PInfo,double X[])
{
	u2 i;

	for(i = 0; i < MAXSYS * 2; i++)
	{
		if(i < 4)
		{
			InitialPara(PInfo, X[i+3], 10000.0, i+3);
		}
		else
		{
			InitialPara(PInfo, X[i-1], 10000.0, i+3);
		}
	}

// 	for(i = 0; i < 3; i++)
// 	{
// 		InitialPara(PInfo, PInfo.Coor[i], 1e-8, i);
// 	}
// 	return;

	if(norm(PInfo.m_X, 3) < 1e-5)
	{
		for(i = 0; i < 3; i++)
		{
			InitialPara(PInfo, X[i], 900.0, i);
		}
	}
	if(PInfo.Mode == PMODE_PPP_KINEMA || PInfo.Mode == PMODE_BL_KINEMA)
	{
		for(i = 0; i < 3; i++)
		{
			InitialPara(PInfo, X[i], 900.0, i);
		}
	}
	
	return;
}
void UpdateISB(GNSSInfo &PInfo, gtime_t NowObsTime)
{
	double DetaT;
	if(fabs(PInfo.m_X[6+3]) < PRECISION)
	{
		InitialPara(PInfo, 0.05, 4.0, 9);
	}
	else
	{
		DetaT = timediff(NowObsTime, PInfo.LastObsTime);
		PInfo.m_P[9 + 9 * PInfo.Nx] += DetaT * 1e-4 * 1e-4;
	}
// 	if(fabs(PInfo.m_X[7+3]) < PRECISION)
// 	{
// 		InitialPara(PInfo, 0.05, 4.0, 10);
// 	}
// 	else
// 	{
// 		DetaT = timediff(NowObsTime, PInfo.LastObsTime);
// 		PInfo.m_P[10 + 10 * PInfo.Nx] += DetaT * 1e-4 * 1e-4;
// 	}
}
/****************************************************************
 *Function  Name:          PPPPreProcess
 *Input Parameter:  GNSSDATA  Obs  观测值
 *Input Parameter:  PPPInfo  PInfo PPP信息
 *Function:         对数据进行预处理（粗差，周跳，模糊度初值）
 ***************************************************************/
void PreProcess(const GNSSDATA Obs, GNSSInfo &PInfo,double X[])
{
	DetectError(Obs, PInfo);

	UpdateClk_Pos(PInfo, X);

	//UpdateISB(PInfo, Obs.ObsTime);

	UpdateTrop(PInfo, Obs.ObsTime);

	DetectCS_GF2(Obs, PInfo.SSat);
	DetectCS_GF2_13(Obs, PInfo.SSat);

	DetectCS_MW12(Obs, PInfo.SSat);
	DetectCS_MW23(Obs, PInfo.SSat);
	
	IniAmbiguity(Obs, PInfo);
	return;
}