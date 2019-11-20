#include"../GNSSProcess/PPPProcess.h"
double P1C1[32]={
 1.552,-1.037,1.774,0.707,1.631,1.960,1.095,-1.140,0.319,
-2.415,-0.043, 1.048, 0.921,-0.369,1.497,-1.029,1.142,-0.745,
-2.149,-1.810,-1.148,-2.516,0.149, 1.535,-0.348,-1.104, 0.098,
-0.678,1.608,-0.104, 1.078,-1.478};
double GetMeasureMent(GNSSDATA Obs, u2 Numb, u2 SysNum, u2 Type, u2 Fre)
{
	double P[3],L[3];
	int i;
	double C1,C2,C3,C4;

	C1 =  SQR(G_Lam[SysNum][1])/(SQR(G_Lam[SysNum][1])-SQR(G_Lam[SysNum][0]));
	C2 = -SQR(G_Lam[SysNum][0])/(SQR(G_Lam[SysNum][1])-SQR(G_Lam[SysNum][0]));
	C3 =  SQR(G_Lam[SysNum][2])/(SQR(G_Lam[SysNum][2])-SQR(G_Lam[SysNum][0]));
	C4 = -SQR(G_Lam[SysNum][0])/(SQR(G_Lam[SysNum][2])-SQR(G_Lam[SysNum][0]));

	for(i = 0; i < 3; i++)
	{
		P[i] = Obs.ObsData[Numb].P[i];
		L[i] = Obs.ObsData[Numb].L[i] * G_Lam[SysNum][i];
	}
	if(Obs.ObsData[Numb].prn < 33)
	{
		P[0] -= P1C1[Obs.ObsData[Numb].prn-1]*0.3;
	}

	if(1)
	{
		if(Type == 1) // pusedo-range
		{
			if(Fre == 0) //L1/B1  L2/B2 ionoshpere-free
			{
				if(fabs(P[0]) < 1e-5 || fabs(P[1]) < 1e-5)
				{
					return 0.0;
				}
				return (C1*P[0]  + C2*P[1]);
			}
			else       //L1/B1  L3/B3 ionoshpere-free
			{
				if(fabs(P[0]) < 1e-5 || fabs(P[2]) < 1e-5)
				{
					return 0.0;
				}
				return (C3*P[0]  + C4*P[2]);		
			}

		}
		else         //carrier phase
		{
			if(Fre == 0) //L1/B1  L2/B2 ionoshpere-free
			{
				if(fabs(L[0]) < 1e-5 || fabs(L[1]) < 1e-5)
				{
					return 0.0;
				}
				return C1*L[0]  + C2*L[1];
			}
			else       //L1/B1  L3/B3 ionoshpere-free
			{
				if(fabs(L[0]) < 1e-5 || fabs(L[2]) < 1e-5)
				{
					return 0.0;
				}
				return C3*L[0]  + C4*L[2];

			}
		}
	}
}
bool GetUDSatStaInfo(GNSSInfo GNSSInf, const GNSSDATA Obs, int Number, const double SatXYZ[], double *rb, double rr[3], UDSatStaInfo &Info)
{
	double pos[3];
	char StrSys;
	int i;
	Prn2Sys(Obs.ObsData[Number].prn, Info.SYS, Info.SysNum, StrSys);
	ecef2pos(rr, pos);
	for(i = 0; i < NEFREQ*2; i++)
	{
		Info.OMC[i] = INVALID_VALUE;
	}
	if ((Info.Geometry = geodist(&SatXYZ[Number*3], rr, Info.e)) <= 1e-5 )
		return false;
	satazel(pos,Info.e, &Info.Azel[0]);
	GetTropDelay(Obs.ObsTime, pos, Info.Azel, Info.Dtrop, Info.Vtrop);
	Info.Dtrop += prectrop(Obs.ObsTime, pos, &Info.Azel[0], &GNSSInf.m_X[IT(0)], Info.TropMFuncR, Info.Vtrop);

	if(NULL != rb)  //the observation is SD between station
	{
		double Geometry,e[3],Dtrop,TropMFunc[3];
		ecef2pos(rb, pos);
		if ((Geometry = geodist(&SatXYZ[(Number+Obs.NumSats)*3], rb, e)) <= 1e-5 )
			return false;
		satazel(pos, e, &Info.Azel[2]);
		Info.Geometry -= Geometry;
		GetTropDelay(Obs.ObsTime, pos, &Info.Azel[2], Dtrop, Info.Vtrop);
		Dtrop += prectrop(Obs.ObsTime, pos, &Info.Azel[2], &GNSSInf.m_X[IT(1)], Info.TropMFuncB, Info.Vtrop);
		Info.Dtrop -=Dtrop;
	}
	int Nf = GNSSInf.Nf;
	for(i = 0; i < Nf*2; i++)
	{
		Info.Observation[i] = GetMeasureMent(Obs, Number,  Info.SysNum, i/Nf, i%Nf);
		if(fabs(Info.Observation[i]) < 1e-5)
			continue;
		Info.OMC[i] = Info.Observation[i] - (Info.Geometry + Info.Dtrop +Info.Diono);
		Info.Weight[i] = varerr(Info.Azel[1], i/Nf, Info.SYS);
		if(i < Nf)
		{
			int Apos = IB(Obs.ObsData[Number].prn, i%Nf);
			Info.OMC[i] += GNSSInf.m_X[Apos];
		}
	}
	return true;	
}
int Get_UnDiff_Res(u2 Iter, GNSSInfo &PInfo, const GNSSDATA Obs, const double SatXYZ[], const double SatClk[], 
				   double Azel[], double EphRMS[], double H[], double V[], double R[])
{
	u2 i,j,f,Nv = 0,Nx = PInfo.Nx,Apos[2];
	PInfo.Nf = 2;
	int Nf = PInfo.Nf;
	double *Var,Disp[3];
	double rb[3],rr[3],pos[3];
	Var = zeros(Obs.NumSats * 2 * 3, 1);
	/*------------------------------------Earth Tide Correction : ROVER-----------------------*/
	for (i=0;i<3;i++) 
	{
		rb[i] = PInfo.Station[PInfo.N_Base].XYZ[i];
		rr[i] = PInfo.m_X[i];
	}
	ecef2pos(rr, pos);
	tidedisp(gpst2utc(Obs.ObsTime), rr, 1, NULL, NULL, Disp);
	for (i=0;i<3;i++) 
	{
		rr[i]+=Disp[i];
	}
	ecef2pos(rb, pos);
	tidedisp(gpst2utc(Obs.ObsTime), rb, 1, NULL, NULL, Disp);
	for (i=0;i<3;i++) 
	{
		rb[i]+=Disp[i];
	}
	/*------------------------------------Earth Tide Correction : ROVER-----------------------*/
	UDSatStaInfo *UNSatInfo;
	UNSatInfo = NULL;
	UNSatInfo = (UDSatStaInfo *)realloc(UNSatInfo, sizeof(UDSatStaInfo)*Obs.NumSats);
	memset(UNSatInfo, 0x00, sizeof(UDSatStaInfo)*Obs.NumSats);
	for(i = 0; i < Obs.NumSats; i++)
	{
		GetUDSatStaInfo(PInfo, Obs, i, SatXYZ, rb, rr, UNSatInfo[i]);
	}
	for(f = 0; f < Nf * 2 + 1; f++)
	{
		int Ref = -1;
		if(f < Nf * 2)
		{
			for(i = 0; i < Obs.NumSats; i++)
			{
				if(fabs(UNSatInfo[i].OMC[f]-INVALID_VALUE) > PRECISION && UNSatInfo[i].SYS == SYS_BDS)
				{
					if(-1 == Ref || UNSatInfo[i].Azel[1] > UNSatInfo[Ref].Azel[1])
					{
						Ref = i;
					}
				}
			}
			if(-1 == Ref)
				continue;
			for(i = 0; i < Obs.NumSats; i++)
			{
				if(fabs(UNSatInfo[i].OMC[f]-INVALID_VALUE) < PRECISION || i == Ref)
				{
					continue;
				}
				for(j = 0; j < 3; j++)
				{
					H[j + Nx * Nv] = -(UNSatInfo[i].e[j] - UNSatInfo[Ref].e[j]);
				}
				H[IT(0)+ Nx * Nv] =   (UNSatInfo[i].TropMFuncR[0] - UNSatInfo[Ref].TropMFuncR[0]);
				H[IT(1)+ Nx * Nv] =  -(UNSatInfo[i].TropMFuncB[0] - UNSatInfo[Ref].TropMFuncB[0]);
				V[Nv] = UNSatInfo[i].OMC[f] - UNSatInfo[Ref].OMC[f];
				Var[Nv] = UNSatInfo[i].Weight[f] + UNSatInfo[Ref].Weight[f];
// 				if(SYS_BDS == UNSatInfo[i].SYS)
// 				{
// 					if(f % Nf == 0)
// 					{
// 						H[6 + 3  + Nv * Nx] = 1.0;
// 						V[Nv] +=PInfo.m_X[6 + 3];
// 					}
// 				}
				if(f < Nf)
				{
					Apos[0] = IB(Obs.ObsData[i].prn, f);
					Apos[1] = IB(Obs.ObsData[Ref].prn, f);
					H[Apos[0] + Nv * Nx] = -1.0;
					H[Apos[1] + Nv * Nx] =  1.0;
				}
				Nv++;
			}
		}	
	}
	return Nv;
}
/*******************************************************/
int  Get_SingleDiff_Res(u2 Iter, GNSSInfo &PInfo, const GNSSDATA Obs, const double SatXYZ[], const double SatClk[], 
				 double Azel[], double EphRMS[], double H[], double V[], double R[])
{
	u2 i,j,f,Nv = 0,Nx = PInfo.Nx,Apos[2];
	PInfo.Nf = 2;
	int Nf = PInfo.Nf;
	double *Var,Disp[3];
	double rb[3],rr[3],pos[3];
	Var = zeros(Obs.NumSats * 2 * 3, 1);
	/*------------------------------------Earth Tide Correction : ROVER-----------------------*/
	for (i=0;i<3;i++) 
	{
		rb[i] = PInfo.Station[PInfo.N_Base].XYZ[i];
		rr[i] = PInfo.m_X[i];
	}
	ecef2pos(rr, pos);
	tidedisp(gpst2utc(Obs.ObsTime), rr, 1, NULL, NULL, Disp);
	for (i=0;i<3;i++) 
	{
		rr[i]+=Disp[i];
	}
	ecef2pos(rb, pos);
	tidedisp(gpst2utc(Obs.ObsTime), rb, 1, NULL, NULL, Disp);
	for (i=0;i<3;i++) 
	{
		rb[i]+=Disp[i];
	}
	/*------------------------------------Earth Tide Correction : ROVER-----------------------*/
	UDSatStaInfo *UNSatInfo;
	UNSatInfo = NULL;
	UNSatInfo = (UDSatStaInfo *)realloc(UNSatInfo, sizeof(UDSatStaInfo)*Obs.NumSats);
	memset(UNSatInfo, 0x00, sizeof(UDSatStaInfo)*Obs.NumSats);
	for(i = 0; i < Obs.NumSats; i++)
	{
		GetUDSatStaInfo(PInfo, Obs, i, SatXYZ, rb, rr, UNSatInfo[i]);
	}
	for(int m = 0; m < 2; m++)
	{
		for(f = 0; f < Nf * 2 + 1; f++)
		{
			int Ref = -1;
			if(f < Nf * 2)
			{
				for(i = 0; i < Obs.NumSats; i++)
				{
					if(UNSatInfo[i].SysNum == m && fabs(UNSatInfo[i].OMC[f]-INVALID_VALUE) > PRECISION)
					{
						if(-1 == Ref || UNSatInfo[i].Azel[1] > UNSatInfo[Ref].Azel[1])
						{
							Ref = i;
						}
					}
				}
				if(-1 == Ref)
					continue;
				for(i = 0; i < Obs.NumSats; i++)
				{
					if(fabs(UNSatInfo[i].OMC[f]-INVALID_VALUE) < PRECISION || UNSatInfo[i].SysNum != m || i == Ref)
					{
						continue;
					}
					for(j = 0; j < 3; j++)
					{
						H[j + Nx * Nv] = -(UNSatInfo[i].e[j] - UNSatInfo[Ref].e[j]);
					}
					H[IT(0)+ Nx * Nv] =   (UNSatInfo[i].TropMFuncR[0] - UNSatInfo[Ref].TropMFuncR[0]);
					H[IT(1)+ Nx * Nv] =  -(UNSatInfo[i].TropMFuncB[0] - UNSatInfo[Ref].TropMFuncB[0]);
					V[Nv] = UNSatInfo[i].OMC[f] - UNSatInfo[Ref].OMC[f];
					Var[Nv] = UNSatInfo[i].Weight[f] + UNSatInfo[Ref].Weight[f];
					if(f < Nf)
					{
						Apos[0] = IB(Obs.ObsData[i].prn, f);
						Apos[1] = IB(Obs.ObsData[Ref].prn, f);
						H[Apos[0] + Nv * Nx] = -1.0;
						H[Apos[1] + Nv * Nx] =  1.0;
					}
					Nv++;
				}
			}	
		}
	}
	for(i = 0; i < Nv; i++)
	{
		for(j = 0; j < Nv; j++)
		{
			R[i+j*Nv]=i==j?Var[i]:0.0;
		}
	}
	free(Var);    free(UNSatInfo);
	return Nv;
}