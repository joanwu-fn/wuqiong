#include "../GNSSProcess/PPPProcess.h"
#include "../Common/GNSSERROR.h"
#include "../Common/AmbPro.h"
#include "../PreProcess/PreProcess.h"

const double AmbCorr[14] = {0, -30.28, -14.37, 3.0, -10, 26.10, 49.26, 42.35, 53.20, 47.75, 37.70, 30.63, 0, 45.30};

bool GetReceiverFCB(double *FCB, int Numb,  double &AveFCB)
{
	AveFCB = FCB[0];
	int RecNum = 0;

	if(Numb < 2)
	{
		return false;
	}

	for(int i = 1; i < Numb; i++)
	{
		while(fabs(FCB[i] - FCB[0]) > 0.5)
		{
			if(FCB[i] - FCB[0] > 0.5)
			{
				FCB[i] -= 1.0;
			}
			else
			{
				FCB[i] += 1.0;
			}
		}
		AveFCB = (AveFCB * RecNum + FCB[i]) / (RecNum + 1);
		RecNum++;
	}
	return true;
}
PPP::PPP()
{
	_DCB13[0] = 7.8086;  _DCB13[1] = -5.5289; _DCB13[2] =  -2.8895; _DCB13[3] = -0.7787;
	_DCB13[4] = -6.1334; _DCB13[5] =  1.9368; _DCB13[6] = 8.0206; _DCB13[7] = 5.4828;
	_DCB13[8] = 1.1097;  _DCB13[9] = 1.0888;  _DCB13[10] = -2.3894; _DCB13[11] = -1.8875; 
	_DCB13[12] =-1.4572;  _DCB13[13] = 1.0909;

	_DCB12[0] = 16.0264; _DCB12[1] = 5.7435; _DCB12[2] = 5.0819; _DCB12[3] = 3.9562;
	_DCB12[4] =  -0.1673; _DCB12[5] =   1.1228; _DCB12[6] =  4.1743; _DCB12[7] = 2.6981; 
	_DCB12[8] = -5.7684; _DCB12[9] = -4.4165; _DCB12[10] = -7.1783; _DCB12[11] = -3.9031;
	_DCB12[12] =-5.7404; _DCB12[13] = -4.0928;
}

int  PPP::ResPPP(u2 Iter, GNSSInfo &PInfo, const GNSSDATA Obs, const double SatXYZ[], const double SatClk[], 
					   double Azel[], double EphRMS[], double H[], double V[], double R[])
{
	u2 i,j,k,Nv = 0;
	double P,Geometric,Dion=0,Vion=0,Dtrop=0,Vtrop=0,Dtr=0;
	double e[3],pos[3],TropMFunc[3];
	int Nf = 2;
	u2 SysNum,SYS,Apos[3];
	char StrSys;
	double AmbVal,*Var,StaXYZ[3],Disp[3];

	Var = zeros(Obs.NumSats * 2 * 3, 1);
    /*-------------------Earth Tide Correction------------------*/
	memcpy(&StaXYZ, &PInfo.m_X[0], sizeof(double)*3);
	tidedisp(gpst2utc(Obs.ObsTime), StaXYZ, 1, NULL, NULL, Disp);
	for (i=0;i<3;i++) 
	{
		StaXYZ[i]+=Disp[i];
	}

	ecef2pos(StaXYZ, pos);


	for(i=0; i < Obs.NumSats; i++)
	{
		/* geometric distance/azimuth/elevation angle */
		if ((Geometric = geodist(&SatXYZ[i*3], StaXYZ, e)) <= 1e-5 || satazel(pos,e, &Azel[i*2]) < 10 * D2R)
			continue;

		Prn2Sys(Obs.ObsData[i].prn, SYS, SysNum, StrSys);

		
		m_C1 =  SQR(G_Lam[SysNum][1])/(SQR(G_Lam[SysNum][1])-SQR(G_Lam[SysNum][0]));
		m_C2 = -SQR(G_Lam[SysNum][0])/(SQR(G_Lam[SysNum][1])-SQR(G_Lam[SysNum][0]));
		m_C3 =  SQR(G_Lam[SysNum][2])/(SQR(G_Lam[SysNum][2])-SQR(G_Lam[SysNum][0]));
		m_C4 = -SQR(G_Lam[SysNum][0])/(SQR(G_Lam[SysNum][2])-SQR(G_Lam[SysNum][0]));

		GetTropDelay(Obs.ObsTime, pos, &Azel[i*2], Dtrop, Vtrop);

		Dtrop += prectrop(Obs.ObsTime, pos, &Azel[i*2], &PInfo.m_X[IT(0)], TropMFunc, Vtrop);

		Apos[0] = IB(Obs.ObsData[i].prn, 0);
		Apos[1] = IB(Obs.ObsData[i].prn, 1);
		Apos[2] = IB(Obs.ObsData[i].prn, 2);

		for(j = 0; j < Nf*2; j++)  // 0 lc1     1 lc2       2 pc1     3 pc2
		{
			if(j % 2 == 1/* && SYS == SYS_GPS*/)
				continue;


			for(k = 0; k < PInfo.Nx; k++)
			{
				H[k + Nv*PInfo.Nx] = 0.0;
			}

			if((P = GetMeasureMent(Obs, i, SysNum, j/Nf, j%Nf)) < 1e-5)
				continue;

			double Dcb = 0;
			if(j % 2 == 1 && SYS == SYS_BDS)
			{
				Dcb = (m_C2 * _DCB12[Obs.ObsData[i].prn - 33] - m_C4 * _DCB13[Obs.ObsData[i].prn - 33])* 0.29998;
			}

			Dtr   = PInfo.m_X[3+SysNum];
			V[Nv] = P - (Geometric + Dtr - CLIGHT*SatClk[i] + Dion + Dtrop + Dcb);
			Var[Nv] = varerr(Azel[1+i*2], 1, SYS);

			for(k = 0; k <3; k++)
			{
				H[k + Nv*PInfo.Nx] = -e[k];
			}
			if(j % 2 == 1) //the B1B3 combination clock is differ from B1B2
			{
				H[SysNum + 3 + MAXSYS + Nv*PInfo.Nx] = 1.0; //rec clock
			}
			else
			{			
				H[SysNum + 3 + Nv*PInfo.Nx] = 1.0; //rec clock
			}				
			H[IT(0) + Nv*PInfo.Nx] = TropMFunc[0];  //tropshpere

			if(j < Nf)  //carrier phase
			{
				Var[Nv] = varerr(Azel[1+i*2], 0, SYS);
				if(0 == j)    // L1/L2  |  B1/B2
				{
					AmbVal =  PInfo.m_X[Apos[0]];
					V[Nv] += AmbVal;
					H[Apos[0] + Nv*PInfo.Nx] = -1.0;
				}
				else if(1 == j)
				{
					AmbVal =  PInfo.m_X[Apos[1]];
					V[Nv] += AmbVal;
					H[Apos[1] + Nv*PInfo.Nx] = -1.0;
				}
			}	
			Nv++;
		} //j = 0 : nf
	}  // for i =0 : satNum

	int Ref = -1;
	for(i = 0; i < Obs.NumSats; i++)
	{
		Prn2Sys(Obs.ObsData[i].prn, SYS, SysNum, StrSys);
		if(SYS != SYS_BDS)
			continue;
		if(1 == PInfo.SSat[Obs.ObsData[i].prn - 1].FixFlag[1])
		{
			Ref = i;
			break;
		}
	}
	for(i = 0; i < Obs.NumSats; i++)
	{
		Prn2Sys(Obs.ObsData[i].prn, SYS, SysNum, StrSys);
		if(SYS != SYS_BDS)
			continue;
		Apos[0] = IB(Obs.ObsData[i].prn, 0);
		Apos[1] = IB(Obs.ObsData[i].prn, 1);
		int RefPos[3];
		RefPos[0] = IB(Obs.ObsData[Ref].prn, 0);
		RefPos[1] = IB(Obs.ObsData[Ref].prn, 1);

		int prn_1 = Obs.ObsData[i].prn-1;
		int prn_ref =  Obs.ObsData[Ref].prn-1;
		if(1 == PInfo.SSat[prn_1].FixFlag[1] && prn_1 != prn_ref)
		{
			Var[Nv] = 1e-6;
			double xs1,xs2,xs3,xs4;
			double lam_w= (G_Lam[1][0]*G_Lam[1][1]) / (G_Lam[1][1]-G_Lam[1][0]);
			double lam_w2 = (G_Lam[1][0]*G_Lam[1][2]) / (G_Lam[1][2]-G_Lam[1][0]);
			xs1 = (G_Lam[1][0] + G_Lam[1][1])/(G_Lam[1][0] * G_Lam[1][1]);
			xs2 = (G_Lam[1][0] + G_Lam[1][2])/(G_Lam[1][0] * G_Lam[1][2]);
			xs3 = lam_w / G_Lam[1][1];
			xs4 = lam_w2 / G_Lam[1][2];
			V[Nv] = (xs3 - xs4) * (PInfo.SSat[prn_1].FixWL[1] - PInfo.SSat[prn_ref].FixWL[1]) + 
				xs4 * (PInfo.SSat[prn_1].FixWL[0] + PInfo.SSat[prn_ref].FixWL[0]);

			H[Apos[0] + Nv*PInfo.Nx] =   xs1;
			H[RefPos[0] + Nv*PInfo.Nx] =  -xs1;
			H[Apos[1] + Nv*PInfo.Nx] =  -xs2;
			H[RefPos[1] + Nv*PInfo.Nx] = xs2;
			V[Nv] -= (xs1 * (PInfo.m_X[Apos[0]] - PInfo.m_X[RefPos[0]])- xs2 * (PInfo.m_X[Apos[1]] - PInfo.m_X[RefPos[1]]));
			double tEMP = AmbCorr[prn_1-32] - AmbCorr[prn_ref-32];
			V[Nv] -=tEMP;
			//if(prn_1 == 34 && prn_ref == 32)
				//printf("%s %3d %10.3f %10.3f %10.3f\n",PInfo.StaName,prn_1,V[Nv], xs4, V[Nv] - xs4);
			Nv++;
		}
	}

	for(i = 0; i < Nv; i++)
	{
		for(j = 0; j < Nv; j++)
		{
			R[i+j*Nv]=i==j?Var[i]:0.0;
		}
	}
	free(Var);
	return Nv;
}
bool PPP::PPPProcess(GNSSInfo &PInfo, double X[],const GNSSDATA Obs,const double SatXYZ[], const double SatClk[], double SatAzEl[], double ephRMS[])
{
	u2 i,j;
	int Nv=0;
	double *H,*V,*R;
	int info;
	double *Xp,*Pp;

	PreProcess(Obs, PInfo, X);

	Fix_EWL_Amb(Obs, PInfo.SSat, SatAzEl,PInfo.Station[PInfo.N_Rove].Name);
	
	
	//return false;

	memcpy(&PInfo.LastObsTime, &Obs.ObsTime, sizeof(gtime_t));

	H = zeros(PInfo.Nx, Obs.NumSats * 2 * 3);
	V = zeros(Obs.NumSats * 2 * 3, 1);
	R = zeros(Obs.NumSats * 2 * 3, Obs.NumSats * 2 * 3);
	Xp = mat(PInfo.Nx, 1);
	Pp = mat(PInfo.Nx, PInfo.Nx);
	for(i=0; i < 1; i++)
	{
		memcpy(Xp, PInfo.m_X, sizeof(double) * PInfo.Nx);
		memcpy(Pp, PInfo.m_P, sizeof(double) * PInfo.Nx * PInfo.Nx);
		Nv = ResPPP(i, PInfo, Obs, SatXYZ, SatClk, SatAzEl, ephRMS, H, V, R);

		if ((info=filter(Xp, Pp, H, V, R, PInfo.Nx, Nv))) 
		{
			break;
		}
		memcpy(PInfo.m_X, Xp, sizeof(double) * PInfo.Nx);
		memcpy(PInfo.m_P, Pp, sizeof(double) * PInfo.Nx * PInfo.Nx);
	}

	free(H);free(V);free(Xp);free(Pp);free(R);
	return true;
}

void PreProcessAmb(GNSSInfo &PInfo, GNSSDATA Obs,AMB_Info &AmbInfo)
{
	/*-------- detect cycle slip ---------*/
	DetectCS_GF2(Obs, PInfo.SSat);
	DetectCS_GF2_13(Obs, PInfo.SSat);
	DetectCS_MW12(Obs, PInfo.SSat);
	DetectCS_MW23(Obs, PInfo.SSat);
	/*-------- detect cycle slip ---------*/
	int k = 0,i;
	/*-------- reset if sat cannot visible ---------*/
	for( i = 1; i <= MAXSAT; i++)
	{
		if(Obs.ObsData[k].prn == i)
		{
			k++;
		}
		else
		{
			AmbInfo.Numb[i - 1] = 0;
		}
	}
	/*-------- reset if sat cannot visible ---------*/
	/*-------- reset if cycle slip occurs ---------*/
	for(i = 0; i < Obs.NumSats; i++)
	{
		int prn_1 = Obs.ObsData[i].prn-1;
		if(PInfo.SSat[prn_1].CSlip[0] || PInfo.SSat[prn_1].CSlip[1]
		|| PInfo.SSat[prn_1].CSlip[2])
		{
			AmbInfo.Numb[prn_1] = 0;
			AddWLInfo(AmbInfo, Obs.ObsData[i].prn, PInfo.SSat[prn_1].SoD, PInfo.LogSoD);
			PInfo.SSat[prn_1].SoD = PInfo.LogSoD;
		}
		else if(Obs.ObsData[i].P[0] < 10000.0 || Obs.ObsData[i].P[1] < 10000.0 || Obs.ObsData[i].P[2] < 10000.0)
		{
			PInfo.SSat[prn_1].Flag = true;
			AmbInfo.Numb[prn_1] = 0;
		}
		else if(fabs(Obs.ObsData[i].L[0]) < PRECISION || fabs(Obs.ObsData[i].L[1]) < PRECISION
			 || fabs(Obs.ObsData[i].L[2]) < PRECISION)
		{
			PInfo.SSat[prn_1].Flag = true;
			AmbInfo.Numb[prn_1] = 0;
		}
	}
	/*-------- reset if cycle slip occurs ---------*/
}
void PPP::ProcessAmb(GNSSInfo &PInfo, GNSSDATA Obs, double *Azel, AMB_Info &AmbInfo)
{
	/*-------------detect cycle slip and reset count number------------*/
	PreProcessAmb(PInfo, Obs, AmbInfo);

	u2 SYS,SysNum;
	char StrSys;	
	double EWL,WL;
	double Lamda[2][3] = {G_EWL[0],CLIGHT/(FREQ1-FREQ5),0.0,
		                  G_EWL[1],CLIGHT/(FREB1-FREB5),0.0},Coeff[4][5];
	Coeff[0][0] = 0.0;    Coeff[0][1] =  0.488;  Coeff[0][2] = 0.5120;  Coeff[0][3] = 0.0;    Coeff[0][4] = 0.0;
	Coeff[1][0] = 0.4889; Coeff[1][1] = 0.1882;  Coeff[1][2] = 0.3229;  Coeff[1][3] = 0.0;    Coeff[1][4] = 0.0;
	Coeff[2][0] = 0.0429; Coeff[2][1] = 0.6273;  Coeff[2][2] = 0.3298;  Coeff[2][3] = 0.0;    Coeff[2][4] = 0.0;
	Coeff[3][0] = 0.5971; Coeff[3][1] = 0.1479;  Coeff[3][2] = 0.2550;  Coeff[3][3] = 0.0;    Coeff[3][4] = 0.0;
	/*-------- get raw EWL/WL ambiguity ****minus FCB***** ---------*/
	for(int i = 0; i < Obs.NumSats; i++)
	{
		int Prn_1 = Obs.ObsData[i].prn - 1;
		int numb = AmbInfo.Numb[Prn_1];
		if(PInfo.SSat[Prn_1].Flag)
		{
			continue;
		}
		Prn2Sys(Obs.ObsData[i].prn, SYS, SysNum, StrSys);
		if(SYS == SYS_GPS)
		{
//			if(GetMWObs(Obs, i, WL, Coeff[1], 1, Lamda))
			{
				AmbInfo.WL[Prn_1][numb]  = WL  - AmbInfo.FCB[Prn_1][1];
				AmbInfo.Var[Prn_1][numb] = var_LC(0, 1, 1, 0.3, SYS) * (1.0 + 1.0/sin(Azel[1 + 2 * i]) / sin(Azel[1 + 2 * i]));
				AmbInfo.Numb[Prn_1]++;
				AmbInfo.SatFloAmb[Prn_1][0] = WL  - AmbInfo.FCB[Prn_1][1];
			}
		}
		else
		{
			if(GetMWObs(Obs, i, EWL, Coeff[2], 0, Lamda[1]) && GetMWObs(Obs, i, WL, Coeff[3], 1, Lamda[1]))
			{
				AmbInfo.EWL[Prn_1][numb] = EWL - AmbInfo.FCB[Prn_1][0];
				AmbInfo.WL[Prn_1][numb]  = WL  - AmbInfo.FCB[Prn_1][1];
				AmbInfo.Var[Prn_1][numb] = var_LC(0, 1, 1, 0.3, SYS) * (1.0 + 1.0/sin(Azel[1 + 2 * i]) / sin(Azel[1 + 2 * i]));
				AmbInfo.Numb[Prn_1]++;
				AmbInfo.SatFloAmb[Prn_1][0] = WL  - AmbInfo.FCB[Prn_1][1];
			}
		}	
	}
	/*-------- get raw EWL/WL ambiguity ****minus FCB***** ---------*/

	ComputeAR(Obs, 40, AmbInfo);
	return;
}
void ComputeUDAmb(GNSSDATA Obs, AMB_Info &AmbInfo, int Sample, int type, double *FloAmb, 
				  double *Var, int *AmbNum, int *sat, int *SatType)
{
	int i,j,LSys,prn_1,RawNum;
	u2 Sys,SysNum;
	char StrSys;
	double lam_WL[2] = {CLIGHT/(FREQ1 - FREQ2), CLIGHT/(FREB1 - FREB2)};
	double lam_EWL[2] = {CLIGHT/(FREQ2 - FREQ5), CLIGHT/(FREB5 - FREB2)};

	AmbNum[0] = AmbNum[1] = 0;
    AmbInfo.AmbProcess = true;
	if(type == 0) // EWL
	{
		for(LSys = 0; LSys < 2; LSys++)
		{
			int    RecNum = 0;
			double RecBias[MAXCHN] = {0},AveFcb;
			for(i = 0; i < Obs.NumSats; i++)
			{
				int Number = 0;
				RawNum = AmbNum[0] + AmbNum[1];
				Prn2Sys(Obs.ObsData[i].prn, Sys, SysNum, StrSys);
				prn_1 = Obs.ObsData[i].prn - 1;
				if(SysNum != LSys || AmbInfo.Numb[prn_1] < Sample)
				{
					continue;
				}
				for(j = AmbInfo.Numb[prn_1] - Sample; j < AmbInfo.Numb[prn_1]; j++)
				{
					FloAmb[RawNum] +=  AmbInfo.EWL[prn_1][j];
					   Var[RawNum] +=  AmbInfo.Var[prn_1][j];
					   Number++;
				}
				FloAmb[RawNum] = FloAmb[RawNum] / Number;
				   Var[RawNum] = sqrt(Var[RawNum] / (lam_EWL[SysNum] * lam_EWL[SysNum] * Number * Number));
				//if(prn_1 > 36 || prn_1 < 32) //not GEO
				{
					RecBias[RecNum] = FloAmb[RawNum] - (int)(FloAmb[RawNum]);
					RecNum++;
				}
				if(prn_1 > 31 && prn_1 < 37)  //GEO
				{
					SatType[0]++;
				}
				else if(prn_1 > 36 && prn_1 < 42)
				{
					SatType[1]++;
				}
				else if(prn_1 > 41)
				{
					SatType[2]++;
				}
				sat[RawNum] = Obs.ObsData[i].prn;
				AmbNum[LSys]++;
			} //Satllite 
			/* ----------------minus receiver FCB-----------------*/
			if(RecNum > 0)
			{
				if(GetReceiverFCB(RecBias, RecNum, AveFcb))
				{
					if(AveFcb - AmbInfo.RcvFcb[LSys][0] > 0.5)
					{
						AveFcb -= 1.0;
					}
					else if(AveFcb - AmbInfo.RcvFcb[LSys][0] < -0.5)
					{
						AveFcb += 1.0;				
					}
					AmbInfo.RcvFcb[LSys][0] = AveFcb;

					if(LSys == 0)
					{
						for(j = 0; j < AmbNum[LSys]; j++)
						{
							FloAmb[j] -= AveFcb;
						}
					}
					else
					{
						for(j = AmbNum[0]; j < AmbNum[0] + AmbNum[LSys]; j++)
						{
							FloAmb[j] -= AveFcb;
						}
					}
				}
			}
			/* ----------------minus receiver FCB-----------------*/
		}
	}
	else
	{
		for(LSys = 0; LSys < 2; LSys++)
		{
			int    RecNum = 0;
			double RecBias[MAXCHN] = {0},AveFcb;
			for(i = 0; i < Obs.NumSats; i++)
			{
				int Number = 0;
				RawNum = AmbNum[0] + AmbNum[1];
				Prn2Sys(Obs.ObsData[i].prn, Sys, SysNum, StrSys);
				prn_1 = Obs.ObsData[i].prn - 1;
				if(SysNum != LSys || AmbInfo.Numb[prn_1] < Sample || !AmbInfo.Flag[prn_1][0])
				{
					continue;
				}
				for(j = AmbInfo.Numb[prn_1] - Sample; j < AmbInfo.Numb[prn_1]; j++)
				{
					FloAmb[RawNum] +=  AmbInfo.WL[prn_1][j];
					Var[RawNum] +=  AmbInfo.Var[prn_1][j];
					Number++;
				}
				//FloAmb[RawNum] += 0.4631 * (AmbInfo.Amb[prn_1][0] + AmbInfo.FCB[prn_1][0]) * G_EWL[1] / G_WL[1];
 				FloAmb[RawNum] = FloAmb[RawNum] / Number;
				AmbInfo.SatFloAmb[prn_1][1] = FloAmb[RawNum];
				Var[RawNum] = sqrt(Var[RawNum] / (lam_WL[SysNum] * lam_WL[SysNum] * Number * Number));
// 				if(prn_1 > 36 || prn_1 < 32) //not GEO
				{
					RecBias[RecNum] = FloAmb[RawNum] - (int)(FloAmb[RawNum]);
					if(fabs(AmbInfo.RcvFcb[LSys][1]) > 0)
					{
						if(RecBias[RecNum] - AmbInfo.RcvFcb[LSys][1] > 0.5)
						{
							RecBias[RecNum] -= 1.0; 
						}
						else if(RecBias[RecNum] - AmbInfo.RcvFcb[LSys][1] < -0.5)
						{
							RecBias[RecNum] += 1.0; 
						}
					}
					else
					{
						if(RecBias[RecNum] - RecBias[0] > 0.5)
						{
							RecBias[RecNum] -= 1.0; 
						}
						else if(RecBias[RecNum] - RecBias[1] < -0.5)
						{
							RecBias[RecNum] += 1.0; 
						}
					}
					RecNum++;
				}
				if(prn_1 > 31 && prn_1 < 37)  //GEO
				{
					SatType[0]++;
				}
				else if(prn_1 > 36 && prn_1 < 42)
				{
					SatType[1]++;
				}
				else if(prn_1 > 41)
				{
					SatType[2]++;
				}
				sat[RawNum] = Obs.ObsData[i].prn;
				AmbNum[LSys]++;
			} //Satllite 
			/* ----------------minus receiver FCB-----------------*/
			if(RecNum > 0)
			{
				double Temp;
				if(DetectRobust(RecBias, RecNum, 0.12, AveFcb, Temp))
// 				if(GetReceiverFCB(RecBias, RecNum, AveFcb))
				{
					if(AveFcb - AmbInfo.RcvFcb[LSys][1] > 0.5)
					{
						AveFcb -= 1.0;
					}
					else if(AveFcb - AmbInfo.RcvFcb[LSys][1] < -0.5)
					{
						AveFcb += 1.0;				
					}
					AmbInfo.RcvFcb[LSys][1] = AveFcb;
					if(LSys == 0)
					{
						for(j = 0; j < AmbNum[LSys]; j++)
						{
							FloAmb[j] -= AveFcb;
						}
					}
					else
					{
						for(j = AmbNum[0]; j < AmbNum[0] + AmbNum[LSys]; j++)
						{
							FloAmb[j] -= AveFcb;
						}
					}
				}
				else
				{
					AmbInfo.AmbProcess = false;
// 					AmbNum[LSys] = 0;
				}
			}
			/* ----------------minus receiver FCB-----------------*/
		}
	}
}
void PPP::ComputeAR(GNSSDATA Obs, int Sample, AMB_Info &AmbInfo)
{
	double Amb[MAXCHN] = {0},Var[MAXCHN] = {0};
	int sat[MAXCHN] = {0},AmbNum[ENABLEGNSS] = {0};
	bool FixFlag;
	
	memset(AmbInfo.Flag, 0x00, sizeof(bool) * MAXSAT * 2);

	for(int type = 0; type < 2; type++)
	{
		int FixedCount = 0, SatTypeNum[3] = {0};
		memset(Amb,    0x00, sizeof(Amb));
		memset(Var,    0x00, sizeof(Var));
        memset(sat,    0x00, sizeof(sat));
        memset(AmbNum, 0x00, sizeof(AmbNum));
		ComputeUDAmb(Obs, AmbInfo, Sample, type, Amb, Var, AmbNum, sat, SatTypeNum);

		for(int i = 0; i < 3; i++)
		{
			AmbInfo.FixNumber[i][type * 2 + 1] += SatTypeNum[i];
		}

		if(AmbInfo.AmbProcess)
		{
			int SatTypeFix[3] = {0};
			for(int i = 0; i < AmbNum[0] + AmbNum[1]; i++)
			{
				double FloatAmb = Amb[i];
				FixFlag = probability(FloatAmb, Var[i]) > 0.9999;
				if(FixFlag)
				{
					if(FloatAmb < 0.0)
					{
						AmbInfo.Amb[sat[i]-1][type] = int(FloatAmb - 0.5);
					}
					else
					{
						AmbInfo.Amb[sat[i]-1][type] = int(FloatAmb + 0.5);
					}		
					AmbInfo.Flag[sat[i]-1][type] = true;
					if(sat[i] > 32 && sat[i] < 38)
					{
						SatTypeFix[0]++;
					}
					else if(sat[i] > 37 && sat[i] < 43)
					{
						SatTypeFix[1]++;
					}
					else if(sat[i] > 42)
					{
						SatTypeFix[2]++;
					}
					FixedCount++;
				}	
			}
			for(int i = 0; i < 3; i++)
			{
				AmbInfo.FixNumber[i][type * 2] += SatTypeFix[i];
			}
		}
	}
	return;
}
