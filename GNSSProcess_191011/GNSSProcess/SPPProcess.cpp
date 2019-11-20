#include "SPPProcess.h"
#include "../Common/BaseFunction.h"
#include "../Common/debug.h"
SPP::SPP()
{
	tec.tec = NULL;
	tec.nt  = 0;
	tec.nmax = 0;
}
bool SPP::SPPProcess(SPPResult &Result,const GNSSDATA Obs,const double SatXYZ[], const double SatClk[], double SatAzEl[], double ephRMS[])
{
	return EstPos(Result,Obs, SatXYZ, SatClk, SatAzEl, ephRMS);
}

bool SPP::EstPos(SPPResult &Result,const GNSSDATA Obs, const double SatXYZ[], const double SatClk[], double Azel[], double EphRMS[])
{
	u2 i,j,k,Nv=0,Nx=0,RecClk[4] = {0};
	double *H,*V,*Var,x[7]={0},dx[7],Q[49],sig;
	int info;

	double DOP_B[100],DOP_E[9],DOP_POS[3],DOP_BB[100],DOP[9];

	H   = zeros(Obs.NumSats, 7); 
	V   = zeros(Obs.NumSats, 1);
	Var = zeros(Obs.NumSats, 1);
	for(i=0; i<7; i++)
	{
		x[i] = Result.x[i];
	}
	for(i = 0; i < MAXITER; i++)
	{
		Nv = ResCode(i, x, Obs, SatXYZ, SatClk, Azel, EphRMS, H, V, Var, Nx, RecClk);
		if(Nv < Nx)
		{
			break;
		}
		if(i == 0)
		{
			for(j = 0; j < Nv; j++)
			{
				for(k = 0; k < 3; k++)
				{
					DOP_B[k + j * 3] = H[k + j * Nx];
				}
			}
// 			ecef2pos(x, DOP_POS);
// 			xyz2enu(DOP_POS, DOP_E);
// 			matmul("NT",3, Nv, 3, 1.0, DOP_E, DOP_B, 0.0, DOP_BB);
			matmul("NT",3, 3, Nv, 1.0, DOP_B, DOP_B, 0.0, DOP);
			matinv(DOP, 3);										
			Result.Q[0] = sqrt(DOP[0]);
			Result.Q[1] = sqrt(DOP[4]);
			Result.Q[2] = sqrt(DOP[8]);
		}

		/* weight by variance */
		for (j = 0; j < Nv; j++) 
		{
			sig=sqrt(Var[j]);
			V[j]/=sig;

			for (k = 0; k < Nx; k++)
			{
				H[k+j*Nx]/=sig;
			}	
		}

		/* least square estimation */
		if ((info=lsq(H, V, Nx, Nv, dx, Q))) 
		{
			break;
		}
		for(j=0; j<Nx; j++)
		{
			x[j] += dx[j];
		}
		if(norm(dx, Nx)<1E-4)
		{
			if(valsol(Azel, V, Nv, Nx) == 1)
			{
				memcpy(&Result.x,&x,sizeof(double)*7);
				//memcpy(&Result.Q, &Q, sizeof(double)*49);
				free(H);free(V);free(Var);
				return true;
			}
			else
			{
				break;
			}
			
		}
		u2 ClkPos = 0;
		for(j = 0; j < 4; j++)
		{
			if(RecClk[j] == 1)
			{
				x[3+j] = x[3 + ClkPos];
				ClkPos++;
			}
		}
		ASSERT((ClkPos + 3) == Nx);
	}

	free(H);free(V);free(Var);
	return false;
}

u2 SPP::ResCode(u2 Iter, double x[], const GNSSDATA Obs, const double SatXYZ[], const double SatClk[], 
				double Azel[], double EphRMS[], double H[], double V[],double Var[], u2 &Nx, u2 RecCLK[])
{
	u2 i,j,nv=0;
	double P,Geometric,Dion=0,Vion,Dtrop,Vtrop,Dtr;
	double pos[3],e[3];
	u2 SYS,SysNum;
	u2 ParaNum=0;
	char StrSys;

	Nx=3;
	ecef2pos(x,pos);
	
	RecCLK[0] = RecCLK[1] = RecCLK[2] = RecCLK[3] = 0;
	for(i = 0; i < Obs.NumSats; i++)
	{
		Azel[2*i] = Azel[1+2*i] = 0.0;
	
		/* geometric distance/azimuth/elevation angle */
		if ((Geometric = geodist(&SatXYZ[i*3], x, e)) <= 1e-5 || satazel(pos,e, &Azel[i*2]) < 10 * D2R)
		{
			continue;
		}
		if(fabs(CLIGHT * SatClk[i]) < 0.1)
		{
			continue;
		}

		Prn2Sys(Obs.ObsData[i].prn, SYS, SysNum, StrSys);
		/* psudorange with code bias correction */
		if ((P = Prange(Obs, i, Azel[i*2+1], SysNum)) <= 1e-5) 
		{
			continue;
		}

		GetTropDelay(Obs.ObsTime, pos, &Azel[i*2], Dtrop, Vtrop);
		if(m_Type == SPPMODE_SINGLE)
		{
			iontec(Obs.ObsTime, pos, tec.tec, tec.nt, &Azel[i*2], 1, &Dion, &Vion);
		}
		
        Dtr = x[3+SysNum];

		V[nv] = P - (Geometric + Dtr - CLIGHT*SatClk[i] + Dion + Dtrop);
		
		Var[nv] = EphRMS[i] + varerr(Azel[1+i*2], 1, SYS) ;

		RecCLK[SysNum] = 1;
		/* design matrix */
		for (j = 0; j < 7; j++) 
		{
			H[j + nv * 7] = j < 3?-e[j]:0.0;
		}
		H[3+ SysNum + nv*7] = 1.0;
		nv++;
	}
	for(i=0; i < 4; i++)
	{
		if(RecCLK[i] == 1)
		{
			Nx++;
		}
	}

	for(i = 0; i< nv; i++)
	{
		ParaNum=0;
		for(j=0; j<3; j++)
		{
			H[j+ i*Nx] = H[j + i*7];
		}
		for(j = 3; j < 7; j++)
		{
			if(RecCLK[j - 3] == 1)
			{
				H[3 + ParaNum + i * Nx] = H[j + i*7];
				ParaNum++;
			}
		}
	}
	return nv;
}

double SPP::Prange(const GNSSDATA Obs, u2 Numb, double Satel, u2 SysNum)
{
	double P1,P2,Bias[4]={0};
	double C1,C2;

	if(m_Type == SPPMODE_IFLC)  //ionoshpere free
	{
		
		C1 =  SQR(G_Lam[SysNum][1])/(SQR(G_Lam[SysNum][1])-SQR(G_Lam[SysNum][0]));
		C2 = -SQR(G_Lam[SysNum][0])/(SQR(G_Lam[SysNum][1])-SQR(G_Lam[SysNum][0]));

		P1 = Obs.ObsData[Numb].P[0];
		P2 = Obs.ObsData[Numb].P[1];

		if(fabs(P1) < 2000000.0 || fabs(P2) < 2000000.0)
		{
			return 0.0;
		}
		return C1*(P1) + C2*(P2);
	}
	else
	{
		P1 = Obs.ObsData[Numb].P[0];
		if(Obs.ObsData[Numb].prn > 32)
		{
			Bias[0] = 1.48 * G_B1B2Default[Obs.ObsData[Numb].prn - 33] * 0.2998;
		}
		if(fabs(P1) < 1e-5)
		{
			return 0.0;
		}
		return P1 + Bias[0];
	}
}