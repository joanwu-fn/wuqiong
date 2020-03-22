#include "../GNSSProcess/PPPProcess.h"
bool BaseLineProcess(GNSSInfo &BaseLineInfo, const GNSSDATA Obs,const double SatXYZ[], const double SatClk[], double SatAzEl[], double ephRMS[])
{
	int Iter=0,i,j;
	int Nv=0,Nx = BaseLineInfo.Nx;
	double *H,*V,*R;
	int info;
	double *Xp,*Pp;
	double X[7]={0};

	memcpy(X, BaseLineInfo.Station[BaseLineInfo.N_Rove].XYZ, sizeof(double)*3);
	PreProcess(Obs, BaseLineInfo, X);

	//Fix_EWL_Amb(Obs, BaseLineInfo.SSat, SatAzEl,BaseLineInfo.StaName);

	memcpy(&BaseLineInfo.LastObsTime, &Obs.ObsTime, sizeof(gtime_t));

	H = zeros(Nx, Obs.NumSats * 2 * 3);                           V = zeros(Obs.NumSats * 2 * 3, 1);
	R = zeros(Obs.NumSats * 2 * 3, Obs.NumSats * 2 * 3);          Xp = mat(Nx, 1);
	Pp = mat (Nx,  Nx);
	for(Iter = 0; Iter < 1; Iter++)
	{
		memcpy(Xp, BaseLineInfo.m_X, sizeof(double) *  Nx);
		memcpy(Pp, BaseLineInfo.m_P, sizeof(double) *  Nx *  Nx);
		Nv = Get_SingleDiff_Res(Iter, BaseLineInfo, Obs, SatXYZ, SatClk, SatAzEl, ephRMS, H, V, R);
		//Nv = Get_UnDiff_Res(Iter, BaseLineInfo, Obs, SatXYZ, SatClk, SatAzEl, ephRMS, H, V, R);

		if ((info=filter(Xp, Pp, H, V, R, Nx, Nv))) 
		{
			break;
		}
		memcpy(BaseLineInfo.m_X, Xp, sizeof(double) * Nx);
		memcpy(BaseLineInfo.m_P, Pp, sizeof(double) * Nx * Nx);
	}
	free(H);    free(V);    free(R);   free(Xp);
	free(Pp);
	return true;
}