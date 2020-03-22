#include "PosDifferential.h"
#include "SPPProcess.h"

void TransEph(eph_t eph_old, EPHN_t *eph)
{
	eph->Cic             = eph_old.cic;
	eph->Cis             = eph_old.cis;
	eph->clock_bias      = eph_old.f0;
	eph->clock_drift     = eph_old.f1;
	eph->clock_driftrate = eph_old.f2;
	eph->Crc			 = eph_old.crc;
	eph->Crs			 = eph_old.crs;
	eph->Cuc			 = eph_old.cuc;
	eph->Cus			 = eph_old.cus;
	eph->Delta_n		 = eph_old.deln;
	eph->e				 = eph_old.e;
	eph->flags			 = 1;
	eph->i0				 = eph_old.i0;
	eph->IDOT			 = eph_old.idot;
	eph->IODC			 = eph_old.iodc;
	eph->IODE			 = eph_old.iode;
	eph->M0				 = eph_old.M0;
	eph->omega			 = eph_old.omg;
	eph->OMEGA0			 = eph_old.OMG0;
	eph->OMEGADOT		 = eph_old.OMGd;
	eph->prn			 = eph_old.sat;
	eph->sqrt_A			 = eph_old.A;
	eph->SVhealth		 = eph_old.svh;
	eph->TGD[0]			 = eph_old.tgd[0];
	eph->TGD[1]			 = eph_old.tgd[1];
	memcpy(&eph->TOC, &eph_old.toc, sizeof(gtime_t));
	memcpy(&eph->TOE, &eph_old.toe, sizeof(gtime_t));
	eph->toes			 = eph_old.toes;
	eph->URAindex		 = eph_old.sva;
	return;
}
int Processing(int nObs, obsd_t *Obs, int nEph, eph_t *eph,  int NumbSate, int ObserveSate[], double BaseCoor[], double RoveCoor[])
{
	if(nObs < 4 || NumbSate < 4 || nEph < 4)   //Satellite Number is not enough
		return 0;
	int i,j,Number=0;
	/*        Transfer Obs to GNSSDATA      */
	GNSSDATA TransObs;
	
	memcpy(&TransObs.ObsTime, &Obs[0].time, sizeof(gtime_t));	
	for(i = 0; i < nObs; i++)
	{
		for(j = 0; j < NumbSate; j++)
		{
			if(Obs[i].sat == ObserveSate[j])
			{
				memcpy(&TransObs.ObsData[Number].L, &Obs[i].L, sizeof(double)*3);
				memcpy(&TransObs.ObsData[Number].P, &Obs[i].P, sizeof(double)*3);
				TransObs.ObsData[Number].prn = Obs[i].sat;
				Number++;
				break;
			}
		}	
	}
	TransObs.NumSats = Number;
	if(TransObs.NumSats < 4)   //Satellite Number is not enough
		return 0;
	/*        Transfer Obs to GNSSDATA      */
	Orb_Clk  m_Orb;
	for(i = 0; i < nEph; i++)
	{
		TransEph(eph[i], &m_Orb.m_eph[eph[i].sat-1]);
	}
	double Corr[3]={0};
	if( PositionDiff(TransObs, m_Orb, BaseCoor, Corr) == 1)
	{
		for(i = 0; i < 3; i++)
		{
			RoveCoor[i] += Corr[i];
		}
		return 1;
	}
	else
	{
		return 2;
	}
}
int PositionDiff(GNSSDATA Obs, Orb_Clk  m_Orb, double BaseCoor[], double Corr[])
{

	SPP      m_Spp;
	SPPResult Result = {0};
	double *SatXYZ,*Clk,*Azel,*EphRMS;

	SatXYZ = zeros(3, Obs.NumSats);
	Clk    = zeros(1, Obs.NumSats);
	Azel   = zeros(2, Obs.NumSats);
	EphRMS = zeros(1, Obs.NumSats);

	m_Orb.m_type = 0; //brdc;
	m_Orb.GetOrb_Clk(Obs, 0.0, SatXYZ, Clk, Azel, EphRMS, NULL, NULL);
	int i;
	for(i = 0; i < 3; i++)
	{
		Result.x[i] = BaseCoor[i];
		Corr[i] = 0;
	}
	if(m_Spp.SPPProcess(Result, Obs, SatXYZ, Clk, Azel, EphRMS) == true)
	{
		for(i = 0; i < 3; i++)
		{
			Corr[i] = BaseCoor[i] - Result.x[i];
		}
		free(SatXYZ);free(Clk);free(Azel);free(EphRMS);
		return 1;
	}
	else
	{
		free(SatXYZ);free(Clk);free(Azel);free(EphRMS);
		return 0;
	}
	
	
}