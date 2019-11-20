#include "../PreProcess/PreProcess.h"
#include "../Common/debug.h"
#include "../Common/GNSSERROR.h"


double SerialMean(double& mean, double& std, double val, int ind)
{
	if (ind > 1)
	{
		std = std * (ind - 2) / (ind - 1) + ((val - mean) * (val - mean) / ind);
	}
	else
	{
		std = 0;
	}
	mean += (val - mean) / ind;

	return mean;
}

double SerialMean(float& mean, float& std, float val, int ind)
{
	if (ind > 1)
	{
		std = std * (ind - 2) / (ind - 1) + ((val - mean) * (val - mean) / ind);
	}
	else
	{
		std = 0;
	}
	mean += (val - mean) / ind;

	return mean;
}

double Lower(double a, double b)
{
	return (a < b ? a : b);
}

double Upper(double a, double b)
{
	return (a > b ? a : b);
}

MedialDetection::MedialDetection(int maxobs, float* stdCutoff)
{
	_nOmc       = 0;
	_nSort      = 0;

	_maxObs     = maxobs;
	_omcList    = new float    [maxobs];
	_omcSort    = new float    [maxobs];
	_omcWeight  = new float    [maxobs];
	_obsLevel   = new OBS_LEVEL[maxobs];

	int i = 0;
	for (i = 0; i < _maxObs; i++)
	{
		_omcWeight[i] = 1.0;
	}
	_average   = 0.0;
	_std       = 99.0;

	if (NULL != stdCutoff)
	{
		_stdCutoff[0] = stdCutoff[0];
		_stdCutoff[1] = stdCutoff[1];
	}
	else
	{
		_stdCutoff[0] = 0.05;
		_stdCutoff[1] = 100.0;
	}
}

MedialDetection::~MedialDetection()
{
	delete[] _omcList;
	_omcList = NULL;
	delete[] _omcSort;
	_omcSort = NULL;
	delete[] _omcWeight;
	_omcWeight = NULL;
	delete[] _obsLevel;
	_obsLevel = NULL;
}

void MedialDetection::Initialize()
{
	_nOmc     = 0;
	_nSort    = 0;
	_nObsCheck= 0;
	_meanBak  = 0.0;
	_varBak   = 0.0;
	int i = 0;
	for (i = 0; i < _maxObs; i++)
	{
		_omcWeight[i] = 1.0;
	}
	_average = 0.0;
	_std     = 99.0;
}


void MedialDetection::UpdateList(float val, OBS_LEVEL obsLevel)
{
	_omcList [_nOmc] = val;
	_obsLevel[_nOmc] = obsLevel;
	_nObsCheck += (OBS_NOCHECK != _obsLevel[_nOmc]);
	_nOmc++;
	SerialMean(_meanBak, _varBak, val, _nOmc);

	int i = 0, j = 0;
	if (OBS_REFER == obsLevel)
	{
		_omcSort[_nSort] = val;
		if (_nSort == 0 || val >= _omcSort[_nSort - 1])
		{
			_nSort++;
		}
		else
		{
			int low = 0, high = _nSort - 1, m;
			while (low <= high)
			{
				m = (low + high) / 2;
				if (val < _omcSort[m] )
				{
					high = m - 1;
				}
				else 
				{
					low = m + 1;
				}
			}
			memmove(&_omcSort[high + 2], &_omcSort[high + 1], (_nSort - 1 - high) * sizeof(float));
			_omcSort[high + 1] = val;
			_nSort++;
		}
	}
}

bool MedialDetection::ReWeight(int& nNormal, int minObs, int enlarge)
{
	nNormal = 0;
	if (0 == _nObsCheck)
	{
		_average = _meanBak;
		_std = sqrt(_varBak);
		nNormal = _nOmc;
		return true;
	}
	int i = 0, j = 0;
	float mean = 0.0, var = 0.0, meanTemp = 0.0, varTemp = 0.0, weightTotal = 0.0;
	float Q[3] = {0.0, 0.0, 0.0};
	bool bConvergence = false;

	_std = sqrt(_std);
	if (_nSort > 4)
	{
		for (i = 0; i < 3; i++)
		{
			Q[i] = (_nSort + 1) * float(i + 1) / 4.0 - 1;
			if ((_nSort - 1) > int(Q[i]))
			{
				Q[i] = _omcSort[int(Q[i])] + (_omcSort[int(Q[i]) + 1] - _omcSort[int(Q[i])]) * (Q[i] - int(Q[i]));
			}
			else
			{
				Q[i] = _omcSort[int(Q[i])];
			}
		}
		meanTemp = Q[1];
		varTemp  = Upper((Q[2] - Q[0]) / 1.349, _stdCutoff[0]);
		varTemp  = Lower(varTemp, _stdCutoff[1]);
		do 
		{
			bConvergence = true;
			for (i = 0; i < _nSort; i++)
			{
				if (fabs((_omcSort[i] - meanTemp) / varTemp) >= 4.0)
				{
					bConvergence = false;

					_nSort--;
					memmove(_omcSort + i, _omcSort + (i + 1), sizeof(float) * (_nSort - i));
					Q[0] = (_nSort + 1) * 0.25 - 1;
					Q[1] = (_nSort + 1) * 0.50 - 1;
					Q[2] = (_nSort + 1) * 0.75 - 1;
					Q[0] = _omcSort[int(Q[0])] + (_omcSort[int(Q[0]) + 1] - _omcSort[int(Q[0])]) * (Q[0] - int(Q[0]));
					Q[1] = _omcSort[int(Q[1])] + (_omcSort[int(Q[1]) + 1] - _omcSort[int(Q[1])]) * (Q[1] - int(Q[1]));
					Q[2] = _omcSort[int(Q[2])] + (_omcSort[int(Q[2]) + 1] - _omcSort[int(Q[2])]) * (Q[2] - int(Q[2]));
					meanTemp = Q[1];
					varTemp  = Upper((Q[2] - Q[0]) / 1.349, _stdCutoff[0]);
					varTemp  = Lower(varTemp, _stdCutoff[1]);
					break;
				}
			}
		} while (!bConvergence);
	}
	else
	{
		meanTemp = _meanBak;
		varTemp  = sqrt(_varBak);
		varTemp  = Upper(varTemp, _stdCutoff[0]);
		varTemp  = Lower(varTemp, _stdCutoff[1]);
	}

	int iTimer   = 0;
	int normal   = 0;
	int factor = 1, lowLevel = 2.5;
	do 
	{
		mean = meanTemp;
		var  = varTemp;

		weightTotal = 0.0;
		for (i = 0; i < _nOmc; i++)
		{
			if (OBS_NOCHECK == _obsLevel[i])
			{
				continue;
			}
			else
			{
				factor = (OBS_LOWLEVEL == _obsLevel[i] ? lowLevel : 1) * enlarge;
				float v = 1e16;
				if (varTemp > 0.0)
				{
					v = fabs((_omcList[i] - meanTemp) / varTemp);
				}
				if (v < (3.0 * factor))
				{
					_omcWeight[i] = 1.0;
				}
				else if (v < (5.0 * factor))
				{
					float diff = (5.0 * factor - v) / (2.0 * factor);
					_omcWeight[i] = pow(diff, 4);
				}
				else
				{
					_omcWeight[i] = 1e-9;
				}
			}
			weightTotal += (OBS_LOWLEVEL == _obsLevel[i] ? _omcWeight[i] / lowLevel : _omcWeight[i]);
		}

		meanTemp = 0.0;
		for (i = 0; i < _nOmc; i++)
		{
			if (OBS_NOCHECK == _obsLevel[i])
			{
				continue;
			}
			float weightTemp = _omcWeight[i];
			if (OBS_LOWLEVEL == _obsLevel[i])
			{
				weightTemp /= lowLevel;
			}
			if (fabs(_omcWeight[i]) < 1e-9)
			{
				weightTemp = 0;
			}
			meanTemp += (_omcList[i] * weightTemp / weightTotal);
		}

		varTemp = 0.0;
		int nTemp = 0;
		normal = 0;
		for (i = 0; i < _nOmc; i++)
		{
			if (OBS_NOCHECK == _obsLevel[i])
			{
				normal++;
				continue;
			}
			float weightTemp = _omcWeight[i];
			if (OBS_LOWLEVEL == _obsLevel[i])
			{
				weightTemp /= lowLevel;
			}
			if (fabs(_omcWeight[i]) < 1e-9)
			{
				weightTemp = 0;
			}
			else
			{
				normal++;
			}
			varTemp += ((_omcList[i] - meanTemp) * weightTemp * (_omcList[i] - meanTemp));
			nTemp++;
		}
		if (nTemp < 2)
		{
			break;
		}
		varTemp = Upper(sqrt(varTemp / (nTemp - 1)), _stdCutoff[0]);
		varTemp = Lower(varTemp, _stdCutoff[1]);

		iTimer++;
		if (iTimer > 20)
		{
			nNormal = 0;
			return false;
		}
	} while(iTimer < 3 || (fabs(var - varTemp) / var > 5e-2) || ((fabs(mean - meanTemp) / meanTemp > 0.5) && fabs(mean - meanTemp) > 1.0));
	_average = meanTemp;

	nNormal = normal;
	if (varTemp > (_stdCutoff[1] - 0.0001))
	{
		normal = 0;
	}
#ifndef MULTNAVSYS
	if (varTemp > (_stdCutoff[1] - 0.0001)  || normal < minObs)
#else
	if (varTemp > (_stdCutoff[1] - 0.0001) || normal < int(_nOmc / 2.0 + 0.5))
#endif
	{
		_average = 0.0;
		for (i = 0; i < _nOmc; i++)
		{
			if (OBS_NOCHECK == _obsLevel[i])
			{
				continue;
			}
			_omcWeight[i]  = 1.0;
			_average += _omcList[i];
		}
		_average /= _nOmc;
		return false;
	}
	_std = varTemp;

	return true;
}

float MedialDetection::GetOMCWei(int ind)
{
	return _omcWeight[ind];
}

float MedialDetection::GetValue(int ind)
{
	return _omcList[ind] - _average;
}

int MedialDetection::GetNumObs()
{
	return _nOmc;
}

float MedialDetection::GetAverage()
{
	return _average;
}

float MedialDetection::GetSTD()
{
	return _std;
}

float MedialDetection::GetPercent(int percent)
{
	float Q = (_nSort + 1) * percent / 100.0 - 1;
	if ((_nSort - 1) > int(Q))
	{
		Q = _omcSort[int(Q)] + (_omcSort[int(Q) + 1] - _omcSort[int(Q)]) * (Q - int(Q));
	}
	else
	{
		Q = _omcSort[int(Q)];
	}

	return Q;
}


void AddOneFcb(GNSSInfo &PInfo, int Prn)
{
	M_FCB *fcbtemp;
	if(PInfo.Nfcb[0] >= PInfo.Nfcb[1])
	{
		PInfo.Nfcb[1] += 100;
		if((fcbtemp = (M_FCB *)realloc(PInfo.Fcb, sizeof(M_FCB)*PInfo.Nfcb[1])) == NULL)
		{
			ASSERT(4 < 0);
		}
		PInfo.Fcb = fcbtemp;
	}
	PInfo.Fcb[PInfo.Nfcb[0]].FCB[0]    = PInfo.SSat[Prn-1].FCB[0];
	PInfo.Fcb[PInfo.Nfcb[0]].FCB[1]    = PInfo.SSat[Prn-1].FCB[1];
	PInfo.Fcb[PInfo.Nfcb[0]].FCBNum[0] = PInfo.SSat[Prn-1].FCBNum[0];
	PInfo.Fcb[PInfo.Nfcb[0]].FCBNum[1] = PInfo.SSat[Prn-1].FCBNum[1];
	PInfo.Fcb[PInfo.Nfcb[0]++].Prn     = Prn;
}
void ComputeNewAmb(const GNSSDATA &Obs, double *Azel, GNSSInfo &PInfo, FILE* ftMW)
{
	double dMW[MAXSAT][2] = {0};
	double P23,P12,P13,L12,L23,L13,EWL,WL,Lam23,Lam12,Lam13;
	double P[3],L[3];
	int Prn_1;
	char StrSys;
	u2 NavSys,SYS;
	
	for(int i = 0; i < Obs.NumSats; i++)
	{
		int Prn = Prn2Sys(Obs.ObsData[i].prn, SYS, NavSys, StrSys);
// 		if(Prn == 15)
// 		{
// 			continue;
// 		}
// 		if(Prn <= 0 || Prn <= 5)
// 		{
// 			continue;
// 		}
		memcpy(P, Obs.ObsData[i].P, sizeof(double) * 3);
		memcpy(L, Obs.ObsData[i].L, sizeof(double) * 3);
		for(int j = 0; j < 3; j++)
		{
			P[j] = P[j] < 1e-9 ? Obs.ObsData[i].C[j] : P[j];
		}
		if(Azel[i * 2 +1] <= 20 * D2R)
		{
			continue;
		}
		if(P[0] < 10000.0 || fabs(L[0]) < 1e-5)
		{
			continue;
		}
// 		if(P[1] < 10000.0 || fabs(L[1]) < 1e-5)
// 		{
// 			continue;
// 		}
		if(P[2] < 10000.0 || fabs(L[2]) < 1e-5)
		{
 			continue;
		}
		
		Prn_1 = Obs.ObsData[i].prn - 1;
		double cof[2] = {0.0};
		if(NavSys == 0) //GPS
		{
			Lam23 = CLIGHT / (FREQ2 - FREQ5);	
			Lam12 = CLIGHT / (FREQ1 - FREQ2);
			cof[1] = (FREQ1 * FREQ1 / (FREQ2 * FREQ5) - 1.0) / (FREQ1 * FREQ1 / (FREQ2 * FREQ2) - 1.0);
			cof[0] = 1 - cof[1];
			P23   = (P[0] * cof[0] + P[1] * cof[1]);
 			P23   = (P[1] * FREQ2 + P[2] * FREQ5) / (FREQ2 + FREQ5);
 			L23   = (L[1] - L[2]);	
			P12   = (P[0] * FREQ1 + P[1] * FREQ2) / (FREQ1 + FREQ2);	
			L12   = (L[0] - L[1]);
		}
		else if(NavSys == 1) //BDS
		{
			Lam23 = CLIGHT / (FREB5 - FREB2);
			Lam12 = CLIGHT / (FREB1 - FREB2);
			Lam13 = CLIGHT / (FREB1 - FREB5);

			cof[1] = (FREB1 * FREB1 / (FREB2 * FREB5) - 1.0) / (FREB1 * FREB1 / (FREB2 * FREB2) - 1.0);
			cof[0] = 1 - cof[1];
			P23   = (P[0] * cof[0] + P[1] * cof[1]);
			P23   = (P[1] * FREB2 + P[2] * FREB5) / (FREB2 + FREB5);
			L23   = (L[2] - L[1]);

 			P12   = (P[0] * FREB1 + P[1] * FREB2) / ((FREB1 + FREB2));	
			L12   = (L[0] - L[1]);

			P13   = (P[0] * FREB1 + P[2] * FREB5) / ((FREB1 + FREB5));	
			L13   = (L[0] - L[2]);
		}	
		else if(NavSys == 2) //GAL
		{
			Lam23 = CLIGHT / (FREE7 - FREE5);
			Lam12 = CLIGHT / (FREE1 - FREE5);
			cof[1] = (FREE1 * FREE1 / (FREE5 * FREE7) - 1.0) / (FREE1 * FREE1 / (FREE5 * FREE5) - 1.0);
			cof[0] = 1 - cof[1];
			P23   = (P[0] * cof[0] + P[1] * cof[1]);

// 			P23   = (P[1] * FREE5 + P[2] * FREE7) / (FREE5 + FREE7);
			P12   = (P[0] * FREE1 + P[1] * FREE5) / (FREE1 + FREE5);	
			L23   = (L[2] - L[1]);
			L12   = (L[0] - L[1]);
		}
		else
		{
			continue;
		}

		EWL  = P23 / Lam23 - L23;
// 		WL   = P12 / Lam12 - L12;
		WL   = P13 / Lam13 - L13;
		dMW[Prn_1][0] = EWL;
		dMW[Prn_1][1] = WL;
		PInfo.SSat[Prn_1].FCB[0] = (PInfo.SSat[Prn_1].FCB[0] * PInfo.SSat[Prn_1].FCBNum[0] + EWL) 
			/ (PInfo.SSat[Prn_1].FCBNum[0] + 1);
		PInfo.SSat[Prn_1].FCB[1] = (PInfo.SSat[Prn_1].FCB[1] * PInfo.SSat[Prn_1].FCBNum[1] +  WL)
			/ (PInfo.SSat[Prn_1].FCBNum[1] + 1);
		PInfo.SSat[Prn_1].FCBNum[0]++;
		PInfo.SSat[Prn_1].FCBNum[1]++;
	}
	if(ftMW != NULL)
	{
		for(Prn_1 = 0; Prn_1 < MAXPRNGPS + MAXPRNBDS; Prn_1++)
		{
			fprintf(ftMW,"%11.3f	", dMW[Prn_1][0]);
		}
		fprintf(ftMW,"\n");
	}
}
void AddGFIF(const GNSSDATA &Obs, double *Azel, GNSSInfo &PInfo, FILE* ftMW)
{
	double dMW[MAXSAT] = {0.0};
	for(int i = 0; i < Obs.NumSats; i++)
	{
		int Prn_1 = Obs.ObsData[i].prn - 1;
		double P[3], L[3];
		memcpy(P, Obs.ObsData[i].P, sizeof(double) * 3);
		memcpy(L, Obs.ObsData[i].L, sizeof(double) * 3);
		if(Azel[i * 2 +1] <= 20 * D2R || Obs.ObsData[i].prn <= NSATGPS || Obs.ObsData[i].prn > NSATGPS + NSATBDS)
		{
			continue;
		}
		if(P[0] < 10000.0 || P[1] < 10000.0 || P[2] < 10000.0)
		{
			continue;
		}
		double coef[4];
		coef[0] = FREB1 * FREB1 / (FREB1 * FREB1 - FREB2 * FREB2);
		coef[1] = 1 - coef[0];
		coef[2] = FREB1 * FREB1 / (FREB1 * FREB1 - FREB5 * FREB5);
		coef[3] = 1 - coef[2];
		double PIFGF = (coef[0] - coef[2]) * P[0] + coef[1] * P[1] - coef[3] * P[2];
		if(PInfo.iPIFGF[Prn_1] > 100 && fabs(PInfo.dPIFGF[Prn_1] - PIFGF) > 4.0)
		{
			continue;
		}
		dMW[Prn_1] = PIFGF;
		PInfo.dPIFGF[Prn_1] = (PInfo.dPIFGF[Prn_1] * PInfo.iPIFGF[Prn_1] + PIFGF) / (PInfo.iPIFGF[Prn_1] + 1);
		PInfo.iPIFGF[Prn_1]++;
	}
// 	if(ftMW != NULL)
// 	{
// 		for(int Prn_1 = 32; Prn_1 < 46; Prn_1++)
// 		{
// 			fprintf(ftMW,"%11.3f	", dMW[Prn_1]);
// 		}
// 		fprintf(ftMW,"\n");
// 	}
}
void EstFCB(const GNSSDATA &Obs, double *Azel, GNSSInfo &PInfo, FILE* ftMW)
{
	int Prn_1,StaNo = PInfo.N_Rove;
	int SatFlag[MAXSAT] = {0};

	/*-------- detect cycle slip ---------*/
	DetectCS_GF2(Obs, PInfo.SSat);
	DetectCS_GF2_13(Obs, PInfo.SSat);
	DetectCS_MW12(Obs, PInfo.SSat);
	DetectCS_MW23(Obs, PInfo.SSat);
	/*-------- detect cycle slip ---------*/
    
	for(int i = 0; i < Obs.NumSats; i++)
	{
		Prn_1 = Obs.ObsData[i].prn - 1;
		int MinArc = 200;
		if(Prn_1 >= NSATGPS && Prn_1 < NSATGPS + 5) //GEO
		{
			MinArc = 600;
		}
		if((Prn_1 >= NSATGPS + 5 && Prn_1 < NSATGPS + 10) || (Prn_1 == NSATGPS + 12)) //IGSO
		{
			MinArc = 400;
		}
		if(Azel[i*2+1] > 20 * D2R)
		{
			SatFlag[Prn_1] = 1;
		}
		if(PInfo.SSat[Prn_1].CSlip[0] == true || PInfo.SSat[Prn_1].CSlip[1] == true || PInfo.SSat[Prn_1].CSlip[2] == true)
		{
			if(PInfo.SSat[Prn_1].FCBNum[0] > MinArc || PInfo.SSat[Prn_1].FCBNum[1] > MinArc)
			{
				AddOneFcb(PInfo, Prn_1+1);
			}
			PInfo.SSat[Prn_1].FCB[0]    = 0;
			PInfo.SSat[Prn_1].FCB[1]    = 0;
			PInfo.SSat[Prn_1].FCBNum[0] = 0;
			PInfo.SSat[Prn_1].FCBNum[1] = 0;
		}
	}
	for(Prn_1 = 0; Prn_1 < MAXSAT; Prn_1++)
	{
		int MinArc = 200;
		if(Prn_1 >= NSATGPS && Prn_1 < NSATGPS + 5) //GEO
		{
			MinArc = 600;
		}
		if((Prn_1 >= NSATGPS + 5 && Prn_1 < NSATGPS + 10) || (Prn_1 == NSATGPS + 12)) //IGSO
		{
			MinArc = 400;
		}
		if(SatFlag[Prn_1] == 0 || PInfo.SSat[Prn_1].FCBNum[1] > 2880)
		{
			if(PInfo.SSat[Prn_1].FCBNum[0] > MinArc || PInfo.SSat[Prn_1].FCBNum[1] > MinArc)
			{
				AddOneFcb(PInfo, Prn_1+1);
				PInfo.SSat[Prn_1].FCB[0]    = 0;
				PInfo.SSat[Prn_1].FCB[1]    = 0;
				PInfo.SSat[Prn_1].FCBNum[0] = 0;
				PInfo.SSat[Prn_1].FCBNum[1] = 0;
			}
			PInfo.SSat[Prn_1].FCB[0]    = 0;
			PInfo.SSat[Prn_1].FCB[1]    = 0;
			PInfo.SSat[Prn_1].FCBNum[0] = 0;
			PInfo.SSat[Prn_1].FCBNum[1] = 0;
		}	
	}
	ComputeNewAmb(Obs, Azel, PInfo, ftMW);
	return;
}
void AddLastTime(GNSSInfo &PInfo)
{
	for(int Prn_1 = 0; Prn_1 < MAXSAT; Prn_1++)
	{
		int MinArc = 200;
		if(Prn_1 >= NSATGPS && Prn_1 < NSATGPS + 5) //GEO
		{
			MinArc = 600;
		}
		if((Prn_1 >= NSATGPS + 5 && Prn_1 < NSATGPS + 10) || (Prn_1 == NSATGPS + 12)) //IGSO
		{
			MinArc = 400;
		}

		if(PInfo.SSat[Prn_1].FCBNum[0] > MinArc || PInfo.SSat[Prn_1].FCBNum[1] > MinArc)
		{
			AddOneFcb(PInfo, Prn_1+1);
			PInfo.SSat[Prn_1].FCB[0]    = 0;
			PInfo.SSat[Prn_1].FCB[1]    = 0;
			PInfo.SSat[Prn_1].FCBNum[0] = 0;
			PInfo.SSat[Prn_1].FCBNum[1] = 0;
		}
	}
	return;
}


void ComputeAverage(UPDInfo* UPDIf, GNSSInfo *PInfo, int StaNum, int Type, FILE *ft, int nSat[2])
{
	int StaNo,i,j,prn_1;
	double FCBONESta[MAXSAT][100] = {0}, Fcb;
	int    FCBNum[MAXSAT] = {0};
	int StartPrn,EndPrn;

	StartPrn = nSat[0];
	EndPrn = nSat[1];
    //================================take average of different arcs for each station
	for(StaNo = 0; StaNo < StaNum; StaNo++)
	{
		memset(FCBNum, 0x00, sizeof(int) * MAXSAT);
		double RefFCB[MAXSAT] = {0}, FCBRMS[MAXSAT] = {0.0};
		for(i = 0; i < PInfo[StaNo].Nfcb[0]; i++)
		{
			prn_1 = PInfo[StaNo].Fcb[i].Prn - 1;
			double Fcbtemp1 = PInfo[StaNo].Fcb[i].FCB[Type] - int(PInfo[StaNo].Fcb[i].FCB[Type]);
			if(FCBNum[prn_1] == 0)
			{
				FCBONESta[prn_1][FCBNum[prn_1]] = Fcbtemp1;
				RefFCB[prn_1] = Fcbtemp1;
				FCBRMS[prn_1] = Fcbtemp1 * Fcbtemp1;
			}
			else
			{
				while(fabs(Fcbtemp1 - RefFCB[prn_1]) > 0.5)
				{
					if(Fcbtemp1 - RefFCB[prn_1] > 0.5)
					{
						Fcbtemp1 -= 1.0;
					}
					else
					{
						Fcbtemp1 += 1.0;
					}
				}
				FCBONESta[prn_1][FCBNum[prn_1]] = Fcbtemp1;
				RefFCB[prn_1] = (RefFCB[prn_1] * FCBNum[prn_1] + Fcbtemp1)/(FCBNum[prn_1] + 1);
				FCBRMS[prn_1] = (FCBRMS[prn_1] * FCBNum[prn_1] + Fcbtemp1 * Fcbtemp1)/(FCBNum[prn_1] + 1);
			}
			FCBNum[prn_1]++;
		}
#ifdef FCB_DETAIAL
		fprintf(ft,"-------------------Raw FCB of each station-----------------------\n");
		fprintf(ft,"%s\n",PInfo[StaNo].Station[StaNo].Name);
#endif
		double Temp = 0;
		for(i = StartPrn; i < EndPrn; i++)
		{
#ifdef FCB_DETAIAL
			if(i <= 31)
			{
				fprintf(ft,"G%02d:  ",i+1);
			}
			else
			{
				fprintf(ft,"C%02d:  ",i-31);
			}
			for(j = 0; j < FCBNum[i]; j++)
			{
				fprintf(ft,"%6.3f ",FCBONESta[i][j]);
			}
			fprintf(ft,"\n");
#endif
			if(DetectRobust(FCBONESta[i], FCBNum[i], 0.3, PInfo[StaNo].UDFcb[i][Type],Temp))//need modification 0.05
			{
				PInfo[StaNo].UDFcbFlag[i][Type] = true;
			}
		}
#ifdef FCB_DETAIAL
		fprintf(ft,"\n\n");
#endif
	}
	//================================take average of different arcs for each station
	UPDModel pUPDModel;
	pUPDModel.ComputeStatistical(UPDIf, PInfo, StaNum, Type, ft, nSat);

// 	for(prn_1 = StartPrn; prn_1 < EndPrn; prn_1++)
// 	{	
// 		if(Type == 0)
// 		{
// 			UPDIf->_dEWLSeri[UPDIf->_nDay * MAXSAT + prn_1] = FinalFCB[prn_1] - FCBTMP;
// 			UPDIf->_bEWLSeri[UPDIf->_nDay * MAXSAT + prn_1] = bFinalFCB[prn_1];
// 		}
// 		else if(Type == 1)
// 		{
// 			UPDIf->_dWLSeri[UPDIf->_nDay * MAXSAT + prn_1] = FinalFCB[prn_1] - FCBTMP;
// 			UPDIf->_bWLSeri[UPDIf->_nDay * MAXSAT + prn_1] = bFinalFCB[prn_1];
// 		}
// 	}
// 	UPDIf->_nDay++;

	return;
}