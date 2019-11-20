#include "../PreProcess/PreProcess.h"
#include "../Common/debug.h"
#include "../Common/GNSSERROR.h"
#define UPDCONSISTENCE  0            //0: all satellites,    1: BDS2 satelltes      2: BDS3 satellites

int sign(double val)
{
	return (val >= 0 ? 1 : -1);
}

void UPDModel::FindStation(int &iRefSite, bool bFlag)
{
	iRefSite = -1;
	int StaNo = 0, prn_1 = 0, MaxSV = 0;
	//============================Select Max SV station for reference station
	for(StaNo = 0; StaNo < _nSite; StaNo++)
	{
		if(_bIniSite[StaNo])
		{
			continue;
		}
		int SVNum = 0;
		for(prn_1 = _iPrnNum[0]; prn_1 < _iPrnNum[1]; prn_1++)
		{
#if (UPDCONSISTENCE == 2)
			if(prn_1 < NSATGPS + 16)
			{
				continue;
			}
#else if (UPDCONSISTENCE == 1)
			if(prn_1 >= NSATGPS + 16)
			{
				continue;
			}
#endif
			if(_PInfo[StaNo].UDFcbFlag[prn_1][_iType] && (bFlag || _bIniSat[prn_1]))
			{
				SVNum++;			
			}
		}
		if(SVNum > MaxSV)
		{
			MaxSV     = SVNum;
			iRefSite = StaNo;
		}
	}
	return;
}
void UPDModel::GetIniUPD()
{
	int iSite = -1;
	FindStation(iSite, true);

	if(-1 != iSite)
	{
		fprintf(_ftOut,"  %s ,", _PInfo[iSite].Station[iSite].Name);
		int prn_1 = 0;
		for(prn_1 = _iPrnNum[0]; prn_1 < _iPrnNum[1]; prn_1++)
		{
			if(_PInfo[iSite].UDFcbFlag[prn_1][_iType])
			{
				_dSatUPD[prn_1] = _PInfo[iSite].UDFcb[prn_1][_iType];
				_bIniSat[prn_1] = true;	
				fprintf(_ftOut,"%7.4f ,", _PInfo[iSite].UDFcb[prn_1][_iType]);
			}
			else
			{
				fprintf(_ftOut,"        ,");
			}
		}
		fprintf(_ftOut,"\n");
		_bIniSite[iSite] = true;
	}
	return;
}

bool UPDModel::ExamineConsistence(MedialDetection* pMedia1, MedialDetection* pMedia2, int iStaNo, int& iMediaNum, int prnlist[MAXSAT])
{
	int nNormal = 0;
	int iValid[2] = {0};
	bool bFlag[MAXSAT][2];
	memset(bFlag, 0x00, sizeof(bFlag));
	if (pMedia1->ReWeight(nNormal, 1))
	{
		for (int i = 0; i < pMedia2->_nOmc; i++)
		{
			if (pMedia1->GetOMCWei(i) > 0.99)
			{
				iValid[0]++;
			}
			else
			{
				bFlag[prnlist[i]][0] = true;
			}
		}
	}
	if (pMedia2->ReWeight(nNormal, 1))
	{
		for (int i = 0; i < pMedia2->_nOmc; i++)
		{
			if (pMedia2->GetOMCWei(i) > 0.99)
			{
				iValid[1]++;
			}
			else
			{
				bFlag[prnlist[i]][1] = true;
			}
		}
	}
	if ((iValid[0] >= iValid[1]) && iValid[0] >= MINAMB && (pMedia1->_nOmc - iValid[0]) < 3) //
	{
		_dRcvUPD[iStaNo] = pMedia1->GetAverage();
		iMediaNum = 0;
		for(int prn_1 = _iPrnNum[0]; prn_1 < _iPrnNum[1]; prn_1++)
		{
			_PInfo[iStaNo].bUPDExclude[prn_1][_iType] = bFlag[prn_1][0];
		}
	}
	else if ((iValid[1] > iValid[0]) && iValid[1] >= MINAMB && (pMedia2->_nOmc - iValid[1]) < 3) //
	{
		_dRcvUPD[iStaNo] = pMedia2->GetAverage();
		iMediaNum = 1;
		for(int prn_1 = _iPrnNum[0]; prn_1 < _iPrnNum[1]; prn_1++)
		{
			_PInfo[iStaNo].bUPDExclude[prn_1][_iType] = bFlag[prn_1][1];
		}
	}
	else
	{
		_bSiteSlip[iStaNo] = true;
		return false;
	}							// receiver UPD is initialized 
	_bSiteSlip[iStaNo] = false;
	return true;
}
void UPDModel::GetNextUPD()
{
	//============================Select Max SV station for reference station
	float stdcoutoff[2] = {0.02, 0.20};
	MedialDetection* pMedia1 = new MedialDetection(MAXSAT, stdcoutoff);			// gross error detection for mean UPD [0.0 1.0]                    
	MedialDetection* pMedia2 = new MedialDetection(MAXSAT, stdcoutoff);			// gross error detection for mean UPD [-0.5 0.5] 

	while(1)
	{
		int iNext = -1, prn_1 = 0;
		int prnlist[MAXSAT] = {0};
		int nPrn = 0;
		FindStation(iNext, false);
		if(iNext == -1)
		{
			break;
		}
		pMedia1->Initialize();
		pMedia2->Initialize();
		_bIniSite[iNext] = true;
		double UPDTmp[2], SatUPDTmp[MAXSAT][2] = {0.0};
		printf("_iSite %d \n", iNext);
		fprintf(_ftOut,"  %s ,", _PInfo[iNext].Station[iNext].Name);
		for(prn_1 = _iPrnNum[0]; prn_1 < _iPrnNum[1]; prn_1++)
		{
			if(_PInfo[iNext].UDFcbFlag[prn_1][_iType] && _bIniSat[prn_1] > 0)
			{
				UPDTmp[0] = UPDTmp[1] = _PInfo[iNext].UDFcb[prn_1][_iType] - _dSatUPD[prn_1];
				SatUPDTmp[prn_1][0] = SatUPDTmp[prn_1][1] = _PInfo[iNext].UDFcb[prn_1][_iType];
				while (UPDTmp[0] < 0.0)
				{
					SatUPDTmp[prn_1][0] += 1.0;
					UPDTmp[0] += 1.0;										// [0.0 1.0]                                                       
				}
				while (UPDTmp[0] > 1.0)
				{
					SatUPDTmp[prn_1][0] -= 1.0;
					UPDTmp[0] -= 1.0;										// [0.0 1.0]                                                       
				}
				while (fabs(UPDTmp[1]) > 0.5)
				{
					SatUPDTmp[prn_1][1] -= sign(UPDTmp[1]);
					UPDTmp[1] -= sign(UPDTmp[1]);							// [-0.5 0.5]                                                      
				}
#if (UPDCONSISTENCE == 2)
				if(prn_1 < NSATGPS + 16)
				{
					continue;
				}
#else if (UPDCONSISTENCE == 1)
				if(prn_1 >= NSATGPS + 16)
				{
					continue;
				}
#endif
				pMedia1->UpdateList(UPDTmp[0]);
				pMedia2->UpdateList(UPDTmp[1]);
				prnlist[nPrn] = prn_1;
				nPrn++;
				printf("%8.3f %8.3f %8.3f\n", UPDTmp[0], UPDTmp[1], _PInfo[iNext].UDFcb[prn_1][_iType]);
			}
		}
		int iMediaNum = 0;
		//----------------------check which UPD should be specified to receiver: iSNext between [0.0 1.0] and [-0.5 0.5]-----------------------
		if(ExamineConsistence(pMedia1, pMedia2, iNext, iMediaNum, prnlist))
		{
			//--------------------------------------------------get new SV UPDs from site: iSNext--------------------------------------------------
			for(prn_1 = _iPrnNum[0]; prn_1 < _iPrnNum[1]; prn_1++)
			{
				if(_PInfo[iNext].UDFcbFlag[prn_1][_iType])
				{
					if(_iIteration == 1)
					{
						_PInfo[iNext].UDFcb[prn_1][_iType] = SatUPDTmp[prn_1][iMediaNum] - _dRcvUPD[iNext];
					}
					else
					{
						_PInfo[iNext].UDFcb[prn_1][_iType] = _PInfo[iNext].UDFcb[prn_1][_iType] - _dRcvUPD[iNext]; 
					}
					if (_bIniSat[prn_1])
					{
						while(_PInfo[iNext].UDFcb[prn_1][_iType] - _dSatUPD[prn_1] > 0.5)
						{
							_PInfo[iNext].UDFcb[prn_1][_iType] -= 1.0;
						}
						while(_PInfo[iNext].UDFcb[prn_1][_iType] - _dSatUPD[prn_1] < -0.5)
						{
							_PInfo[iNext].UDFcb[prn_1][_iType] += 1.0;
						}
					}
					fprintf(_ftOut,"%7.4f ,", _PInfo[iNext].UDFcb[prn_1][_iType]);
					if (!_bIniSat[prn_1])		// new SV UPD generated by the new site     
					{
						double tmp = _PInfo[iNext].UDFcb[prn_1][_iType];
						_dSatUPD[prn_1] = tmp;
						_bIniSat[prn_1] = true;										// SV have already known initial UPDs                        
					}
				}
				else
				{
					fprintf(_ftOut,"        ,");
				}
			}
		}	
		fprintf(_ftOut,"\n");
	}
	delete pMedia1;
	delete pMedia2;
	return;
}
void UPDModel::GetAverageUPD()
{
	float stdcoutoff[2] = {0.02, 0.15};
	MedialDetection* pMedia1 = new MedialDetection(MAXSTATION, stdcoutoff);			// gross error detection for mean UPD [0.0 1.0]   
	int prn_1 = 0, StaNo = 0;
	double SatSTD[MAXSAT] = {0.0};
	fprintf(_ftOut,"AveUPD ,");
	for(prn_1 = _iPrnNum[0]; prn_1 < _iPrnNum[1]; prn_1++)
	{
		pMedia1->Initialize();
		int iSite = 0;
		for(StaNo = 0; StaNo < _nSite; StaNo++)
		{		
			if(!_bSiteSlip[StaNo] && _PInfo[StaNo].UDFcbFlag[prn_1][_iType] && !_PInfo[StaNo].bUPDExclude[prn_1][_iType])
			{
				pMedia1->UpdateList(_PInfo[StaNo].UDFcb[prn_1][_iType]);
				iSite++;
			}
		}
		if(iSite < 1)
		{
			fprintf(_ftOut,"        ,");
			continue;
		}
		int nNormal = 0;
		if (pMedia1->ReWeight(nNormal, 1))
		{
			_bIniSat[prn_1] = true;
			_dSatUPD[prn_1] = pMedia1->GetAverage();
			SatSTD[prn_1] = pMedia1->GetSTD();
			fprintf(_ftOut,"%7.4f ,", _dSatUPD[prn_1]);
		}
		else
		{
			fprintf(_ftOut,"        ,");
			printf("Prn %d is not Valid %d\n",prn_1, iSite);
			_bIniSat[prn_1] = false;
		}
	}
	fprintf(_ftOut,"\n");

	if(_iIteration > 1)
	{
		for(prn_1 = _iPrnNum[0]; prn_1 < _iPrnNum[1]; prn_1++)
		{
			if(!_bIniSat[prn_1])
			{
				continue;
			}
			if(prn_1 < 32)
			{
				fprintf(_ftOut, "G%02d %8.4f %8.4f\n", prn_1 + 1, _dSatUPD[prn_1], SatSTD[prn_1]);
			}
			else
			{
				fprintf(_ftOut, "C%02d %8.4f %8.4f\n", prn_1 - 31, _dSatUPD[prn_1], SatSTD[prn_1]);
			}
			for(StaNo = 0; StaNo < _nSite; StaNo++)
			{
				if(!_bSiteSlip[StaNo] && _PInfo[StaNo].UDFcbFlag[prn_1][_iType] && !_PInfo[StaNo].bUPDExclude[prn_1][_iType])
				{
					fprintf(_pUPD->_ftRes, "prn %d _iSite %d %6.3f\n", prn_1 + 1, StaNo, _PInfo[StaNo].UDFcb[prn_1][_iType] - _dSatUPD[prn_1]);
				}
			}
		}	
	}	
	delete pMedia1;
}
void UPDModel::InitialUPD(UPDInfo* UPDIf, GNSSInfo *PInfo, int StaNum, int Type, FILE *ft, int nSat[2])
{
	_ftOut = ft;
	_PInfo = PInfo;
	_pUPD  = UPDIf;
	_nSite = StaNum;
	_iType = Type;
	_iPrnNum[0] = nSat[0];
	_iPrnNum[1] = nSat[1];
	_iIteration = 0;
	memset(_bIniSat,   0x00, sizeof(_bIniSat));
	memset(_dSatUPD,   0x00, sizeof(_dSatUPD));
	memset(_bIniSite,  0x00, sizeof(_bIniSite));
	memset(_bSiteSlip, 0x00, sizeof(_bSiteSlip));
	memset(_dRcvUPD,   0x00, sizeof(_dRcvUPD));
	return;
}

void UPDModel::OutPutUPD(char chData[255])
{
#ifdef FCB_DETAIAL	
	fprintf(_ftOut,"---------------------------- %s ------------------------------\nSite   ,", chData);
	for(int prn_1 = _iPrnNum[0]; prn_1 < _iPrnNum[1]; prn_1++)
	{
		if(prn_1 <= 31)
		{
			fprintf(_ftOut,"  G%02d,   ",prn_1+1);
		}
		else
		{
			fprintf(_ftOut,"  C%02d,   ",prn_1-31);
		}
	}
	fprintf(_ftOut,"\n");
	for(int StaNo = 0; StaNo < _nSite; StaNo++)
	{
		fprintf(_ftOut,"%d_%s ,", _bSiteSlip[StaNo], _PInfo[StaNo].Station[StaNo].Name);
		for(int prn_1 = _iPrnNum[0]; prn_1 < _iPrnNum[1]; prn_1++)
		{
			if(_PInfo[StaNo].UDFcbFlag[prn_1][_iType])
			{
				fprintf(_ftOut,"%7.4f ,", _PInfo[StaNo].UDFcb[prn_1][_iType]);
			}
			else
			{
				fprintf(_ftOut,"    -   ,");
			}
		}
		fprintf(_ftOut,"\n");
	}
	fprintf(_ftOut,"---------------------------------------------------------------------\n");
#endif
}
void UPDModel::ComputeStatistical(UPDInfo* UPDIf, GNSSInfo *PInfo, int StaNum, int Type, FILE *ft, int nSat[2])
{
	InitialUPD(UPDIf, PInfo, StaNum, Type, ft, nSat);

	char chData[255];
	memset(chData, 0x00, sizeof(chData));
	sprintf(chData,"Before Minus Bias");
	OutPutUPD(chData);

	GetIniUPD();
	//=====================================
	_iIteration++;
	GetNextUPD();
	memset(chData, 0x00, sizeof(chData));
	sprintf(chData,"After Minus Bias1");
	OutPutUPD(chData);
	GetAverageUPD();

	memset(_bSiteSlip, 0x00, sizeof(_bSiteSlip));
	memset(_bIniSite, 0x00, sizeof(_bIniSite));

	//=====================================
	_iIteration++;
	GetNextUPD();
	memset(chData, 0x00, sizeof(chData));	
	sprintf(chData,"After Minus Bias2");
	OutPutUPD(chData);
	GetAverageUPD();
	return;
}