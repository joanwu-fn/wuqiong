#include "../Common/BaseFunction.h"

const int nDayUPD = 1;

int LoadUPD(char chFile[], double *UPD, bool *bUPD)
{
	FILE* ft;
	ft = fopen(chFile, "r");
	if(ft == NULL)
	{
		return 0;
	}
	char chLineTemp[1000];
	while (fgets(chLineTemp, sizeof(chLineTemp), ft))
	{
		char chPrn[4];
		double UPDValue[30] = {0};
		u2 Sys_Sys;
		sscanf(chLineTemp, "%s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
			chPrn, &UPDValue[0], &UPDValue[1], &UPDValue[2], &UPDValue[3], &UPDValue[4], &UPDValue[5], &UPDValue[6], &UPDValue[7],
			&UPDValue[8], &UPDValue[9], &UPDValue[10], &UPDValue[11], &UPDValue[12], &UPDValue[13], &UPDValue[14], &UPDValue[15],
			&UPDValue[16], &UPDValue[17], &UPDValue[18], &UPDValue[19], &UPDValue[20], &UPDValue[21], &UPDValue[22], &UPDValue[23],
			&UPDValue[24], &UPDValue[25], &UPDValue[26], &UPDValue[27], &UPDValue[28], &UPDValue[29]);
		if((JudgeSys(chLineTemp[0], SYS_ALL, &Sys_Sys)) == -1)
		{
			continue;
		}
		int prn = (int)str2num(chPrn,1,2);
		prn     = SatNo(prn, Sys_Sys);
		if(prn < 1 || prn > MAXSAT)
		{
			continue;
		}
		for(int i = 0; i < nDayUPD; i++)
		{
			UPD[(prn - 1) * nDayUPD + i] = UPDValue[i];
			if(fabs(UPDValue[i]) > 0.0)
			{
				bUPD[(prn - 1) * nDayUPD + i] = true;
			}
			else
			{
				bUPD[(prn - 1) * nDayUPD + i] = false;
			}
		}
	}
	fclose(ft);
	return 1;
}

void RePutUPD()
{
	double *EWLUPD, *WLUPD;
	bool   *bEWLUPD, *bWLUPD;
	EWLUPD = new double[MAXSAT * nDayUPD];
	WLUPD  = new double[MAXSAT * nDayUPD];
	bEWLUPD = new bool[MAXSAT * nDayUPD];
	bWLUPD  = new bool[MAXSAT * nDayUPD];
	memset(EWLUPD,  0x00, sizeof(double) * MAXSAT * nDayUPD);
	memset(WLUPD,   0x00, sizeof(double) * MAXSAT * nDayUPD);
	memset(bEWLUPD, 0x00, sizeof(bool) * MAXSAT * nDayUPD);
	memset(bWLUPD,  0x00, sizeof(bool) * MAXSAT * nDayUPD);

	LoadUPD("D:/EWL.UPD", EWLUPD, bEWLUPD);
	LoadUPD("D:/WL.UPD",  WLUPD,  bWLUPD);

	int year = 2019, month = 4, day = 20, doy = 110;

	for(int i = 0; i < 7; i++)
	{
		FILE *ftEWL, *ftWL, *ftNL;
		char chEWLUPD[255], chWLUPD[255], chNLUPD[255];
		sprintf(chEWLUPD,"D:/project/ZDPos/ambarc/EWL%03d.%02dupd", i + doy, year % 100);
		sprintf(chWLUPD, "D:/project/ZDPos/ambarc/WL%03d.%02dupd", i + doy, year % 100);
		sprintf(chNLUPD, "D:/project/ZDPos/ambarc/NL%03d.%02dupd", i + doy, year % 100);

		ftEWL = fopen(chEWLUPD, "w");
		ftWL  = fopen(chWLUPD, "w");
		ftNL  = fopen(chNLUPD, "w");

		int nNum[2] = {0};
		for(int prn_1 = 0; prn_1 < MAXSAT; prn_1++)
		{
			if(bEWLUPD[prn_1 * nDayUPD])
			{
				nNum[0]++;
			}
			if(bWLUPD[prn_1]) // * nDayUPD
			{
				nNum[1]++;
			}
		}
		for(int sod = 0; sod < 86400; sod += 30)
		{
			int hour = sod / 3600;
			int minute = (sod - hour * 3600) / 60;
			double sec = sod - hour * 3600 - minute * 60;
			fprintf(ftEWL,"* %04d %02d %02d %02d %02d %4.1f %3d\n", year, month, day + i, hour, minute, sec, nNum[0]);
			fprintf(ftWL,"* %04d %02d %02d %02d %02d %4.1f %3d\n", year, month, day + i, hour, minute, sec, nNum[1]);
			for(int prn_1 = 0; prn_1 < MAXSAT; prn_1++)
			{
				char chSat[4];
				if(prn_1 < NSATGPS)
				{
					sprintf(chSat,"G%02d", prn_1 + 1);
				}
				else if(prn_1 < NSATGPS + NSATBDS)
				{
					sprintf(chSat,"C%02d", prn_1 - NSATGPS + 1);
				}
				else if(prn_1 < NSATGPS + NSATBDS + NSATGAL)
				{
					sprintf(chSat,"E%02d", prn_1 - NSATGPS - NSATBDS+ 1);
				}
				if(bEWLUPD[prn_1 * nDayUPD])
				{
					fprintf(ftEWL,"%s %8.4f %8.3f\n", chSat, EWLUPD[prn_1 * nDayUPD], 0.02); 
				}
				if(bWLUPD[prn_1])//nDayUPD
				{
					fprintf(ftWL,"%s %8.4f %8.3f\n", chSat, WLUPD[prn_1 ], 0.02); //* nDayUPD
				}
			}
		}
		fclose(ftEWL);
		fclose(ftWL);
		fclose(ftNL);
	}
	delete[] EWLUPD;
	delete[] WLUPD;
	delete[] bEWLUPD;
	delete[] bWLUPD;
}
