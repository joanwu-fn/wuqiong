#include "../Common/SVEphemeris.h"
#include "../Common/BaseFunction.h"
#include "../Common/GNSSERROR.h"
#include "../Common/Common.h"

/* variance by ura ephemeris (ref [1] 20.3.3.3.1.1) --------------------------*/
double var_uraeph(int ura)
{
	const double ura_value[]={   
		2.4,3.4,4.85,6.85,9.65,13.65,24.0,48.0,96.0,192.0,384.0,768.0,1536.0,
		3072.0,6144.0
	};
	return ura<0||15<ura?6144.0:SQR(ura_value[ura]);
}

/* polynomial interpolation by Neville's algorithm ---------------------------*/
double interppol(const double *x, double *y, int n)
{
	int i,j;

	for (j=1;j<n;j++) {
		for (i=0;i<n-j;i++) {
			y[i]=(x[i+j]*y[i]-x[i]*y[i+1])/(x[i+j]-x[i]);
		}
	}
	return y[0];
}

Orb_Clk::Orb_Clk()
{
	m_nclk = 0;
	m_nsp3 = 0;
	m_type = 0;
}
Orb_Clk::~Orb_Clk()
{
	;
}
/* satellite antenna phase center offset ---------------------------------------
* compute satellite antenna phase center offset in ecef
* args   : gtime_t time       I   time (gpst)
*          double *rs         I   satellite position and velocity (ecef)
*                                 {x,y,z,vx,vy,vz} (m|m/s)
*          int    sat         I   satellite number
*          nav_t  *nav        I   navigation data
*          double *dant       I   satellite antenna phase center offset (ecef)
*                                 {dx,dy,dz} (m) (iono-free LC value)
* return : none
*-----------------------------------------------------------------------------*/
void Orb_Clk::satantoff(gtime_t time, const double *rs,  Pcv_t *pcv,double *dant)
{
	double ex[3],ey[3],ez[3],es[3],r[3],rsun[3],gmst,erpv[5]={0};
	double gamma,C1,C2,dant1,dant2;
	int i,j=0,k=1;
	u2 SYS,SysNum;
	char StrSys;

	/* sun position in ecef */
	sunmoonpos(gpst2utc(time),erpv,rsun,NULL,&gmst);

	/* unit vectors of satellite fixed coordinates */
	for (i=0;i<3;i++) r[i]=-rs[i];
	if (!normv3(r,ez)) return;
	for (i=0;i<3;i++) r[i]=rsun[i]-rs[i];
	if (!normv3(r,es)) return;
	cross3(ez,es,r);
	if (!normv3(r,ey)) return;
	cross3(ey,ez,ex);

	if(Prn2Sys(pcv->sat, SYS, SysNum, StrSys) <= 0 )
	{
		dant[0]=dant[1]=dant[2] = 0.0;
		return;
	}
	gamma=SQR(G_Lam[SysNum][k])/SQR(G_Lam[SysNum][j]);
	C1=gamma/(gamma-1.0);
	C2=-1.0 /(gamma-1.0);

	/* iono-free LC */
	for (i=0;i<3;i++) 
	{
		dant1=pcv->off[j][0]*ex[i]+pcv->off[j][1]*ey[i]+pcv->off[j][2]*ez[i];
		dant2=pcv->off[k][0]*ex[i]+pcv->off[k][1]*ey[i]+pcv->off[k][2]*ez[i];
		dant[i]=C1*dant1+C2*dant2;
	}
}
bool Orb_Clk::Get_Orb_Igs(gtime_t time, u2 prn, double SatXYZ[],double &ephRMS)
{
	int i,j,k,index;
	if (m_nsp3 < NMAX+1|| timediff(time,m_sp3[0].time) < -m_SampleRate[0] || timediff(time,m_sp3[m_nsp3 -1].time) > m_SampleRate[0]) 
	{
		return false;
	}
	/* binary search */
	for (i=0,j=m_nsp3-1; i<j;) 
	{
		k=(i+j)/2;
		if (timediff(m_sp3[k].time,time) < 0.0) 
			i=k+1; 
		else 
			j=k;
	}
	index= i<=0? 0:i-1;
	/* polynomial interpolation for orbit */
	i=index-(NMAX+1)/2;
	if (i < 0) 
	{
		i=0; 
	}
	else if (i+NMAX >= m_nsp3) 
	{
		i=m_nsp3 - NMAX-1;
	}
	double t[NMAX+1],p[3][NMAX+1];

	for (j=0;j<=NMAX;j++) {
		t[j]=timediff(m_sp3[i+j].time,time);
		if (norm(m_sp3[i+j].Orbit[prn-1],3) <= 0.0 || m_sp3[i+j].prn[prn-1] != prn) 
		{
			return false;
		}
	}
	for (j=0;j<=NMAX;j++) 
	{
		p[0][j]=m_sp3[i+j].Orbit[prn-1][0];
		p[1][j]=m_sp3[i+j].Orbit[prn-1][1];
		p[2][j]=m_sp3[i+j].Orbit[prn-1][2];
	}
	for (i=0;i<3;i++) {
		SatXYZ[i] = interppol(t,p[i],NMAX+1);
	}
	return true;
}
bool Orb_Clk::Get_Clk_Igs(gtime_t time, u2 prn, double &SatClk)
{
	int i,j,k,index;
	if (m_nclk < 2 || timediff(time,m_clk[0].time) < (-m_SampleRate[1] - 1)||
		timediff(time,m_clk[m_nclk-1].time) > (m_SampleRate[1] + 1))
	{
		return false;
	}
	/* binary search */
	for (i=0,j=m_nclk-1; i<j;) 
	{
		k=(i+j)/2;
		if (timediff(m_clk[k].time,time) < 0.0) 
			i=k+1; 
		else 
			j=k;
	}
	index= i<=0? 0:i-1;
	if(m_clk[index  ].prn[prn-1] != prn || m_clk[index + 1].prn[prn-1] != prn)
		return false;
	double t[2],c[2];
	/* linear interpolation for clock */
	t[0]=timediff(time, m_clk[index  ].time);
	t[1]=timediff(time, m_clk[index+1].time);
	c[0]=m_clk[index  ].SatClk[prn-1];
	c[1]=m_clk[index+1].SatClk[prn-1];

	if (t[0] <= -(m_SampleRate[1]+1.0)) {
		if ((SatClk = c[0]) == 0.0) 
			return false;
	}
	else if (t[1] >= (m_SampleRate[1] + 1.0)) {
		if ((SatClk = c[1])==0.0) 
			return false;
	}
	else if (fabs(c[0]*CLIGHT) > 1e-5 && fabs(c[1]*CLIGHT) > 1e-5 ) {
		SatClk=(c[1]*t[0]-c[0]*t[1])/(t[0]-t[1]);
		return true;
	}
	else
	{
		return false;
	}
	return false;
}
bool Orb_Clk::GetOrb_Clk_Igs(gtime_t time, u2 prn, double SatXYZ[], double &SatClk, double SatAzEl[], double &ephRMS)
{
	gtime_t NowTime;
	u2 i;
	double dant[3]={0};
	if(Get_Clk_Igs(time, prn, SatClk))
	{
		NowTime = timeadd(time, -SatClk);
		if(Get_Orb_Igs(NowTime, prn, SatXYZ, ephRMS))
		{
			NowTime = timeadd(NowTime, -0.001);
			double SatPre[3],Velocity[3];
			Get_Orb_Igs(NowTime, prn, SatPre, ephRMS);
			satantoff(NowTime, SatXYZ, &m_Pcv[prn-1], dant);
			for(i = 0; i < 3; i++)
			{
				Velocity[i] = (SatXYZ[i] - SatPre[i])/0.001;
				SatXYZ[i] += dant[i];
			}
			 SatClk = SatClk-2.0*dot(SatXYZ, Velocity, 3)/CLIGHT/CLIGHT;
			 ephRMS = 4.0;
			 return true;
		}
		else
		{
			return false;
		}
	}
	else
	{
		return false;
	}
}
void LoadIFCBCoef(char *chFile, int MJDN, double dSin[MAXSAT][4], double dCos[MAXSAT][4])
{
	FILE* ft = NULL;
	ft = fopen(chFile, "r");
	if(ft == NULL)
	{
		return;
	}
	char chLine[MAXLENGTH];
	while (fgets(chLine, MAXLENGTH, ft) != NULL)
	{
		int iMJDN = 0, iprn = 0;
		char chTmp[20];
		double coef[9];
		sscanf(chLine,"%d %s %d %s %lf %lf %lf %lf %lf %lf %lf %lf %lf", &iMJDN, chTmp, &iprn, chTmp, &coef[0],
			&coef[1], &coef[2], &coef[3], &coef[4], &coef[5], &coef[6], &coef[7], &coef[8]);
		if(iMJDN != MJDN)
		{
			continue;
		}
		int prnList[12] = {1, 3, 6, 8, 9, 10, 24, 25, 26, 27, 30, 32};
		int prn_1 = 0;
		if(iprn <= 12)
		{
			prn_1 = prnList[iprn-1] - 1;
		}
		else if(iprn != 14)
		{
			prn_1 = NSATGPS + iprn - 13;
		}
		
		for(int i = 0; i < 4; i++)
		{
			dCos[prn_1][i] = coef[1 + i] * 1000;
			dSin[prn_1][i] = coef[5 + i] * 1000;
		}
	}
	fclose(ft);
}

void Orb_Clk::Initialize(int MJDN)
{
	memset(_dCos, 0x00, sizeof(_dCos));
    memset(_dSin, 0x00, sizeof(_dSin));

	LoadIFCBCoef("D:/project/data/bias/IFCB/IFCBCoef58119.txt", MJDN, _dSin, _dCos);
	return;
}

u2 Orb_Clk::GetOrb_Clk(GNSSDATA &Obs, double RecClk, double SatXYZ[], double SatClk[], double SatAzEl[], double ephRMS[],double rb[3], double pos[3])
{
	u1 i,j;
	gtime_t Nowtime={0};
	u2 SateNumer=0;
	bool Flag = false;
	double dt;
	GNSSDATA ObsTemp = {0};
	double beta[MAXSAT] = {0}, mu = 0.0;

	GetBetaAngle(Obs.ObsTime, beta);
    TTime TMJDN = rtktime2gltime(Obs.ObsTime);

	ObsTemp.iSite = Obs.iSite;
	memcpy(&ObsTemp.ObsTime, &Obs.ObsTime, sizeof(gtime_t));
	int UseNum = 0;
	for(i = 0; i < Obs.NumSats; i++)
	{
		if(Obs.ObsData[i].prn == 0)
		{
			continue;
		}
		SatXYZ[i*3] = SatXYZ[i*3+1] =SatXYZ[i*3+2] = 0.0;
		SatClk[i] = ephRMS[i] = 0;
		SatAzEl[i*2] = SatAzEl[i*2+1]=0;
		for(j = 0; j < 3; j++)
		{
			if(Obs.ObsData[i].P[j] > 1000.0)
			{
				dt =  Obs.ObsData[i].P[j] / CLIGHT + RecClk;
				Nowtime = timeadd(Obs.ObsTime, -dt);
				break;
			}
		}
		if(j >= 3)
		{
			continue;
		}
		if(0 == m_type)
		{
			Flag = GetOrb_Clk_Brdc(Nowtime, Obs.ObsData[i].prn, &SatXYZ[i*3], SatClk[i], &SatAzEl[i*2], ephRMS[i]);
		}
		else
		{
			Flag = GetOrb_Clk_Igs(Nowtime, Obs.ObsData[i].prn, &SatXYZ[i*3], SatClk[i], &SatAzEl[i*2], ephRMS[i]);
		}
		if(Flag == true)
		{
			SateNumer++;
		}
		if(Flag == true && rb != NULL && pos != NULL)
		{
			u2 SysNum,SYS;
			char StrSys;
			double e[3];
			Prn2Sys(Obs.ObsData[i].prn, SYS, SysNum, StrSys);
			geodist(&SatXYZ[i*3], rb, e);
			satazel(pos, e, &SatAzEl[i*2]);
			/* -------------- 改正北斗伪距随高度角变化的码偏差 ---------------*/
			if(SYS_BDS == SYS)
			{
				double Bias;
				for(int k = 0; k < 3; k++)
				{
					if(Obs.ObsData[i].P[k] > PRECISION)
					{
						GetBDSCodeBiasVariations(Obs.ObsData[i].prn, SatAzEl[i*2+1], k, Bias);
						Obs.ObsData[i].P[k] += Bias;
					}										
				}
			}
			/* -------------- 改正北斗伪距随高度角变化的码偏差 ---------------*/
			if(SYS_GPS == SYS && fabs(Obs.ObsData[i].L[2]) > 0.0)
			{
				double BLH[3], coef = G_Lam[0][2] * (FREQ5 * FREQ5) / (FREQ1 * FREQ1 - FREQ5 * FREQ5);
				ecef2pos(&SatXYZ[i*3], BLH);
				mu = BLH[1] + (TMJDN.SoD * 15.0 / 3600.0) * D2R;

                double IFCB = GetIFCBCorr('G', Obs.ObsData[i].prn, TMJDN.MJDN, beta[Obs.ObsData[i].prn-1], mu);

				double IFCB1 = GetIFCBCorr(_dCos[Obs.ObsData[i].prn - 1], _dSin[Obs.ObsData[i].prn - 1], mu);

				Obs.ObsData[i].L[2] += IFCB * 1e-3 / coef;
			}
			else if(SYS_BDS == SYS && fabs(Obs.ObsData[i].L[2]) > 0.0)
			{
				double BLH[3], coef = G_Lam[0][2] * (FREQ5 * FREQ5) / (FREQ1 * FREQ1 - FREQ5 * FREQ5);
				ecef2pos(&SatXYZ[i*3], BLH);
				mu = BLH[1] + (TMJDN.SoD * 15.0 / 3600.0) * D2R;
				int prn = Obs.ObsData[i].prn - MAXPRNGPS;

// 				double IFCB = GetIFCBCorr('C', prn, TMJDN.MJDN, beta[Obs.ObsData[i].prn-1], mu);

				double IFCB = GetIFCBCorr(_dCos[Obs.ObsData[i].prn - 1], _dSin[Obs.ObsData[i].prn - 1], mu);

				Obs.ObsData[i].L[2] += IFCB * 1e-3 / coef;
			}
		}	
	}
	return SateNumer;
}
void Orb_Clk::GetBetaAngle(gtime_t time, double Beta[MAXSAT])
{
	double XSun[3];
	TTime tt;
	tt = rtktime2gltime(time);
	sunxyz(tt, XSun);
	memset(Beta, 0x00, sizeof(double) * MAXSAT);

	for(int prn_1 = 0; prn_1 < MAXSAT; prn_1++)
	{
		int prn = prn_1 + 1;
		if(m_eph[prn_1].prn != prn)
		{
			continue;
		}	
		double SatXYZ[3] = {0.0}, SatXYZ1[3] = {0.0}, ephRMS, SatClk, VSat[3] = {0.0};
		gtime_t gtmp = time;
		GetSat_Pos(gtmp, prn, SatXYZ, ephRMS, SatClk);

        gtmp = timeadd(time,1E-3);       
		GetSat_Pos(gtmp, prn, SatXYZ1, ephRMS, SatClk);

		for(int i = 0; i < 3; i++)
		{
			VSat[i] = (SatXYZ1[i] - SatXYZ[i]) / 1E-3;
		}
        double dist = 0;
		for (int j = 0; j < 3; j++)
		{
			dist += (SatXYZ[j]* SatXYZ[j]);
		}
		dist = sqrt(dist);

		/*-------------------------------------------------------apply correction of earth rotation----------------------------------------------------*/
		double rotateDegree = dist * G_OMGE[0] / CLIGHT;
		double x0 = SatXYZ[0];
		double x1 = SatXYZ[1];
		SatXYZ[0] = cos(rotateDegree) * x0 + sin(rotateDegree) * x1;
		SatXYZ[1] = -sin(rotateDegree) * x0 + cos(rotateDegree) * x1;
		double v0 = VSat[0];
		double v1 = VSat[1];
		VSat[0] = v0 - SatXYZ[1] * G_OMGE[0];
		VSat[1] = v1 + SatXYZ[0] * G_OMGE[0];
		/*-------------------------------------------------------apply correction of earth rotation----------------------------------------------------*/
		Coordinate RVec(SatXYZ[0], SatXYZ[1], SatXYZ[2]);
		Coordinate coorTmp(VSat[0], VSat[1], VSat[2]);
		Coordinate orbnorm = unit(coorTmp % RVec);
		Coordinate sunTmp = unit(Coordinate(XSun[0], XSun[1], XSun[2]));
		double tmp1[3] = { orbnorm.X, orbnorm.Y, orbnorm.Z };
		double tmp2[3] = { sunTmp.X, sunTmp.Y, sunTmp.Z };
		Beta[prn_1] = asin(Dot(3, tmp1, tmp2)) * R2D;
	}
	return;
}
bool Orb_Clk::GetOrb_Clk_Brdc(gtime_t time, u2 prn, double SatXYZ[], double &SatClk, double SatAzEl[], double &ephRMS)
{
	double detaToe,detaToc;
	double BrdcTime[MAXSYS] = {7201.0, 7201.0, 7201.0, 7201.0};
	u2 SYS,SysNum;
	char StrSys;
	if(m_eph[prn-1].prn != prn)
	{
		return false;
	}
	detaToe = timediff(time, m_eph[prn-1].TOE);
	detaToc = timediff(time, m_eph[prn-1].TOC);
	
	Prn2Sys(prn, SYS, SysNum, StrSys);
	if(fabs(detaToe) < BrdcTime[SysNum] && fabs(detaToc) < BrdcTime[SysNum])
	{
		SatClk =  m_eph[prn-1].clock_bias + m_eph[prn-1].clock_drift * detaToc + m_eph[prn-1].clock_driftrate * (detaToc*detaToc);
		if(SYS == SYS_BDS) //transfer the B3 to B1B2 ionoshpere-free combination
		{
// 			if(fabs(m_eph[prn-1].TGD[0]*CLIGHT) > 1e-5 && fabs(m_eph[prn-1].TGD[1]*CLIGHT) > 1e-5)
// 			{
// 				SatClk = SatClk + (2.48 * m_eph[prn-1].TGD[0] - 1.48 * m_eph[prn-1].TGD[1]);
// 			}
// 			else
			{
				SatClk = SatClk - (1.48 * G_B1B2Default[prn-1-NSATGPS] + G_B1B3Default[prn-1-NSATGPS])*1e-9;
			}
		}
		time = timeadd(time, -SatClk); // minus satellite clock
		GetSat_Pos(time, prn, SatXYZ, ephRMS, SatClk);
		return true;
	}
	else
	{
		return false;
	}	
}
void Orb_Clk::GetSat_Pos(gtime_t time, u2 prn, double SatXYZ[], double &ephRMS, double &SatClk)
{
	double tk,M,E,Ek,sinE,cosE,u,r,i,O,sin2u,cos2u,x,y,sinO,cosO,cosi,mu;
	double OMGE_TEMP,XYZ1[3],XYZ2[3];
	u2 sys,sat,Sysnum;
	char StrSys;

	sat = Prn2Sys(prn, sys, Sysnum, StrSys);
	if (m_eph[prn-1].sqrt_A <= 0.0 || sat <= 0) {
		SatXYZ[0] = SatXYZ[1] = SatXYZ[2] = ephRMS =0.0;
		return;
	}
	tk=timediff(time,m_eph[prn-1].TOE);

	
	mu        =  G_MU[Sysnum];
	OMGE_TEMP =  G_OMGE[Sysnum];

	M=m_eph[prn-1].M0+(sqrt(mu/(m_eph[prn-1].sqrt_A*m_eph[prn-1].sqrt_A*m_eph[prn-1].sqrt_A))+m_eph[prn-1].Delta_n)*tk;
	int icount = 0;
	for (E=M,sinE=Ek=0.0;fabs(E-Ek)>1E-12;) 
	{
		Ek=E; 
		sinE=sin(Ek); 
		E=M+m_eph[prn-1].e*sinE;
		icount++;
		if(icount > 50)
		{
			break;
		}
	}
	cosE=cos(E);
	u=atan2(sqrt(1.0-m_eph[prn-1].e*m_eph[prn-1].e)*sinE,cosE-m_eph[prn-1].e)+m_eph[prn-1].omega;
	r=m_eph[prn-1].sqrt_A*(1.0-m_eph[prn-1].e*cosE);
	i=m_eph[prn-1].i0+m_eph[prn-1].IDOT*tk;
	sin2u=sin(2.0*u); cos2u=cos(2.0*u);
	u+=m_eph[prn-1].Cus*sin2u+m_eph[prn-1].Cuc*cos2u;
	r+=m_eph[prn-1].Crs*sin2u+m_eph[prn-1].Crc*cos2u;
	i+=m_eph[prn-1].Cis*sin2u+m_eph[prn-1].Cic*cos2u;

	
	if(sys == SYS_BDS && sat <= 5)
	{
		O = m_eph[prn-1].OMEGA0 + m_eph[prn-1].OMEGADOT*tk - OMGE_TEMP*(m_eph[prn-1].toes); //BeiDou GEO
	    x=r*cos(u); y=r*sin(u); sinO=sin(O); cosO=cos(O); cosi=cos(i);

		XYZ1[0]=x*cosO-y*cosi*sinO;
		XYZ1[1]=x*sinO+y*cosi*cosO;
		XYZ1[2]=y*sin(i);

		//加两次旋转
		XYZ2[0] = XYZ1[0];
		XYZ2[1] =cos(-5.0*D2R)*XYZ1[1]+sin(-5.0*D2R)*XYZ1[2];
		XYZ2[2] =sin(5.0*D2R)*XYZ1[1]+cos(-5.0*D2R)*XYZ1[2]; 

		SatXYZ[0]=XYZ2[0]*cos(OMGE_TEMP*(tk))+XYZ2[1]*sin(OMGE_TEMP*(tk));
		SatXYZ[1]=-1.0*XYZ2[0]*sin(OMGE_TEMP*(tk))+XYZ2[1]*cos(OMGE_TEMP*(tk));
		SatXYZ[2]=XYZ2[2];
	}
	else
	{
		O=m_eph[prn-1].OMEGA0+(m_eph[prn-1].OMEGADOT-OMGE_TEMP)*tk-OMGE_TEMP*m_eph[prn-1].toes;  //?????
		x=r*cos(u); y=r*sin(u); sinO=sin(O); cosO=cos(O); cosi=cos(i);
		SatXYZ[0]=x*cosO-y*cosi*sinO;
		SatXYZ[1]=x*sinO+y*cosi*cosO;
		SatXYZ[2]=y*sin(i);
	}
	
	/* relativity correction */
	SatClk -= 2.0*sqrt(mu*m_eph[prn-1].sqrt_A)*m_eph[prn-1].e*sinE/SQR(CLIGHT);

	/* position and clock error variance */
	ephRMS = var_uraeph(m_eph[prn-1].URAindex);
}