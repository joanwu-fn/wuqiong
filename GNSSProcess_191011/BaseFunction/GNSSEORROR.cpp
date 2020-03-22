#include "../Common/GNSSERROR.h"
#define ERR_SAAS    0.3         /* saastamoinen model error std (m) */
#define REL_HUMI    0.7         /* relative humidity for saastamoinen model */
#define SQRT(x)     ((x)<0.0?0.0:sqrt(x))
#define VAR_NOTEC   SQR(30.0)   /* variance of no tec */
#define MIN_EL      0.0         /* min elevation angle (rad) */
#define MIN_HGT     -1000.0     /* min user height (m) */

const double chisqr[100]={      /* chi-sqr(n) (alpha=0.001) */
	10.8,13.8,16.3,18.5,20.5,22.5,24.3,26.1,27.9,29.6,
	31.3,32.9,34.5,36.1,37.7,39.3,40.8,42.3,43.8,45.3,
	46.8,48.3,49.7,51.2,52.6,54.1,55.5,56.9,58.3,59.7,
	61.1,62.5,63.9,65.2,66.6,68.0,69.3,70.7,72.1,73.4,
	74.7,76.0,77.3,78.6,80.0,81.3,82.6,84.0,85.4,86.7,
	88.0,89.3,90.6,91.9,93.3,94.7,96.0,97.4,98.7,100 ,
	101 ,102 ,103 ,104 ,105 ,107 ,108 ,109 ,110 ,112 ,
	113 ,114 ,115 ,116 ,118 ,119 ,120 ,122 ,123 ,125 ,
	126 ,127 ,128 ,129 ,131 ,132 ,133 ,134 ,135 ,137 ,
	138 ,139 ,140 ,142 ,143 ,144 ,145 ,147 ,148 ,149
};
const double IGSOCorrection[3][10]={{-0.55,-0.40,-0.34,-0.23,-0.15,-0.04,0.09,0.19,0.27,0.35},
	{-0.71,-0.36,-0.33,-0.19,-0.14,-0.03,0.08,0.17,0.24,0.33},
	{-0.27,-0.23,-0.21,-0.15,-0.11,-0.04,0.05,0.14,0.19,0.32}};
const double MEOCorrection[3][10]={ {-0.47,-0.38,-0.32,-0.23,-0.11,0.06,0.34,0.69,0.97,1.05},
	{-0.40,-0.31,-0.26,-0.18,-0.06,0.09,0.28,0.48,0.64,0.69},
	{-0.22,-0.15,-0.13,-0.10,-0.04,0.05,0.14,0.27,0.36,0.47}};

/* single-difference noise variance ------------------------------------------*/
double SD_var(double var, double el)
{							
	double sinel=sin(el);	
	return 2.0*(var+var/sinel/sinel);
}							
/* noise variance of LC (m) 0.40--------------------------------------------------*/
double var_LC(double i, double j, double k, double sig,int sys)
{							
	double f1,f2,f5;		
	if(sys == SYS_GPS)		
	{						
		f1=FREQ1;f2=FREQ2;f5=FREQ5;
	}
	else if(sys == SYS_BDS)
	{
		f1=FREB1;f2=FREB2;f5=FREB5;
	}
	return (SQR(i*f1)+SQR(j*f2)+SQR(k*f5))/SQR(i*f1+j*f2+k*f5)*SQR(sig);
}
/* pseudorange measurement error variance ------------------------------------*/
double varerr(double el, int type, u2 sys)
{
	double fact=1.0,varr=1.0;
	double a,b;

	if (type==1)  // pseu
	{
		if(sys == SYS_GPS)
			fact *= 100.0;
		else
			fact *= 150;
	}
	a=fact*0.003;
	b=fact*0.003;
	return a*a+b*b/(sin(el)*sin(el));
}

void dops(int ns, const double *azel, double elmin, double *dop)
{
	double H[4*MAXSAT],Q[16],cosel,sinel;
	int i,n;

	for (i=0;i<4;i++) dop[i]=0.0;
	for (i=n=0;i<ns&&i<MAXSAT;i++) {
		if (azel[1+i*2]<elmin||azel[1+i*2]<=0.0) continue;
		cosel=cos(azel[1+i*2]);
		sinel=sin(azel[1+i*2]);
		H[  4*n]=cosel*sin(azel[i*2]);
		H[1+4*n]=cosel*cos(azel[i*2]);
		H[2+4*n]=sinel;
		H[3+4*n++]=1.0;
	}
	if (n<4) return;

	matmul("NT",4,4,n,1.0,H,H,0.0,Q);
	if (!matinv(Q,4)) {
		dop[0]=SQRT(Q[0]+Q[5]+Q[10]+Q[15]); /* GDOP */
		dop[1]=SQRT(Q[0]+Q[5]+Q[10]);       /* PDOP */
		dop[2]=SQRT(Q[0]+Q[5]);             /* HDOP */
		dop[3]=SQRT(Q[10]);                 /* VDOP */
	}
}
/* validate solution ---------------------------------------------------------*/
int valsol(const double *azel, const double *v, int nv, int nx)
{
	double vv;

	/* chi-square validation of residuals */
	vv = 3 * dot(v,v,nv);

	if (nv > nx && vv > chisqr[nv-nx-1]) 
	{
		return 0;
	}
	return 1;
}
/* troposphere mapping function (NMF) ------------------------------------------
* compute tropospheric mapping function by NMF
* args   : gtime_t t        I   time
*          double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          double *mapfw    IO  wet mapping function (NULL: not output)
* return : dry mapping function
* note   : see ref [5]
*          original JGR paper of [5] has bugs in eq.(4) and (5). the corrected
*          paper is obtained from:
*          ftp://web.haystack.edu/pub/aen/nmf/NMF_JGR.pdf
*-----------------------------------------------------------------------------*/
double interpc(const double coef[], double lat)
{
	int i=(int)(lat/15.0);
	if (i<1) return coef[0]; else if (i>4) return coef[4];
	return coef[i-1]*(1.0-lat/15.0+i)+coef[i]*(lat/15.0-i);
}
double mapf(double el, double a, double b, double c)
{
	double sinel=sin(el);
	return (1.0+a/(1.0+b/(1.0+c)))/(sinel+(a/(sinel+b/(sinel+c))));
}

double tropmapf(gtime_t time, const double pos[], const double azel[],
				double *mapfw)
{
	/* ref [5] table 3 */
	/* hydro-ave-a,b,c, hydro-amp-a,b,c, wet-a,b,c at latitude 15,30,45,60,75 */
	const double coef[][5]={
		{ 1.2769934E-3, 1.2683230E-3, 1.2465397E-3, 1.2196049E-3, 1.2045996E-3},
		{ 2.9153695E-3, 2.9152299E-3, 2.9288445E-3, 2.9022565E-3, 2.9024912E-3},
		{ 62.610505E-3, 62.837393E-3, 63.721774E-3, 63.824265E-3, 64.258455E-3},

		{ 0.0000000E-0, 1.2709626E-5, 2.6523662E-5, 3.4000452E-5, 4.1202191E-5},
		{ 0.0000000E-0, 2.1414979E-5, 3.0160779E-5, 7.2562722E-5, 11.723375E-5},
		{ 0.0000000E-0, 9.0128400E-5, 4.3497037E-5, 84.795348E-5, 170.37206E-5},

		{ 5.8021897E-4, 5.6794847E-4, 5.8118019E-4, 5.9727542E-4, 6.1641693E-4},
		{ 1.4275268E-3, 1.5138625E-3, 1.4572752E-3, 1.5007428E-3, 1.7599082E-3},
		{ 4.3472961E-2, 4.6729510E-2, 4.3908931E-2, 4.4626982E-2, 5.4736038E-2}
	};
	const double aht[]={ 2.53E-5, 5.49E-3, 1.14E-3}; /* height correction */

	double y,cosy,ah[3],aw[3],dm;
	double az=azel[0],el=azel[1],lat=pos[0]*R2D,lon=pos[1]*R2D,hgt=pos[2];
	int i;

	if (el<=0.0) {
		if (mapfw) *mapfw=0.0;
		return 0.0;
	}
	/* year from doy 28, added half a year for southern latitudes */
	y=(time2doy(time)-28.0)/365.25+(lat<0.0?0.5:0.0);

	cosy=cos(2.0*PI*y);
	lat=fabs(lat);

	for (i=0;i<3;i++) {
		ah[i]=interpc(coef[i  ],lat)-interpc(coef[i+3],lat)*cosy;
		aw[i]=interpc(coef[i+6],lat);
	}
	/* ellipsoidal height is used instead of height above sea level */
	dm=(1.0/sin(el)-mapf(el,aht[0],aht[1],aht[2]))*hgt/1E3;

	if (mapfw) *mapfw=mapf(el,aw[0],aw[1],aw[2]);

	return mapf(el,ah[0],ah[1],ah[2])+dm;
}

double tropmodel(gtime_t time, const double pos[], const double azel[],
				 double humi)
{
	const double temp0=15.0; /* temparature at sea level */
	double hgt,pres,temp,e,z,trph,trpw;

	if (pos[2]<-100.0||1E4<pos[2]||azel[1]<=0) return 0.0;

	/* standard atmosphere */
	hgt=pos[2]<0.0?0.0:pos[2];

	pres=1013.25*pow(1.0-2.2557E-5*hgt,5.2568);
	temp=temp0-6.5E-3*hgt+273.16;
	e=6.108*humi*exp((17.15*temp-4684.0)/(temp-38.45));

	/* saastamoninen model */
	z=PI/2.0-azel[1];
	trph=0.0022768*pres/(1.0-0.00266*cos(2.0*pos[0])-0.00028*hgt/1E3)/cos(z);
	trpw=0.002277*(1255.0/temp+0.05)*e/cos(z);
	return trph+trpw;
}
/* data index (i:lat,j:lon,k:hgt) --------------------------------------------*/
int dataindex(int i, int j, int k, const int *ndata)
{
	if (i<0||ndata[0]<=i||j<0||ndata[1]<=j||k<0||ndata[2]<=k) return -1;
	return i+ndata[0]*(j+ndata[1]*k);
}
/* ionosphere model ------------------------------------------------------------
* compute ionospheric delay by broadcast ionosphere model (klobuchar model)
* args   : gtime_t t        I   time (gpst)
*          double *ion      I   iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3}
*          double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
* return : ionospheric delay (L1) (m)
*-----------------------------------------------------------------------------*/
double ionmodel(gtime_t t, const double ion[], const double pos[],const double azel[])
{
	const double ion_default[]={ /* 2004/1/1 */
		0.1118E-07,-0.7451E-08,-0.5961E-07, 0.1192E-06,
		0.1167E+06,-0.2294E+06,-0.1311E+06, 0.1049E+07
	};
	double tt,f,psi,phi,lam,amp,per,x;
	int week;

	if (pos[2]<-1E3||azel[1]<=0) return 0.0;
	if (norm(ion,8)<=0.0) 
		ion=ion_default;

	/* earth centered angle (semi-circle) */
	psi=0.0137/(azel[1]/PI+0.11)-0.022;

	/* subionospheric latitude/longitude (semi-circle) */
	phi=pos[0]/PI+psi*cos(azel[0]);
	if      (phi> 0.416) phi= 0.416;
	else if (phi<-0.416) phi=-0.416;
	lam=pos[1]/PI+psi*sin(azel[0])/cos(phi*PI);

	/* geomagnetic latitude (semi-circle) */
	phi+=0.064*cos((lam-1.617)*PI);

	/* local time (s) */
	tt=43200.0*lam+time2gpst(t,&week);
	tt-=floor(tt/86400.0)*86400.0; /* 0<=tt<86400 */

	/* slant factor */
	f=1.0+16.0*pow(0.53-azel[1]/PI,3.0);

	/* ionospheric delay */
	amp=ion[0]+phi*(ion[1]+phi*(ion[2]+phi*ion[3]));
	per=ion[4]+phi*(ion[5]+phi*(ion[6]+phi*ion[7]));
	amp=amp<    0.0?    0.0:amp;
	per=per<72000.0?72000.0:per;
	x=2.0*PI*(tt-50400.0)/per;

	return CLIGHT*f*(fabs(x)<1.57?5E-9+amp*(1.0+x*x*(-0.5+x*x/24.0)):5E-9);
}
/* ionosphere mapping function -------------------------------------------------
* compute ionospheric delay mapping function by single layer model
* args   : double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
* return : ionospheric mapping function
*-----------------------------------------------------------------------------*/
extern double ionmapf(const double *pos, const double *azel)
{
	if (pos[2]>=HION) return 1.0;
	return 1.0/cos(asin((RE_WGS84+pos[2])/(RE_WGS84+HION)*sin(PI/2.0-azel[1])));
}
/* ionospheric pierce point position -------------------------------------------
* compute ionospheric pierce point (ipp) position and slant factor
* args   : double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          double re        I   earth radius (km)
*          double hion      I   altitude of ionosphere (km)
*          double *posp     O   pierce point position {lat,lon,h} (rad,m)
* return : slant factor
* notes  : see ref [2], only valid on the earth surface
*          fixing bug on ref [2] A.4.4.10.1 A-22,23
*-----------------------------------------------------------------------------*/
double ionppp(const double *pos, const double *azel, double re,
					 double hion, double *posp)
{
	double cosaz,rp,ap,sinap,tanap;

	rp=re/(re+hion)*cos(azel[1]);
	ap=PI/2.0-azel[1]-asin(rp);
	sinap=sin(ap);
	tanap=tan(ap);
	cosaz=cos(azel[0]);
	posp[0]=asin(sin(pos[0])*cos(ap)+cos(pos[0])*sinap*cosaz);

	if ((pos[0]> 70.0*D2R&& tanap*cosaz>tan(PI/2.0-pos[0]))||
		(pos[0]<-70.0*D2R&&-tanap*cosaz>tan(PI/2.0+pos[0]))) {
			posp[1]=pos[1]+PI-asin(sinap*sin(azel[0])/cos(posp[0]));
	}
	else {
		posp[1]=pos[1]+asin(sinap*sin(azel[0])/cos(posp[0]));
	}
	return 1.0/sqrt(1.0-rp*rp);
}
/* interpolate tec grid data -------------------------------------------------*/
int interptec(const tec_t *tec, int k, const double *posp, double *value,
					 double *rms)
{
	double dlat,dlon,a,b,d[4]={0},r[4]={0};
	int i,j,n,index;

	if (tec->lats[2]==0.0||tec->lons[2]==0.0) return 0;

	dlat=posp[0]*R2D-tec->lats[0];
	dlon=posp[1]*R2D-tec->lons[0];
	if (tec->lons[2]>0.0) dlon-=floor( dlon/360)*360.0; /*  0<=dlon<360 */
	else                  dlon+=floor(-dlon/360)*360.0; /* -360<dlon<=0 */

	a=dlat/tec->lats[2];
	b=dlon/tec->lons[2];
	i=(int)floor(a); a-=i;
	j=(int)floor(b); b-=j;

	/* get gridded tec data */
	for (n=0;n<4;n++) {
		if ((index=dataindex(i+(n%2),j+(n<2?0:1),k,tec->ndata))<0)
		{
			continue;
		}
		d[n]=tec->data[index];
		r[n]=tec->rms [index];
	}
	if (d[0]>0.0&&d[1]>0.0&&d[2]>0.0&&d[3]>0.0) 
	{
		/* bilinear interpolation (inside of grid) */
		*value=(1.0-a)*(1.0-b)*d[0]+a*(1.0-b)*d[1]+(1.0-a)*b*d[2]+a*b*d[3];
		*rms  =(1.0-a)*(1.0-b)*r[0]+a*(1.0-b)*r[1]+(1.0-a)*b*r[2]+a*b*r[3];
	}
	/* nearest-neighbour extrapolation (outside of grid) */
	else if (a<=0.5&&b<=0.5&&d[0]>0.0) {*value=d[0]; *rms=r[0];}
	else if (a> 0.5&&b<=0.5&&d[1]>0.0) {*value=d[1]; *rms=r[1];}
	else if (a<=0.5&&b> 0.5&&d[2]>0.0) {*value=d[2]; *rms=r[2];}
	else if (a> 0.5&&b> 0.5&&d[3]>0.0) {*value=d[3]; *rms=r[3];}
	else return 0;

	return 1;
}
/* ionosphere delay by tec grid data -----------------------------------------*/
int iondelay(gtime_t time, const tec_t *tec, const double *pos,
					const double *azel, int opt, double *delay, double *var)
{
	const double fact=40.30E16/FREQ1/FREQ1; /* tecu->L1 iono (m) */
	double fs,posp[3]={0},vtec,rms,hion,rp;
	int i;

	*delay=*var=0.0;

	for (i=0;i<tec->ndata[2];i++) 
	{ /* for a layer */
		hion=tec->hgts[0]+tec->hgts[2]*i;
		/* ionospheric pierce point position */
		fs=ionppp(pos,azel,tec->rb,hion,posp);

		if (opt&2) 
		{
			/* modified single layer mapping function (M-SLM) ref [2] */
			rp=tec->rb/(tec->rb+hion)*sin(0.9782*(PI/2.0-azel[1]));
			fs=1.0/sqrt(1.0-rp*rp);
		}
		if (opt&1) 
		{
			/* earth rotation correction (sun-fixed coordinate) */
			posp[1]+=2.0*PI*timediff(time,tec->time)/86400.0;
		}
		/* interpolate tec grid data */
		if (!interptec(tec,i,posp,&vtec,&rms)) 
		{
			return 0;
		}

		*delay+=fact*fs*vtec;
		*var+=fact*fact*fs*fs*rms*rms;
	}
	return 1;
}
/* ionosphere model by tec grid data -------------------------------------------
* compute ionospheric delay by tec grid data
* args   : gtime_t time     I   time (gpst)
*          nav_t  *nav      I   navigation data
*          double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          int    opt       I   model option
*                                bit0: 0:earth-fixed,1:sun-fixed
*                                bit1: 0:single-layer,1:modified single-layer
*          double *delay    O   ionospheric delay (L1) (m)
*          double *var      O   ionospheric dealy (L1) variance (m^2)
* return : status (1:ok,0:error)
* notes  : before calling the function, read tec grid data by calling readtec()
*          return ok with delay=0 and var=VAR_NOTEC if el<MIN_EL or h<MIN_HGT
*-----------------------------------------------------------------------------*/
int iontec(gtime_t time, const double *pos, tec_t *tec, int nt,
				  const double *azel, int opt, double *delay, double *var)
{
	double dels[2],vars[2],a,tt;
	int i,stat[2];

	if (azel[1] < MIN_EL || pos[2] < MIN_HGT) 
	{
		*delay=0.0;
		*var=VAR_NOTEC;
		return 1;
	}
	
	for (i = 0; i < nt; i++) 
	{
		if (timediff(tec[i].time,time)>0.0) 
		{
			break;
		}
	}
	if (i == 0 || i >= nt) 
	{
		return 0;
	}
	if ((tt=timediff(tec[i].time, tec[i-1].time))==0.0) 
	{
		return 0;
	}
	stat[0]=iondelay(time, &tec[i-1],pos,azel,opt,dels  ,vars  );
	stat[1]=iondelay(time, &tec[i]  ,pos,azel,opt,dels+1,vars+1);

	if (!stat[0]&&!stat[1]) 
	{
		return 0;
	}
	if (stat[0]&&stat[1]) { 
		a=timediff(time,tec[i-1].time)/tt;
		*delay=dels[0]*(1.0-a)+dels[1]*a;
		*var  =vars[0]*(1.0-a)+vars[1]*a;
	}
	else if (stat[0]) 
	{ 
		*delay=dels[0];
		*var  =vars[0];
	}
	else 
	{
		*delay=dels[1];
		*var  =vars[1];
	}
	
	return 1;
}
/* precise tropspheric model -------------------------------------------------*/
double prectrop(gtime_t time, double pos[], double azel[],
				double x[], double dtdx[], double &var)
{
	double m_w=0.0,cotz,grad_n,grad_e;

	/* wet mapping function */
	tropmapf(time,pos,azel,&m_w);

	if (0 && azel[1]>0.0) {

		/* m_w=m_0+m_0*cot(el)*(Gn*cos(az)+Ge*sin(az)): ref [6] */
		cotz=1.0/tan(azel[1]);
		grad_n=m_w*cotz*cos(azel[0]);
		grad_e=m_w*cotz*sin(azel[0]);
		m_w+=grad_n*x[1]+grad_e*x[2];
		dtdx[1]=grad_n*x[0];
		dtdx[2]=grad_e*x[0];
	}
	else 
		dtdx[1]=dtdx[2]=0.0;
	dtdx[0]=m_w;
	return m_w*x[0];
}

void GetTropDelay(gtime_t time,double pos[], double azel[],double &Trop, double &Var)
{
	double zhd,zazel[]={0.0,90.0*D2R};
	//Trop =tropmodel(time,pos,azel,REL_HUMI);
	//Var=SQR(ERR_SAAS/(sin(azel[1])+0.1));

	zhd=tropmodel(time,pos,zazel,0.0);
	Trop = tropmapf(time, pos, azel, NULL)*zhd;
	return;
}
double liner(double x[2],double y[2],double z)
{
	double value = y[0] + (y[1] -y[0])*(z-x[0])/(x[1]-x[0]);
	return value;
}
inline double GetBDSIGSO_Cbias(int iFre, double Elevation)
{
	const double B1[10] = {-0.55, -0.40, -0.34, -0.23, -0.15, -0.04, 0.09, 0.19, 0.27, 0.35};
	const double B2[10] = {-0.71, -0.36, -0.33, -0.19, -0.14, -0.03, 0.08, 0.17, 0.24, 0.33};
	const double B3[10] = {-0.27, -0.23, -0.21, -0.15, -0.11, -0.04, 0.05, 0.14, 0.19, 0.32};
	int ind = 0;
	double k = 0.0, b = 0.0;
	Elevation *= R2D;
	if (Elevation < 0 || Elevation > 90)
	{
		return 0;
	}
	ind = floor(Elevation / 10);
	if (1 == iFre)
	{
		k = (B1[ind + 1] - B1[ind]) / 10;
		b = B1[ind] - k * ind * 10;
	}
	else if (2 == iFre)
	{
		k = (B2[ind + 1] - B2[ind]) / 10;
		b = B2[ind] - k * ind * 10;
	}
	else
	{
		k = (B3[ind + 1] - B3[ind]) / 10;
		b = B3[ind] - k * ind * 10;
	}
	return (k * Elevation + b);
}
inline double GetBDSMEO_Cbias(int iFre, double Elevation)
{
	const double B1[10] = {-0.47, -0.38, -0.32, -0.23, -0.11, 0.06, 0.34, 0.69, 0.97, 1.05};
	const double B2[10] = {-0.40, -0.31, -0.26, -0.18, -0.06, 0.09, 0.28, 0.48, 0.64, 0.69};
	const double B3[10] = {-0.22, -0.15, -0.13, -0.10, -0.04, 0.05, 0.14, 0.27, 0.36, 0.47};
	int ind = 0;
	double k = 0.0, b = 0.0;

	Elevation *= R2D;
	if (Elevation < 0 || Elevation > 90)
	{
		return 0;
	}
	ind = floor(Elevation / 10);
	if (1 == iFre)
	{
		k = (B1[ind + 1] - B1[ind]) / 10;
		b = B1[ind] - k * ind * 10;
	}
	else if (2 == iFre)
	{
		k = (B2[ind + 1] - B2[ind]) / 10;
		b = B2[ind] - k * ind * 10;
	}
	else
	{
		k = (B3[ind + 1] - B3[ind]) / 10;
		b = B3[ind] - k * ind * 10;
	}
	return /*-1**/(k * Elevation + b);
}
void GetBDSCodeBiasVariations(int prn, double Elevation, u2 Fre, double &bias)
{
	bias = 0.0;
	u2 SYS,SysNum;
	u2 Sat;
	char StrSys;
	Sat = Prn2Sys(prn, SYS, SysNum, StrSys);
	if (SYS == SYS_BDS)
	{
		if (Sat > 10 && Sat < 17 && Sat != 13 && Sat != 16)
		{
			bias = GetBDSMEO_Cbias(Fre + 1, Elevation);
		}
		else if (Sat <= 10 || Sat == 13 || Sat == 16)
		{
			bias = GetBDSIGSO_Cbias(Fre + 1, Elevation);
		}
	}
	return;
}
void GetBDSBias(double Satel, u2 Prn, u2 Fre, double &Bias)
{
	u2 SYS,SysNum;
	u2 Sat;
	char StrSys;

	Sat = Prn2Sys(Prn, SYS, SysNum, StrSys);
	if(SYS != SYS_BDS || Satel < 1e-5)
	{
		Bias = 0.0;
		return;
	}
	else 
	{
		double IGSO[3][4] = {{-0.590,1.624,-0.645,-0.281},{-0.257,0.995,-0.381, -0.278},{-0.102,0.748,-0.307, -0.248}};
		double MEO[3][4] = {{-0.946,2.158,-0.642, -0.465},{-0.598,1.635,-0.556, -0.380},{-0.177,0.652,-0.178, -0.252}};
		if(Sat <= 10 || Sat == 15) //GEO or IGSO
		{
			Bias =  IGSO[Fre][0]*Satel+IGSO[Fre][1]*Satel*Satel+IGSO[Fre][2]*Satel*Satel*Satel;//+IGSO[fre][3];
	     }
		else  //MEO
		{
			Bias =  MEO[Fre][0]*Satel+MEO[Fre][1]*Satel*Satel+MEO[Fre][2]*Satel*Satel*Satel;//+MEO[fre][3];
		}
		return;
	}
	return;
}
/* geometric distance ----------------------------------------------------------
* compute geometric distance and receiver-to-satellite unit vector
* args   : double *rs       I   satellilte position (ecef at transmission) (m)
*          double *rr       I   receiver position (ecef at reception) (m)
*          double *e        O   line-of-sight vector (ecef)
* return : geometric distance (m) (0>:error/no satellite position)
* notes  : distance includes sagnac effect correction
*-----------------------------------------------------------------------------*/
double geodist(const double *rs, const double *rr, double *e)
{
	double r;
	int i;

	if (norm(rs,3) < RE_WGS84) return -1.0;
	for (i=0;i<3;i++) e[i]=rs[i]-rr[i];
	r=norm(e,3);
	for (i=0;i<3;i++) e[i]/=r;
	return r+G_OMGE[0]*(rs[0]*rr[1]-rs[1]*rr[0])/CLIGHT;
}

/* satellite azimuth/elevation angle -------------------------------------------
* compute satellite azimuth/elevation angle
* args   : double *pos      I   geodetic position {lat,lon,h} (rad,m)
*          double *e        I   receiver-to-satellilte unit vevtor (ecef)
*          double *azel     IO  azimuth/elevation {az,el} (rad) (NULL: no output)
*                               (0.0<=azel[0]<2*pi,-pi/2<=azel[1]<=pi/2)
* return : elevation angle (rad)
*-----------------------------------------------------------------------------*/
double satazel(const double *pos, const double *e, double *azel)
{
	double az=0.0,el=PI/2.0,enu[3];

	if (pos[2]>-RE_WGS84) {
		ecef2enu(pos,e,enu);
		az=dot(enu,enu,2)<1E-12?0.0:atan2(enu[0],enu[1]);
		if (az<0.0) az+=2*PI;
		el=asin(enu[2]);
	}
	if (azel) {azel[0]=az; azel[1]=el;}
	return el;
}
double GetCombination(const GNSSDATA Obs, u2 Numb, double el, double Pse[3], double Carr[3])
{
	double L1,L2,L3,P1,P2,P3;

	L1 = Obs.ObsData[Numb].L[0];
	L2 = Obs.ObsData[Numb].L[1];
	L3 = Obs.ObsData[Numb].L[2];

	P1 = Obs.ObsData[Numb].P[0];
	P2 = Obs.ObsData[Numb].P[1];
	P3 = Obs.ObsData[Numb].P[2];
	if ((fabs(Pse[0]) > 1E-5 && (P1 < 1e-5)) || (fabs(Pse[1]) > 1E-5 && (P2 < 1e-5)) || (fabs(Pse[2]) > 1E-5 && (P2 < 1e-5))) 
	{
		return 0.0;
	}
	if ((fabs(Carr[0]) > 1E-5 && (fabs(L1) < 1e-5)) || (fabs(Carr[1]) > 1E-5 && (fabs(L2) < 1e-5)) || (fabs(Carr[2]) > 1E-5 && (fabs(L3) < 1e-5))) 
	{
		return 0.0;
	}
	u2 SYS,Sys_Num;
	char StrSys;
	Prn2Sys(Obs.ObsData[Numb].prn, SYS, Sys_Num, StrSys);
// 	if(SYS == SYS_BDS)
// 	{
// 		double Bias[3] = {0};
// 		GetBDSBias(el, Obs.ObsData[Numb].prn ,0, Bias[0]);
// 		GetBDSBias(el, Obs.ObsData[Numb].prn ,1, Bias[1]);
// 		GetBDSBias(el, Obs.ObsData[Numb].prn ,2, Bias[2]);
// 		P1 += Bias[0];
// 		P2 += Bias[1];
// 		P3 += Bias[2];
// 	}

	double PC,LC;
	PC = Pse[0] * P1 + Pse[1] * P2 + Pse[2] * P3;
	LC = Carr[0] * L1 + Carr[1] * L2 + Carr[2] * L3;
	return (PC-LC);
}
double  MWObs(const GNSSDATA Obs, u2 Numb, double el, int i, int j)
{
	double L1,L2,P1,P2;

	L1 = Obs.ObsData[Numb].L[i-1];
	L2 = Obs.ObsData[Numb].L[j-1];
	P1 = Obs.ObsData[Numb].P[i-1];
	P2 = Obs.ObsData[Numb].P[j-1];

	if(fabs(L1) < 1e-5 || fabs(L2) < 1e-5 || fabs(P1) < 1e-5 || fabs(P2) < 1e-5)
	{
		return 0.0;
	}
	else
	{
		u2 SYS,Sys_Num;
		char StrSys;
		double Pmw;
		u2 prn = Prn2Sys(Obs.ObsData[Numb].prn, SYS, Sys_Num, StrSys);
		if(Sys_Num == 1 && prn > 16) //1 2     3  2
		{
			double Lambda2 = CLIGHT / FREB5I;
			double freq = FREB5I;
			if(i == 1 && j == 2)
			{
				Pmw = (P1 * G_Fre[Sys_Num][i-1] + P2 * freq) *(G_Fre[Sys_Num][i-1] - freq) /((G_Fre[Sys_Num][i-1] + freq) * CLIGHT);
			}
			else if(i == 3 && j == 2)
			{
				Pmw = (P1 * G_Fre[Sys_Num][i-1] + P2 * freq) *(G_Fre[Sys_Num][i-1] - freq) /((G_Fre[Sys_Num][i-1] + freq) * CLIGHT);
			}
		}
		else
		{
			Pmw = (P1*G_Fre[Sys_Num][i-1] + P2*G_Fre[Sys_Num][j-1]) *(G_Fre[Sys_Num][i-1] - G_Fre[Sys_Num][j-1]) /((G_Fre[Sys_Num][i-1]+G_Fre[Sys_Num][j-1]) * CLIGHT);
		}
		return Pmw - (L1 - L2);
	}
}
//Coeff[0]:P1  Coeff[1]:P2 Coeff[2]:P3  Coeff[4]:PEWL  Coeff[5]:PWL 
bool GetMWObs(const GNSSDATA Obs, u2 Numb, double &MW, double Coeff[5], int Type, double Lamda[3])
{
	double L1,L2,L3,P1,P2,P3;

	L1 = Obs.ObsData[Numb].L[0];
	L2 = Obs.ObsData[Numb].L[1];
	L3 = Obs.ObsData[Numb].L[2];
	P1 = Obs.ObsData[Numb].P[0];
	P2 = Obs.ObsData[Numb].P[1];
	P3 = Obs.ObsData[Numb].P[2];

	//quick return
	if(P1 < 10000.0 || P2 < 10000.0 || fabs(L1) < PRECISION || fabs(L2) < PRECISION)
	{
		return false;
	}
	else if(Type == 0 && (P3 < 10000.0 || fabs(L3) < PRECISION))
	{
		return false;
	}
	else if((P3 < 10000.0 || fabs(L3) < PRECISION) && Type == 1 && fabs(Coeff[3]) > PRECISION)
	{
		return false;
	}

	double Ptemp;
	if(Type == 0) //EWL 
	{
		Ptemp = (Coeff[0] * P1 + Coeff[1] * P2 + Coeff[2] * P3);
		if(Obs.ObsData[Numb].prn > 32)  //BDS
		{
			MW = Ptemp / Lamda[0] - (L3 - L2);
		}
		else
		{
			MW = Ptemp / Lamda[0] - (L2 - L3);
		}
		return true;
	}
	else if(Type == 1) //WL
	{
		Ptemp = (Coeff[0] * P1 + Coeff[1] * P2 + Coeff[2] * P3);
		if(Obs.ObsData[Numb].prn > 32)  //BDS
		{
			Ptemp += Coeff[3] * (L3 - L2) * Lamda[0];	
		}
		else
		{
			Ptemp += Coeff[3] * (L2 - L3) * Lamda[0];
		}
		MW = Ptemp / Lamda[1] - (L1 - L3);
		return true;
	}
	else if(Type == 2) //NL
	{
		return false;  //not support now
	}
	else
	{
		return false;
	}
}
bool DetectRobust(double *Val, int number, double std, double &Ave, double &stdval)
{
	int *flag,Cycle=0;
	
	Ave = 0.0;
	stdval = 10.0;

	if(number < 1)
	{
		return false;
	}
	flag = imat(number, 1);
	memset(flag, 0x00, sizeof(int) * number);
	int Num = 0;
	while(Cycle < 3 && stdval > std)
	{
		double RMS = 0;
		Num = 0;
		Ave = 0;
		for(int i = 0; i < number; i++)
		{
			if(flag[i] == 1)
			{
				continue;
			}
			Ave += Val[i];
			RMS += Val[i] * Val[i];
			Num++;
		}
		if(Num > 0)
		{
			Ave = Ave / Num;
			RMS = RMS / Num;
			stdval = sqrt(RMS - Ave * Ave);
			for(int i = 0; i < number; i++)
			{
				if(flag[i] == 1)
				{
					continue;
				}
				if(fabs(Val[i] - Ave) > 3 * stdval || (fabs(Val[i] - Ave) > stdval && stdval > std))
				{
					flag[i] = 1;
				}
			}
		}
		else
		{
			Ave = 0.0;
			free(flag);
			return false;
		}
		Cycle++;
	}
	free(flag);
	if(Num >= 1 && stdval < std) //need modification
	{
// 		printf("Num %d  std %10.3f \n",Num, stdval);
		return true;
	}
	else
	{
		return false;
	}
}