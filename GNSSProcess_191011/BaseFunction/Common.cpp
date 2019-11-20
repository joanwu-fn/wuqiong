#include "../Common/BaseFunction.h"
const static double leaps[][7]={ /* leap seconds {y,m,d,h,m,s,utc-gpst,...} */
	{2017,1,1,0,0,0,-18},
	{2015,7,1,0,0,0,-17},
	{2012,7,1,0,0,0,-16},
	{2009,1,1,0,0,0,-15},
	{2006,1,1,0,0,0,-14},
	{1999,1,1,0,0,0,-13},
	{1997,7,1,0,0,0,-12},
	{1996,1,1,0,0,0,-11},
	{1994,7,1,0,0,0,-10},
	{1993,7,1,0,0,0, -9},
	{1992,7,1,0,0,0, -8},
	{1991,1,1,0,0,0, -7},
	{1990,1,1,0,0,0, -6},
	{1988,1,1,0,0,0, -5},
	{1985,7,1,0,0,0, -4},
	{1983,7,1,0,0,0, -3},
	{1982,7,1,0,0,0, -2},
	{1981,7,1,0,0,0, -1}
};
int cmpWLInfo(const void *p1, const void *p2)
{
	WL_Info *q1=(WL_Info *)p1,*q2=(WL_Info *)p2;
	return (int)q1->SSod - (int)q2->SSod;
}
int cmpobs(const void *p1, const void *p2)
{
	ObsData_t *q1=(ObsData_t *)p1,*q2=(ObsData_t *)p2;
	return (int)q1->prn-(int)q2->prn;
}
double str2num(const char *s, int i, int n)
{
	double value;
	char str[256],*p=str;

	if (i<0||(int)strlen(s)<i||(int)sizeof(str)-1<n) return 0.0;
	for (s+=i;*s&&--n>=0;s++) *p++=*s=='d'||*s=='D'?'E':*s; *p='\0';
	return sscanf_s(str,"%lf",&value)==1?value:0.0;
}

/* inner product ---------------------------------------------------------------
* inner product of vectors
* args   : double *a,*b     I   vector a,b (n x 1)
*          int    n         I   size of vector a,b
* return : a'*b
*-----------------------------------------------------------------------------*/
double dot(const double *a, const double *b, int n)
{
	double c=0.0;

	while (--n>=0) c+=a[n]*b[n];
	return c;
}
/* euclid norm -----------------------------------------------------------------
* euclid norm of vector
* args   : double *a        I   vector a (n x 1)
*          int    n         I   size of vector a
* return : || a ||
*-----------------------------------------------------------------------------*/
double norm(const double *a, int n)
{
	return sqrt(dot(a,a,n));
}
/* outer product of 3d vectors -------------------------------------------------
* outer product of 3d vectors 
* args   : double *a,*b     I   vector a,b (3 x 1)
*          double *c        O   outer product (a x b) (3 x 1)
* return : none
*-----------------------------------------------------------------------------*/
void cross3(const double *a, const double *b, double *c)
{
	c[0]=a[1]*b[2]-a[2]*b[1];
	c[1]=a[2]*b[0]-a[0]*b[2];
	c[2]=a[0]*b[1]-a[1]*b[0];
}
/* normalize 3d vector ---------------------------------------------------------
* normalize 3d vector
* args   : double *a        I   vector a (3 x 1)
*          double *b        O   normlized vector (3 x 1) || b || = 1
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
int normv3(const double *a, double *b)
{
	double r;
	if ((r=norm(a,3))<=0.0) return 0;
	b[0]=a[0]/r;
	b[1]=a[1]/r;
	b[2]=a[2]/r;
	return 1;
}
//取字符串的子串
void xstrmid(char *szDest, char *szSrc, int nPos, int nCount)
{
	int		i;
	char	*str;
	char	c;
	int		iCount;
	int	iTemp;

	str	= szSrc + nPos;

	iTemp	= strlen(szSrc) - nPos + 1;
	iCount	= (nCount < iTemp ? nCount : iTemp);

	for (i=0; i < iCount; i++) {
		c	= str[i];
		if (c) {
			szDest[i]	= c;
		}
		else {
			szDest[i]	= '\0';
			break;
		}
	}

	szDest[(0 > iCount ? 0 : iCount)]	= '\0';
}

//去除字符串的前导空格和尾部空格及换行
void xstrtrim (char* pszDest, char* pszSrc)
{
	int	len;
	int	i;

	len	= strlen (pszSrc);

	//去除前导的空格
	for (i=0; i<len; i++)
	{
		if (pszSrc[i] == ' ')
		{
		}
		else
		{
			strcpy(pszDest, pszSrc+i);
			len	-= i;
			break;
		}
	}
	if (len == i)
	{
		return;
	}

	//去除尾部空格及换行
	for (i=len-1; i>0; i--)
	{
		if (pszDest[i] == ' ' || pszDest[i] == '\n' || pszDest[i] == '\r')
		{
		}
		else
		{
			pszDest[i+1]	= 0;
			break;
		}
	}
}

//扩展字符串到给定长度
void xstrStretch (char* pszSrc, unsigned int lenEnd)
{
	size_t	len;
	size_t	i;

	len	= strlen(pszSrc);

	for (i = len; i <= lenEnd; i++)
	{
		pszSrc[i] = ' ';
	}
	pszSrc[lenEnd + 1] = '\0';
}

int xstrReplace (char szSrc[], char ch2Re, char chRe, unsigned int maxCh)
{
	size_t len = strlen(szSrc);
	size_t i = 0;

	unsigned int nCh = 0;
	for (i = 0; nCh < maxCh && i < len; i++)
	{
		if (ch2Re == szSrc[i])
		{
			szSrc[i] = chRe;
			nCh++;
		}
	}

	return nCh;
}
int GetHighRateTime(FILE *ft, gtime_t &gtime, int type, char line1[MAXLENGTH])
{
	char    	aux[100];
	struct tm	tm;
	double seconds;
	char 		line[MAXLENGTH];
	int 		len = 0;
	TTime       Temp;
	if(type == 0)
	{
		sprintf(line, "%s",line1);
		getstr(aux,line,3,4);
		tm.tm_year = atoi(aux)-1900;
		getstr(aux,line,8,2);
		tm.tm_mon  = atoi(aux)-1;
		getstr(aux,line,11,2);
		tm.tm_mday = atoi(aux);
		getstr(aux,line,14,2);
		tm.tm_hour = atoi(aux);
		getstr(aux,line,17,2);
		tm.tm_min  = atoi(aux);
		getstr(aux,line,20,11);
		tm.tm_sec  = atoi(aux);
		seconds = atof(aux);
		Temp.MJDN = MJDN(&tm);
		Temp.SoD = tm.tm_hour*3600 + tm.tm_min*60 + seconds;
		gtime = gltime2rtktime(Temp);
		return 0; //get
	}
	else
	{
		while(fgets(line, MAXLENGTH, ft) != NULL)
		{
			if (line[0]=='*') 
			{ // Epoch starting
				// *  2008 08 27 00 00 00.00000000
				getstr(aux,line,3,4);
				tm.tm_year = atoi(aux)-1900;
				getstr(aux,line,8,2);
				tm.tm_mon  = atoi(aux)-1;
				getstr(aux,line,11,2);
				tm.tm_mday = atoi(aux);
				getstr(aux,line,14,2);
				tm.tm_hour = atoi(aux);
				getstr(aux,line,17,2);
				tm.tm_min  = atoi(aux);
				getstr(aux,line,20,11);
				tm.tm_sec  = atoi(aux);
				seconds = atof(aux);
				Temp.MJDN = MJDN(&tm);
				Temp.SoD = tm.tm_hour*3600 + tm.tm_min*60 + seconds;
				gtime = gltime2rtktime(Temp);
				return 0; //get
			} 
		}	
	}
	return -1;
}
void GetHighRateProducts(FILE *ft, gtime_t GNowTime, IGS_Orbit_t *Orb, IGS_Clk_t *Clk)
{
	char buff[MAXLENGTH],StrSys;
	double ep[6] = {0},pos[3];
	gtime_t NowTime= {0};
	int sys;
	u2 nPrn=0, temp,prn,NAVSYS = SYS_GPS;
	while(fgets(buff, MAXLENGTH, ft) != NULL)
	{
		if(buff[0] != '*' && buff[0] != 'P')
		{
			continue;
		}
		if(buff[0] == '*')
		{	
			GetHighRateTime(ft, GNowTime, 0, buff);
			return;
		}
		else if(buff[0] == 'P')
		{
			StrSys =   buff[1];
			if((sys = JudgeSys(StrSys, SYS_GPS| SYS_BDS, &NAVSYS)) == -1)
			{
				continue;
			}
			temp = (int)str2num(buff, 2, 2);
			prn  = SatNo(temp, NAVSYS);
			Orb->prn[prn-1]       = prn;
			Orb->Orbit[prn-1][0]  = str2num(buff,  4, 14) * 1000.0;
			Orb->Orbit[prn-1][1]  = str2num(buff, 18, 14) * 1000.0;
			Orb->Orbit[prn-1][2]  = str2num(buff, 32, 14) * 1000.0;

			Clk->prn[prn-1]       = prn;
			Clk->SatClk[prn-1]    = str2num(buff, 46, 14) * 1000.0;
		}
	}
	return;
}
//------------------- --------------------------------------
double enorm8(double Vector[])
{
	return sqrt(Vector[0] * Vector[0] + Vector[1] * Vector[1] + Vector[2] * Vector[2]);
}

double enorm8(Coordinate coor)
{
	return sqrt(coor.X * coor.X + coor.Y * coor.Y + coor.Z * coor.Z);
}
Coordinate operator %(Coordinate coor1,Coordinate coor2)
{
	return Coordinate(coor1.Y * coor2.Z - coor1.Z * coor2.Y, coor1.Z * coor2.X - coor1.X * coor2.Z, coor1.X * coor2.Y - coor1.Y * coor2.X);	
}

Coordinate unit(Coordinate coor)
{
	double doubleTemp = 1.0 / enorm8(coor);

	return Coordinate(coor.X * doubleTemp, coor.Y * doubleTemp, coor.Z * doubleTemp);
}
/* purpose  : dot product of two vectors       */
/* parameter: n -- dimension of vectors        */
/*            v1, v2 -- input vectors          */
/*            dot -- dot product of v1 and v2  */
double Dot(int n, const double v1[], const double v2[])
{
	double dot = 0.0;
	int i = 0;
	for (i = 0; i < n; i++)
	{
		dot += v1[i] * v2[i];
	}

	return dot;
}
/*-------------------------------------------------------------*/
/* convert tsec in GPS to tsec in UTC                          */
/* GPS is ahead of UTC                                         */
/* UTC is behind GPS                                           */
/* GPSleap() is (so far) positice (and increasing)             */
/* so, must subtract GPSleap from GPS to get UTC               */
/*-------------------------------------------------------------*/
double GPS2UTC(TTime& mjd)
{
	int i;
	double ileaps = -18;
	for (i=0;i<(int)sizeof(leaps)/(int)sizeof(*leaps);i++) 
	{
		TTime tt = rtktime2gltime(epoch2time(leaps[i]));
		double diff = (mjd.MJDN - tt.MJDN) * 86400 + mjd.SoD - tt.SoD;
		if (diff > 0)
		{
			ileaps = leaps[i][6];
			break;
		}
	}
	return mjd.SoD - ileaps;
}
/*-------------------------------------------------------------*/
/* convert tsec in GPS to tsec in TT                           */
/*-------------------------------------------------------------*/
double GPS2TT(double tsec)
{
	return tsec + 51.184;			   /* fixed offset */
}
/*----------------------------------------------------------------*/
/* rotate coordinate axes about 3 axis by angle of theta radians  */
/*                              x, y, z transformed into u, v, w  */
/*----------------------------------------------------------------*/
void rot3(double theta, double x, double y, double z, double &u, double &v, double &w)
{
	double s = sin(theta);
	double c = cos(theta);

	u = c * x + s * y;
	v = c * y - s * x;
	w = z;
}
/*----------------------------------------------------------------*/
/* convert mjd in GPS time to Greenwich hour angle (in radians)   */
/*----------------------------------------------------------------*/
void getghar(TTime& mjd, double &ghar)
{
	/* need UT to get sidereal time */
	double tsecgps = mjd.SoD;										/* GPS time (sec of day)           */
	double tsecutc = GPS2UTC(mjd);													/* UTC time (sec of day)           */
	double fmjdutc = tsecutc / 86400.0;												/* UTC time (fract. day)           */

	double d = (mjd.MJDN - 51544) + (fmjdutc - 0.50);								/* days since J2000                */

	/* Greenwich hour angle for J2000  (12:00:00 on 1 Jan 2000) */
	double ghad = 280.46061837504 + 360.9856473662862 * d;							/* corrn.   (+digits)			   */

	/* normalize to 0-360 and convert to radians */
	int  i    = (int)(ghad / 360.0);
	ghar      = (ghad - i * 360.0) * D2R;
	while (ghar >= 2 * PI)
	{
		ghar = ghar - 2 * PI;
	}
	while (ghar < 0.0)
	{
		ghar = ghar + 2 * PI;
	}
}
/*----------------------------------------------------------------*/
/* get low-precision, geocentric coordinates for sun in ECEF      */
/* input, mjd is Modified Julian Date in GPS time                 */
/* output, rs[], is geocentric solar position vector [m] in ECEF  */
/*----------------------------------------------------------------*/
void sunxyz(TTime& mjd, double rs[])
{
	/* mean elements for year 2000, sun ecliptic orbit wrt. Earth */
	double obe     = 23.43929111 * D2R;												/* obliquity of the J2000 ecliptic */
	double sobe    = sin(obe);
	double cobe    = cos(obe);
	double opod    = 282.94;														/* RAAN + arg.peri. (deg.)         */

	/* use TT for solar ephemerides */
	double tsecgps = mjd.SoD;										/* GPS time (sec of day)           */
	double tsectt  = GPS2TT(tsecgps);												/* TT  time (sec of day)           */
	double fmjdtt  = tsectt / 86400.0;												/* TT  time (fract. day)           */

	/* julian centuries since 1.5 january 2000 (J2000) */
	/*    (note also low precision use of mjd --> tjd) */
	double tjdtt   = mjd.MJDN + fmjdtt + 2400000.5;									/* Julian Date, TT                 */
	double t       = (tjdtt - 2451545.0) / 36525.0;									/* Julian centuries, TT            */
	double emdeg   = 357.5256 + 35999.049 * t;										/* degrees                         */
	double em      = emdeg * D2R;													/* radians                         */
	double em2     = em + em;														/* radians                         */

	/* series expansions in mean anomaly, em */
	double r       = (149.619 - 2.499 * cos(em) - 0.021 * cos(em2)) * 1.0e9;		/* m.                              */
	double slond   = opod + emdeg + (6892.0 * sin(em) + 72.0 * sin(em2)) / 3600.0;

	/* precession of equinox */
	slond         += 1.3972 * t;													/* degrees                         */

	/* position vector of sun (mean equinox & ecliptic of J2000) (EME2000, ICRF) */
	/*                                    (plus long. advance due to precession) */
	double slon    = slond * D2R;													/* radians                         */
	double sslon   = sin(slon);
	double cslon   = cos(slon);

	double rs1     = r * cslon;														/* meters                          */
	double rs2     = r * sslon * cobe;												/* meters                          */
	double rs3     = r * sslon * sobe;												/* meters                          */

	/* convert position vector of sun to ECEF (ignore polar motion / LOD) */
	double ghar = 0.0;
	getghar(mjd, ghar);
	rot3(ghar, rs1, rs2, rs3, rs[0], rs[1], rs[2]);
}