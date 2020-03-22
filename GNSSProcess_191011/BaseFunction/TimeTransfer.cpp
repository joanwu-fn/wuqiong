#include "../Common/BaseFunction.h"

const static double gpst0[]={1980,1, 6,0,0,0}; /* gps time reference */
const static double gst0 []={1999,8,22,0,0,0}; /* galileo system time reference */
const static double bdst0[]={2006,1, 1,0,0,0}; /* bds system time reference */
const static double leaps[][7]={ /* leap seconds {y,m,d,h,m,s,utc-gpst,...} */
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

/*****************************************************************************
* Name        : modulo
* Description : Does the modulo mod of a number
* Parameters  :
* Name                           |Da|Unit|Description
* double  a                       I  N/A  Number
* double  mod		               I  N/A  Modulo base
* Returned value (double)         O  N/A  Modulo
*****************************************************************************/
double modulo (double a, double mod) {
	return a - ((int)(a/mod))*mod;
}
/* add time --------------------------------------------------------------------
* add time to gtime_t struct
* args   : gtime_t t        I   gtime_t struct
*          double sec       I   time to add (s)
* return : gtime_t struct (t+sec)
*-----------------------------------------------------------------------------*/
gtime_t timeadd(gtime_t t, double sec)
{
	double tt;

	t.sec+=sec; tt=floor(t.sec); t.time+=(int)tt; t.sec-=tt;
	return t;
}
/* time difference -------------------------------------------------------------
* difference between gtime_t structs
* args   : gtime_t t1,t2    I   gtime_t structs
* return : time difference (t1-t2) (s)
*-----------------------------------------------------------------------------*/
double timediff(gtime_t t1, gtime_t t2)
{
	return difftime(t1.time,t2.time)+t1.sec-t2.sec;
}

gtime_t epoch2time(const double *ep)
{
	const int doy[]={1,32,60,91,121,152,182,213,244,274,305,335};
	gtime_t time={0};
	int days,sec,year=(int)ep[0],mon=(int)ep[1],day=(int)ep[2];

	if (year<1970||2099<year||mon<1||12<mon) return time;

	/* leap year if year%4==0 in 1901-2099 */
	days=(year-1970)*365 + (year-1969)/4 + doy[mon-1]+day-2+( year%4 == 0 && mon>=3?1:0);
	sec=(int)floor(ep[5]);
	time.time=(time_t)days*86400+(int)ep[3]*3600+(int)ep[4]*60+sec;
	time.sec=ep[5]-sec;
	return time;
}
/* time to calendar day/time ---------------------------------------------------
* convert gtime_t struct to calendar day/time
* args   : gtime_t t        I   gtime_t struct
*          double *ep       O   day/time {year,month,day,hour,min,sec}
* return : none
* notes  : proper in 1970-2037 or 1970-2099 (64bit time_t)
*-----------------------------------------------------------------------------*/
void time2epoch(gtime_t t, double *ep)
{
	const int mday[]={ /* # of days in a month */
		31,28,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,
		31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31
	};
	int days,sec,mon,day;

	/* leap year if year%4==0 in 1901-2099 */
	days=(int)(t.time/86400);
	sec=(int)(t.time-(time_t)days*86400);
	for (day=days%1461,mon=0;mon<48;mon++) {
		if (day>=mday[mon]) day-=mday[mon]; else break;
	}
	ep[0]=1970+days/1461*4+mon/12; ep[1]=mon%12+1; ep[2]=day+1;
	ep[3]=sec/3600; ep[4]=sec%3600/60; ep[5]=sec%60+t.sec;
}
gtime_t gpst2time(int week, double sec)
{
	gtime_t t=epoch2time(gpst0);

	if (sec<-1E9||1E9<sec) sec=0.0;
	t.time+=86400*7*week+(int)sec;
	t.sec=sec-(int)sec;
	return t;
}
/* time to gps time ------------------------------------------------------------
* convert gtime_t struct to week and tow in gps time
* args   : gtime_t t        I   gtime_t struct
*          int    *week     IO  week number in gps time (NULL: no output)
* return : time of week in gps time (s)
*-----------------------------------------------------------------------------*/
double time2gpst(gtime_t t, int *week)
{
	gtime_t t0=epoch2time(gpst0);
	time_t sec=t.time-t0.time;
	int w=(int)(sec/(86400*7));

	if (week) *week=w;
	return (double)(sec-w*86400*7)+t.sec;
}
/* time to day of year ---------------------------------------------------------
* convert time to day of year
* args   : gtime_t t        I   gtime_t struct
* return : day of year (days)
*-----------------------------------------------------------------------------*/
double time2doy(gtime_t t)
{
	double ep[6];

	time2epoch(t,ep);
	ep[1]=ep[2]=1.0; ep[3]=ep[4]=ep[5]=0.0;
	return timediff(t,epoch2time(ep))/86400.0+1.0;
}
/*****************************************************************************
* Name        : MJDN
* Description : Get the modified julian day from YY/MM/DD
* Parameters  :
* Name                           |Da|Unit|Description
* struct tm  tm                   I  N/A  Date
* Returned value (int)            O  N/A  Integer days in Modified Julian Date
*****************************************************************************/
int MJDN(struct tm *tm) {
	int y,m;

	y = tm->tm_year + 1900;
	m = tm->tm_mon + 1;
	if (m<=2) {
		y--;
		m+=12;
	}

	return (int)(365.25*y) + (int)(30.6001*(m+1)) + tm->tm_mday - 679019;
}
double Ttime2GPS(TTime T, int *week)
{
	*week = (T.MJDN - 44244) / 7;
	return (T.SoD + ((T.MJDN - 44244)%7)*86400);
}
/*****************************************************************************
* Name        : cal2t
* Description : Transforms between Calendar time (year, month, day of month,
*               hour, minute and second) to TTime structure
* Parameters  :
* Name                           |Da|Unit|Description
* int  year                       I  y    Year (can be of 4 or 2 digits)
* int  month                      I  m    Month
* int  day                        I  d    Day of month
* int  hour                       I  h    Hours of day
* int  minute                     I  m    Minutes in hour
* double  second                  I  s    Seconds of minute
* Returned value (TTime)          O  N/A  TTime structure
*****************************************************************************/
TTime cal2t (int year, int month, int day, int hour, int minute, double second) {
	struct tm	tm;
	TTime		t;

	if (year<100) { // 2 digits
		if (year<=70) { // 20XX
			tm.tm_year = year + 100;
		} else { // 19XX
			tm.tm_year = year;
		}
	} else { // 4 digits
		tm.tm_year = year - 1900;
	}

	tm.tm_mon = month - 1;
	tm.tm_mday = day;
	tm.tm_hour = hour;
	tm.tm_min = minute;
	tm.tm_sec = (int)second;
	t.MJDN = MJDN (&tm);
	t.SoD = tm.tm_hour*3600 + tm.tm_min*60 + second;

	return t;
}


/*****************************************************************************
* Name        : t2doy
* Description : Get the Year/doy from a TTime structure
* Parameters  :
* Name                           |Da|Unit|Description
* TTime  *t                       I  N/A  TTime structure
* int  *year                      O  N/A  Year
* int  *doy                       O  N/A  Day of year
*****************************************************************************/
void t2doy(TTime *t, int *year, double *doy) {
	const int MJD_1980 = 44239;
	int     day;
	int     ref = 1980;
	int     day_aux;
	int     iy;

	day = t->MJDN - MJD_1980;

	day_aux = modulo(day,1461);
	if (day_aux>365) {
		day_aux -= 366;
		*doy = modulo (day_aux,365) + 1;
		iy = day_aux/365 + 1;
	} else {
		*doy = day_aux + 1;
		iy = 0;
	}
	*year = ref + 4*(day/1461) + iy;
	*doy += t->SoD/86400;
}

/*****************************************************************************
* Name        : t2tm
* Description : Get the calendar time from a TTime structure
* Parameters  :
* Name                           |Da|Unit|Description
* TTime  *t                       I  N/A  TTime structure
* tm  *tm                         O  N/A  tm structure
* double  *seconds                O  N/A  Seconds of minute in double format
*****************************************************************************/
void t2tm(TTime *t, struct tm *tm, double *seconds) {
	int		year;
	double	doy;
	time_t	tt;

	t2doy(t,&year,&doy);
	tm->tm_year = year - 1900;
	tm->tm_mday = (int)(doy);
	tm->tm_mon = 0;
	tm->tm_hour = (int)(t->SoD/3600);
	tm->tm_min = (int)(modulo(t->SoD,3600)/60);
	tm->tm_sec = (int)(modulo(t->SoD,60));
	tm->tm_isdst = -1;
	*seconds =  modulo(t->SoD,60);
	tt = mktime(tm);
	//localtime_r(&tt,tm);
	memcpy(tm,localtime(&tt),sizeof(tm));
}
/* convert glab time to rtklib time gtime_t-------------------------------------
* args   : glabtime t    I   TTime struct
* return : rtklibtime gt O   gtime_t sruct
*-----------------------------------------------------------------------------*/
gtime_t gltime2rtktime(TTime t){
	struct tm	tm;
	double second;
	double ep[6];
	t2tm(&t,&tm,&second);
	ep[0] = tm.tm_year+1900;
	ep[1] = tm.tm_mon+1;
	ep[2] = tm.tm_mday;
	ep[3] = tm.tm_hour;
	ep[4] = tm.tm_min;
	ep[5] = second;
	return epoch2time(ep);
}
/*****************************************************************************
* Name        : rtktime2gltime
* Description : translate gtime_t time t to TTime time
* Parameters  :
* Name                           |Da|Unit|Description
* gtime_t t                       I  N/A   time t
*****************************************************************************/
TTime rtktime2gltime(gtime_t t){
	TTime tt;
	double ep[6];
	time2epoch(t, ep);
	tt=cal2t(ep[0],ep[1],ep[2],ep[3],ep[4],ep[5]);
	return tt;
}
/* gpstime to utc --------------------------------------------------------------
* convert gpstime to utc considering leap seconds
* args   : gtime_t t        I   time expressed in gpstime
* return : time expressed in utc
* notes  : ignore slight time offset under 100 ns
*-----------------------------------------------------------------------------*/
gtime_t gpst2utc(gtime_t t)
{
	gtime_t tu;
	int i;

	for (i=0;i<(int)sizeof(leaps)/(int)sizeof(*leaps);i++) {
		tu=timeadd(t,leaps[i][6]);
		if (timediff(tu,epoch2time(leaps[i]))>=0.0) return tu;
	}
	return t;
}
/* utc to gpstime --------------------------------------------------------------
* convert utc to gpstime considering leap seconds
* args   : gtime_t t        I   time expressed in utc
* return : time expressed in gpstime
* notes  : ignore slight time offset under 100 ns
*-----------------------------------------------------------------------------*/
gtime_t utc2gpst(gtime_t t)
{
	int i;

	for (i=0;i<(int)sizeof(leaps)/(int)sizeof(*leaps);i++) {
		if (timediff(t,epoch2time(leaps[i]))>=0.0) return timeadd(t,-leaps[i][6]);
	}
	return t;
}
/* time to day and sec -------------------------------------------------------*/
double time2sec(gtime_t time, gtime_t *day)
{
	double ep[6],sec;
	time2epoch(time,ep);
	sec=ep[3]*3600.0+ep[4]*60.0+ep[5];
	ep[3]=ep[4]=ep[5]=0.0;
	*day=epoch2time(ep);
	return sec;
}
/* utc to gmst -----------------------------------------------------------------
* convert utc to gmst (Greenwich mean sidereal time)
* args   : gtime_t t        I   time expressed in utc
*          double ut1_utc   I   UT1-UTC (s)
* return : gmst (rad)
*-----------------------------------------------------------------------------*/
double utc2gmst(gtime_t t, double ut1_utc)
{
	const double ep2000[]={2000,1,1,12,0,0};
	gtime_t tut,tut0;
	double ut,t1,t2,t3,gmst0,gmst;

	tut=timeadd(t,ut1_utc);
	ut=time2sec(tut,&tut0);
	t1=timediff(tut0,epoch2time(ep2000))/86400.0/36525.0;
	t2=t1*t1; t3=t2*t1;
	gmst0=24110.54841+8640184.812866*t1+0.093104*t2-6.2E-6*t3;
	gmst=gmst0+1.002737909350795*ut;

	return fmod(gmst,86400.0)*PI/43200.0; /* 0 <= gmst <= 2*PI */
}
/* string to time --------------------------------------------------------------
* convert substring in string to gtime_t struct
* args   : char   *s        I   string ("... yyyy mm dd hh mm ss ...")
*          int    i,n       I   substring position and width
*          gtime_t *t       O   gtime_t struct
* return : status (0:ok,0>:error)
*-----------------------------------------------------------------------------*/
int str2time(const char *s, int i, int n, gtime_t *t)
{
	double ep[6];
	char str[256],*p=str;

	if (i<0||(int)strlen(s)<i||(int)sizeof(str)-1<i) return -1;
	for (s+=i;*s&&--n>=0;) *p++=*s++; *p='\0';
	if (sscanf(str,"%lf %lf %lf %lf %lf %lf",ep,ep+1,ep+2,ep+3,ep+4,ep+5)<6)
		return -1;
	if (ep[0]<100.0) ep[0]+=ep[0]<80.0?2000.0:1900.0;
	*t=epoch2time(ep);
	return 0;
}

void TTimeToCommonTime(const TTime* pmjd, COMMONTIME* pct)
{
	int a, b, c, d, e;

	a = pmjd->MJDN  + 2400001;
	b = a + 1537;
	c = static_cast<int>((b - 122.1) / 365.25);
	d = static_cast<int>(365.25 * c);
	e = static_cast<int>((b - d) / 30.6001);

	pct->day = b - d - static_cast<int>(30.6001 * e);
	pct->month = e - 1 - 12 * (e / 14);
	pct->year = c - 4715 - ((7 + pct->month) / 10);
	pct->hour = (int)(pmjd->SoD / 3600);
	pct->minute = ((int)pmjd->SoD % 3600) / 60;
	pct->second = ((int)pmjd->SoD % 3600) % 60;
}