#ifndef _COMMON_HEADER_20150118_
#define _COMMON_HEADER_20150118_
#include <time.h>
#include <string.h>
#include <stdio.h>
typedef unsigned char            u1;
typedef unsigned short           u2;
typedef unsigned long            u4;

//宏定义
#define SQR(x)      ((x)*(x))
#define      INVALID_VALUE       -999.999
#define      ENABLEGNSS           4
#define      NEFREQ               3
#define      MAXSTATION           3000

#define      MAXLENGTH             500               //字符数组最大值
#define      MAXCHN               (20*ENABLEGNSS)     // MAX visible number os satellite
#define      MAXENTRY_NUMBER       30                 //max types of obs


#define      MAXSYS               4
#define      SYS_NONE             0x00
#define      SYS_GPS              0x01
#define      SYS_BDS              0x02
#define      SYS_GAL              0x04
#define      SYS_GLO              0x08
#define      SYS_ALL              (SYS_GPS|SYS_BDS|SYS_GAL|SYS_GLO)

#define      MAXSAT            (NSATGPS+NSATBDS+NSATGAL+NSATGLO)

#define      MINPRNGPS             1                   /* min satellite PRN number of GPS */
#define      MAXPRNGPS             32                  /* max satellite PRN number of GPS */
#define      NSATGPS               32                  /* number of GPS satellites */

#define      MINPRNBDS             1                   /* min satellite sat number of Compass */
#define      MAXPRNBDS             16                  /* max satellite sat number of Compass */
#define      NSATBDS               MAXPRNBDS           /* number of Compass satellites */

#define      MINPRNGAL             1
#define      MAXPRNGAL             30
#define      NSATGAL               30

#define      MINPRNGLO             0
#define      MAXPRNGLO             0
#define      NSATGLO               0

#define FREQ1       1.57542E9           /* L1/E1  frequency (Hz) */
#define FREQ2       1.22760E9           /* L2     frequency (Hz) */
#define FREQ5       1.17645E9           /* L5/E5a frequency (Hz) */

#define FREB1       1.561098E9           /* L1/E1  frequency (Hz) */
#define FREB2       1.207140E9           /* L2     frequency (Hz) */
#define FREB5       1.268520E9           /* L3     frequency (Hz) */
#define FREB5I      1.17645E9

#define FREE1       1.57542E9
#define FREE5       1.17645E9
#define FREE7       1.207140E9

#define RE_WGS84    6378137.0           /* earth semimajor axis (WGS84) (m) */
#define FE_WGS84    (1.0/298.257223563) /* earth flattening (WGS84) */
#define HION        350000.0            /* ionosphere height (m) */
#define AU          149597870691.0      /* 1 AU (m) */
#define AS2R        (D2R/3600.0)        /* arc sec to radian */
#define M2NS						(1e9 / 299792458.0)		    
#define NS2M						(299792458.0 / 1e9)							
#define MAXANT      64                  /* max length of station name/antenna type */
#define PRECISION   (1e-7)

//常量
const double  PI      =   3.1415926535897932;    //pi
const double  CLIGHT  =   299792458.0;           //光速
const double  D2R     =   (PI/180.0);          /* deg to rad */
const double  R2D     =   (180.0/PI);          /* rad to deg */
const double  G_MU[4]   =  {3.9860050E14, 3.986004418e14, 3.986004418E14, 3.9860044E14}; //GPS BDS GLALIEO GLONASS
const double  G_OMGE[4] =  {7.2921151467E-5 , 7.2921150000e-5, 7.2921151467E-5 , 7.292115E-5};   //GPS BDS GLALIEO GLONASS
const u2      MAXITER   =  10;
const double  G_Fre[4][3]={{FREQ1, FREQ2, FREQ5},{FREB1, FREB2, FREB5},{0,0,0},{0,0,0}};
const double  G_Lam[4][3]={{CLIGHT/FREQ1, CLIGHT/FREQ2, CLIGHT/FREQ5},{CLIGHT/FREB1, CLIGHT/FREB2, CLIGHT/FREB5},{CLIGHT/FREE1,CLIGHT/FREE5,CLIGHT/FREE7},{0,0,0}};
const double  G_WL[4] ={ CLIGHT/(FREQ1-FREQ2), CLIGHT/(FREB1-FREB2), 0, 0};
const double  G_EWL[4]={ CLIGHT/(FREQ2-FREQ5), CLIGHT/(FREB5-FREB2), 0, 0};
const double  G_B1B3Default[14] = {7.8086,  -5.5289,  -2.8895, -0.7787, -6.1334,  1.9368, 8.0206, 5.4828, 1.1097, 1.0888, -2.3894, -1.8875,-1.4572, 1.0909};
const double  G_B1B2Default[14] = {16.0264, 5.7435, 5.0819, 3.9562,  -0.1673, 1.1228, 4.1743, 2.6981, -5.7684, -4.4165, -7.1783, -3.9031,-5.7404, -4.0928};

enum PROCESS_MODE
{
	PMODE_SINGLE       = 0,
	PMODE_PPP_KINEMA   = 1,
	PMODE_PPP_STATIC   = 2,
	PMODE_PPP_FIXHOL   = 3,
	PMODE_BL_KINEMA    = 4,
	PMODE_BL_STATIC    = 5,
	PMODE_BL_FIXHOL    = 6
};

enum SPPTYPE
{
	SPPMODE_SINGLE = 0,
	SPPMODE_IFLC   = 1
};
enum MPTYPE
{
	MP1   =  0,
	MP2   =  1,
	MP3   =  2
};

//结构体

typedef struct {        // time struct
	time_t time;        // time (s) expressed by standard time_t
	double sec;         // fraction of second under 1 s 
} gtime_t;

// Time handling structure
typedef struct {
	int				MJDN;
	double			SoD;
} TTime;

typedef struct ST_COORDINATE 
{
	ST_COORDINATE()
	{
		X = 0.0;
		Y = 0.0;
		Z = 0.0;
	}
	ST_COORDINATE(double a, double b, double c)
	{
		X = a;
		Y = b;
		Z = c;
	}

	double X;
	double Y;
	double Z;
} Coordinate, CoorGeo;

/////////////////////////////////////////////////////////////////////////////////////////////////
//structures for time
typedef struct ST_PRECISETIME
{
	ST_PRECISETIME()
	{
		sn = 0;
		tos = 0.0;
	}

	long   sn;
	double tos;
}PRECISETIME;


typedef struct ST_MODIFYJULIANDAY
{
	ST_MODIFYJULIANDAY()
	{
		day = 0;
	}

	long        day;
	PRECISETIME tod;
}MODIFYJULIANDAY;

//for gregorian calendar time structure
typedef struct ST_COMMONTIME
{
	ST_COMMONTIME()
	{
		year   = 0;
		month  = 0;
		day    = 0;
		hour   = 0;
		minute = 0;
		second = 0.0;
	}
	ST_COMMONTIME(int ye, int mo, int d, int ho, int min, double sec)
	{
		year   = ye;
		month  = mo;
		day    = d;
		hour   = ho;
		minute = min;
		second = sec;
	}

	int   year;
	int  month;
	int  day;
	int  hour;
	int  minute;
	double second;
}COMMONTIME;

struct ObsData_t{
	u2     prn;
	double L[NEFREQ];
	double P[NEFREQ];
	double C[NEFREQ];
};
typedef struct GNSSDATA
{
	gtime_t     ObsTime;    //obs time
	u2          iSite;      // site number
	u2          NumSats;    //satellite number
	ObsData_t   ObsData[MAXCHN];
};

typedef struct 
{
	int    flags;            /* GPSEPHF_xxx */
	int    prn;	             /*  SV ID   ICD-GPS data position */
	int    IODE;             /*          [s2w3b01-08]              */
	int    URAindex;         /*  [1..15] [s1w3b13-16]              */
	int    SVhealth;         /*          [s1w3b17-22]              */
	int    IODC;             /*          [s1w3b23-32,w8b01-08]     */
	double toes;
	gtime_t TOE;
	gtime_t TOC;
	double clock_bias;       /*  [s]     [s1w10b1-22, af0]         */
	double clock_drift;      /*  [s/s]   [s1w9b09-24, af1]         */
	double clock_driftrate;  /*  [s/s^2] [s1w9b01-08, af2]         */
	double Crs;              /*  [m]     [s2w3b09-24]              */
	double Delta_n;          /*  [rad/s] [s2w4b01-16 * Pi]         */
	double M0;               /*  [rad]   [s2w4b17-24,w5b01-24 * Pi]*/
	double Cuc;              /*  [rad]   [s2w6b01-16]              */
	double e;                /*          [s2w6b17-24,w6b01-24]     */
	double Cus;              /*  [rad]   [s2w8b01-16]              */
	double sqrt_A;           /*  [m^0.5] [s2w8b16-24,w9b01-24]     */
	double Cic;              /*  [rad]   [s3w3b01-16]              */
	double OMEGA0;           /*  [rad]   [s3w3b17-24,w4b01-24 * Pi]*/
	double Cis;              /*  [rad]   [s3w5b01-16]              */
	double i0;               /*  [rad]   [s3w5b17-24,w6b01-24 * Pi]*/
	double Crc;              /*  [m]     [s3w701-16]               */
	double omega;            /*  [rad]   [s3w7b17-24,w8b01-24 * Pi]*/
	double OMEGADOT;         /*  [rad/s] [s3w9b01-24 * Pi]         */
	double IDOT;             /*  [rad/s] [s3w10b9-22 * Pi]         */
	double TGD[2];              /*  [s]     [s1w7b17-24]              */
}EPHN_t;

struct IGS_Clk_t
{
	gtime_t   time;
	u1        prn[MAXSAT];
	double    SatClk[MAXSAT];
};
struct IGS_Orbit_t
{
	gtime_t  time;
	u1       prn[MAXSAT];
	double   Orbit[MAXSAT][3];
};

struct GNSSBrdc_t
{
	GNSSBrdc_t()
	{
		n=max=0;
		eph = NULL;
	}
	u2 n,max;
	EPHN_t *eph;
};

typedef struct 
{        /* antenna parameter type */
	int sat;            /* satellite number (0:receiver) */
	char type[MAXANT];  /* antenna type */
	char code[MAXANT];  /* serial number or satellite code */
	gtime_t ts,te;      /* valid time start and end */
	double off[NEFREQ][ 3]; /* phase center offset e/n/u or x/y/z (m) */
	double var[NEFREQ][19]; /* phase center variation (m) */
	/* el=90,85,...,0 or nadir=0,1,2,3,... (deg) */
} Pcv_t;

struct GNSSIGS_t
{
	GNSSIGS_t()
	{
		maxclk=maxsp3=nsp3=nclk=0;
		clk = NULL;
		sp3 = NULL;
	}
	u2 nclk,maxclk,nsp3,maxsp3;
	u2 SamplRate[2];     //采样率：轨道  钟差
	IGS_Clk_t    *clk;
	IGS_Orbit_t  *sp3;
};

struct GNSSEPH_t
{
	GNSSBrdc_t brdc[MAXSYS];
	GNSSIGS_t  IGS;
};


typedef struct {        /* earth rotation parameter data type */
	double mjd;         /* mjd (days) */
	double xp,yp;       /* pole offset (rad) */
	double xpr,ypr;     /* pole offset rate (rad/day) */
	double ut1_utc;     /* ut1-utc (s) */
	double lod;         /* length of day (s/day) */
} erpd_t;


typedef struct {        /* TEC grid type */
	gtime_t time;       /* epoch time (GPST) */
	int ndata[3];       /* TEC grid data size {nlat,nlon,nhgt} */
	double rb;          /* earth radius (km) */
	double lats[3];     /* latitude start/interval (deg) */
	double lons[3];     /* longitude start/interval (deg) */
	double hgts[3];     /* heights start/interval (km) */
	double *data;       /* TEC grid data (tecu) */
	float *rms;         /* RMS values (tecu) */
} tec_t;

typedef struct
{
	tec_t  *tec;
	int    nt;
	int    nmax;
}tec_tt;

typedef struct {        /* earth rotation parameter type */
	int n,nmax;         /* number and max number of data */
	erpd_t *data;       /* earth rotation parameter data */
} erp_t;

struct SSat_t
{
	double    Gf[4];
	double    CSMW[6];
	double    Mw[2];
	double    Var[2];
	int       Number[2];
	int       FixWL[2];
	int       FixFlag[2];
	int       FixCount;
	int       TrackNum;
	int       CSNum[3];
	bool      Flag;
	bool      CSlip[3];
	/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	double    DetaDCB;                            /*     计算DCB值     */
	double    DCBRMS;
	int       DCBNum;                             /*参与计算DCB历元个数*/
	/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	double    FCB[2];
	int       FCBNum[2];
	/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	int       SoD;
};

struct M_FCB
{
	double FCB[2];
	int    FCBNum[2];
	int    Prn;
};
struct M_Station
{
	char          Name[10];
	double        XYZ[3];
	double        Pos[3];
	u2            _NO;
	FILE          *ft_Result;
	FILE          *ft_mw,*ft_ewl;
	PROCESS_MODE  Pro_Mode;
	int           RecType;
};
struct FileDir
{
	char           ProDir[MAXLENGTH];   //directory of home address
	char           NavDir[MAXLENGTH];   //directory of Nav address
	char           ObsDir[MAXLENGTH];   //directory of Obs address
	char           ObsRefDir[MAXLENGTH];
	char           Sp3Dir[MAXLENGTH];   //directory of Sp3 address
	char           ClkDir[MAXLENGTH];   //directory of Clk address
	char           TableDir[MAXLENGTH]; //directory of Table address
	char           DCBDir[MAXLENGTH];   //directory of DCB address
	char           TECDir[MAXLENGTH];
	char           MassDir[MAXLENGTH];
};
struct ProcessOption
{
	FileDir        Directory;
	M_Station      *Station;
	TTime          StartTime;
	TTime          EndTime;
	u2             StartMJDN;
	u2             EndMJDN;
	int            StartSoD;
	int            EndSoD;
	u2             StaNum;           // total number of station	
	u2             Mode;
	double         Interval;
	u2             NavSys;
	FILE           *ft_Static[4];
	int            ProcessDoY;
	int            ProcessYear;
	bool           HighRate;
	ProcessOption()
	{
		StaNum = 0;
		Station = NULL;
	}
};
struct WL_Info
{
	int Prn;
	int SSod;
	int ESod;
	int EWL;
	int WL;
};
struct MP_Info
{
	double *MP_1;                   //MP Value
	double *MP_2;
	double Std[8];
	int    LastEpoch[2];
};
struct AMB_Info
{
	double *EWL[MAXSAT];
	double *WL[MAXSAT]; 
	double *Var[MAXSAT];
	int    Numb[MAXSAT];      // raw amb number
//----------------------------------------------------
	double FCB[MAXSAT][2];    //FCB EWL/WL
	double P_DCB[MAXSAT][2];
//----------------------------------------------------
	double Amb[MAXSAT][2];   //fixed amb value EWL/WL
	bool   Flag[MAXSAT][2];  //fixed amb FLAG EWL/WL
	int    FixNumber[3][4];     //fixed amb number EWL/WL
	double RcvFcb[MAXSYS][2];
	double SatFloAmb[MAXSAT][2];
	WL_Info *WLInfo;
	int     NWLInfo[2];
//------------------------------------------------
	double DCBValue[MAXSAT][2];
	int    DCBNumb[MAXSAT];
	//
	bool   AmbProcess;
};
struct GNSSInfo
{
	double        *m_X;
	double        *m_P;
	u2            Nx;
	u2            Nf;
	SSat_t        SSat[MAXSAT];          //Satellite Information
	gtime_t       LastObsTime;
	M_Station     *Station;   //Station Information
	PROCESS_MODE  Mode;
	u2            N_Base;
	u2            N_Rove;
	bool          lastepoch;
	int           LogSoD;
	/*     -----------------    TEQC   -----------------  */
	MP_Info       MP_S[MAXSAT];
	int           EpoNum;
	int           UseFulNum;
	/*     -----------------    TEQC   -----------------  */
	double        Ave[3];
	double        RMS[3];
	int           Numb;
	int           Visible[MAXSAT][2];
	/*++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	M_FCB         *Fcb;
	M_FCB         *FcbDay;
	int           Nfcb[2];
	int           NfcbDay[2];
	double        UDFcb[MAXSAT][2];
	bool          UDFcbFlag[MAXSAT][2];
	bool          bUPDExclude[MAXSAT][2];
	/*++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	double        dPIFGF[MAXSAT];
	int           iPIFGF[MAXSAT];
};

struct UPDInfo
{
	double _EWLRefFCB[MAXSAT];
	double _WLRefFCB[MAXSAT];
	bool   _bEWLRef[MAXSAT];
	bool   _bWLRef[MAXSAT];
	double  *_dEWLSeri;
	double  *_dWLSeri;
	bool    *_bEWLSeri;
	bool    *_bWLSeri;
	int      _nDay;
	FILE    *_ftRes;
};

struct GNSSCodeBias
{
	double    P1C1[MAXSAT];
	double    P2C2[MAXSAT];
	double    P3C3[MAXSAT];
	double    BDSRecBias[10][MAXSAT][3];
	int       nRecType;
};
#endif