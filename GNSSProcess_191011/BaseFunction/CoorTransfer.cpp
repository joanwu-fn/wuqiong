#include "../Common/BaseFunction.h"

/* transform ecef to geodetic postion ------------------------------------------
* transform ecef position to geodetic position
* args   : double *r        I   ecef position {x,y,z} (m)
*          double *pos      O   geodetic position {lat,lon,h} (rad,m)
* return : none
* notes  : WGS84, ellipsoidal height
*-----------------------------------------------------------------------------*/
void ecef2pos(const double *r, double *pos)
{
	double e2=FE_WGS84*(2.0-FE_WGS84),r2=dot(r,r,2),z,zk,v=RE_WGS84,sinp;

	for (z=r[2],zk=0.0;fabs(z-zk)>=1E-4;) {
		zk=z;
		sinp=z/sqrt(r2+z*z);
		v=RE_WGS84/sqrt(1.0-e2*sinp*sinp);
		z=r[2]+v*e2*sinp;
	}
	pos[0]=r2>1E-12?atan(z/sqrt(r2)):(r[2]>0.0?PI/2.0:-PI/2.0);
	pos[1]=r2>1E-12?atan2(r[1],r[0]):0.0;
	pos[2]=sqrt(r2+z*z)-v;
}
/* ecef to local coordinate transfromation matrix ------------------------------
* compute ecef to local coordinate transfromation matrix
* args   : double *pos      I   geodetic position {lat,lon} (rad)
*          double *E        O   ecef to local coord transformation matrix (3x3)
* return : none
* notes  : matirix stored by column-major order (fortran convention)
*-----------------------------------------------------------------------------*/
void xyz2enu(const double *pos, double *E)
{
	double sinp=sin(pos[0]),cosp=cos(pos[0]),sinl=sin(pos[1]),cosl=cos(pos[1]);

	E[0]=-sinl;      E[3]=cosl;       E[6]=0.0;
	E[1]=-sinp*cosl; E[4]=-sinp*sinl; E[7]=cosp;
	E[2]=cosp*cosl;  E[5]=cosp*sinl;  E[8]=sinp;
}
/* transform ecef vector to local tangental coordinate -------------------------
* transform ecef vector to local tangental coordinate
* args   : double *pos      I   geodetic position {lat,lon} (rad)
*          double *r        I   vector in ecef coordinate {x,y,z}
*          double *e        O   vector in local tangental coordinate {e,n,u}
* return : none
*-----------------------------------------------------------------------------*/
void ecef2enu(const double *pos, const double *r, double *e)
{
	double E[9];

	xyz2enu(pos,E);
	matmul("NN",3,1,3,1.0,E,r,0.0,e);
}