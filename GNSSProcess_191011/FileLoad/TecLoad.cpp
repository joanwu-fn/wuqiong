#include "../Common/FileLoad.h"

/* get index -----------------------------------------------------------------*/
int getindex(double value, const double *range)
{
	if (range[2]==0.0) return 0;
	if (range[1]>0.0&&(value<range[0]||range[1]<value)) return -1;
	if (range[1]<0.0&&(value<range[1]||range[0]<value)) return -1;
	return (int)floor((value-range[0])/range[2]+0.5);
}
/* get number of items -------------------------------------------------------*/
int nitem(const double *range)
{
	return getindex(range[1],range)+1;
}
/* data index (i:lat,j:lon,k:hgt) --------------------------------------------*/
int dataindex(int i, int j, int k, const int *ndata)
{
	if (i<0||ndata[0]<=i||j<0||ndata[1]<=j||k<0||ndata[2]<=k) return -1;
	return i+ndata[0]*(j+ndata[1]*k);
}
/* add tec data to navigation data -------------------------------------------*/
tec_t *addtec(const double *lats, const double *lons, const double *hgts,
					 double rb, tec_tt *tec)
{
	tec_t *p,*nav_tec;
	gtime_t time0={0};
	int i,n,ndata[3];

	ndata[0]=nitem(lats);
	ndata[1]=nitem(lons);
	ndata[2]=nitem(hgts);
	if (ndata[0]<=1||ndata[1]<=1||ndata[2]<=0)
	{
		return NULL;
	}

	
	if (tec->nt >= tec->nmax) {
		tec->nmax += 256;
		if (!(nav_tec=(tec_t *)realloc(tec->tec,sizeof(tec_t)*tec->nmax))) 
		{
			free(tec); tec=NULL; tec->nt = tec->nmax=0;
			return NULL;
		}
		tec->tec=nav_tec;
	}
	p = tec->tec + tec->nt;
	p->time=time0;
	p->rb=rb;
	for (i=0;i<3;i++) 
	{
		p->ndata[i]=ndata[i];
		p->lats[i]=lats[i];
		p->lons[i]=lons[i];
		p->hgts[i]=hgts[i];
	}
	n=ndata[0]*ndata[1]*ndata[2];

	if (!(p->data= (double *)malloc(sizeof(double)*n))||
		!(p->rms = (float *)malloc(sizeof(float )*n))) 
	{
		return NULL;
	}
	for (i=0;i<n;i++) 
	{
		p->data[i]=0.0;
		p->rms [i]=0.0f;
	}
	tec->nt++;
	return p;
}
/* read ionex header ---------------------------------------------------------*/
double readionexh(FILE *fp, double *lats, double *lons, double *hgts,
						 double *rb, double *nexp)
{
	double ver=0.0;
	char buff[1024],*label;

	while (fgets(buff,sizeof(buff),fp)) {

		if (strlen(buff)<60) continue;
		label=buff+60;

		if (strstr(label,"IONEX VERSION / TYPE")==label) {
			if (buff[20]=='I') ver=str2num(buff,0,8);
		}
		else if (strstr(label,"BASE RADIUS")==label) {
			*rb=str2num(buff,0,8);
		}
		else if (strstr(label,"HGT1 / HGT2 / DHGT")==label) {
			hgts[0]=str2num(buff, 2,6);
			hgts[1]=str2num(buff, 8,6);
			hgts[2]=str2num(buff,14,6);
		}
		else if (strstr(label,"LAT1 / LAT2 / DLAT")==label) {
			lats[0]=str2num(buff, 2,6);
			lats[1]=str2num(buff, 8,6);
			lats[2]=str2num(buff,14,6);
		}
		else if (strstr(label,"LON1 / LON2 / DLON")==label) {
			lons[0]=str2num(buff, 2,6);
			lons[1]=str2num(buff, 8,6);
			lons[2]=str2num(buff,14,6);
		}
		else if (strstr(label,"EXPONENT")==label) {
			*nexp=str2num(buff,0,6);
		}
		else if (strstr(label,"END OF HEADER")==label) {
			return ver;
		}
	}
	return 0.0;
}
/* read ionex body -----------------------------------------------------------*/
int readionexb(FILE *fp, const double *lats, const double *lons,
					  const double *hgts, double rb, double nexp, tec_tt *tec)
{
	tec_t *p=NULL;
	gtime_t time={0};
	double lat,lon[3],hgt,x;
	int i,j,k,n,m,index,type=0;
	char buff[1024],*label=buff+60;

	while (fgets(buff,sizeof(buff),fp)) {

		if (strlen(buff)<60) continue;

		if (strstr(label,"START OF TEC MAP")==label) {
			if ((p=addtec(lats,lons,hgts,rb, tec))) 
				type=1;
		}
		else if (strstr(label,"END OF TEC MAP")==label) {
			type=0;
			p=NULL;
		}
		else if (strstr(label,"START OF RMS MAP")==label) {
			type=2;
			p=NULL;
		}
		else if (strstr(label,"END OF RMS MAP")==label) {
			type=0;
			p=NULL;
		}
		else if (strstr(label,"EPOCH OF CURRENT MAP")==label) {
			if (str2time(buff,0,36,&time)) 
			{
				continue;
			}
			if (type==2) {
				for (i=tec->nt-1;i>=0;i--) 
				{
					if (fabs(timediff(time,tec->tec[i].time))>=1.0) 
					{
						continue;
					}
					p = tec->tec+i;
					break;
				}
			}
			else if (p) p->time=time;
		}
		else if (strstr(label,"LAT/LON1/LON2/DLON/H")==label&&p) {
			lat   =str2num(buff, 2,6);
			lon[0]=str2num(buff, 8,6);
			lon[1]=str2num(buff,14,6);
			lon[2]=str2num(buff,20,6);
			hgt   =str2num(buff,26,6);

			i=getindex(lat,p->lats);
			k=getindex(hgt,p->hgts);
			n=nitem(lon);

			for (m=0;m<n;m++) {
				if (m%16==0&&!fgets(buff,sizeof(buff),fp))
				{
					break;
				}
				j=getindex(lon[0]+lon[2]*m,p->lons);
				if ((index=dataindex(i,j,k,p->ndata))<0) 
				{
					continue;
				}
				if ((x=str2num(buff,m%16*5,5))==9999.0) 
				{
					continue;
				}

				if (type==1) 
				{
					p->data[index]=x*pow(10.0,nexp);
				}
				else 
				{
					p->rms[index]=(float)(x*pow(10.0,nexp));
				}
			}
		}
	}
	return 1;
}
/* combine tec grid data -----------------------------------------------------*/
void combtec(tec_tt *tec)
{
	tec_t tmp;
	int i,j,n=0;

	for (i = 0; i < tec->nt-1; i++) 
	{
		for (j= i+1; j < tec->nt; j++) 
		{
			if (timediff(tec->tec[j].time,tec->tec[i].time)<0.0) 
			{
				tmp = tec->tec[i];
				tec->tec[i] = tec->tec[j];
				tec->tec[j] = tmp;
			}
		}
	}
	for (i=0;i < tec->nt;i++) 
	{
		if (i>0&&timediff(tec->tec[i].time,tec->tec[n-1].time)==0.0) {
			free(tec->tec[n-1].data);
			free(tec->tec[n-1].rms );
			tec[n-1]=tec[i];
			continue;
		}
		tec->tec[n++] = tec->tec[i];
	}
	tec->nt = n;
}

bool ReadTecFile(char *file, tec_tt * tec)
{
	FILE *fp = NULL;
	double lats[3]={0},lons[3]={0},hgts[3]={0},rb=0.0,nexp=-1.0;
	int i,n;
	fp = fopen(file, "r");
	tec->nt = 0;
	tec->nmax = 0;
	if(fp == NULL)
	{
		return false;
	}
	if (readionexh(fp,lats,lons,hgts,&rb,&nexp)<=0.0) 
	{
		fclose(fp);
		return false;
	}
	/* read ionex body */
	readionexb(fp,lats,lons,hgts,rb,nexp,tec);
	fclose(fp);

	/* combine tec grid data */
	if (tec->nt>0) 
	{
		combtec(tec);
	}
	

}