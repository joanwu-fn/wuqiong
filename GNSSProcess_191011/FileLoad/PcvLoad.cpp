#include "../Common/FileLoad.h"

/* decode antenna parameter field --------------------------------------------*/
int decodef(char *p, int n, double *v)
{
	int i;

	for (i=0;i<n;i++) v[i]=0.0;
	for (i=0,p=strtok(p," ");p&&i<n;p=strtok(NULL," ")) {
		v[i++]=atof(p)*1E-3;
	}
	return i;
}
/* add antenna parameter -----------------------------------------------------*/
void addpcv(const Pcv_t PcvT, Pcv_t *PcvS)
{
	if(PcvT.sat <= 0 || PcvT.sat > MAXSAT)
		return;
	int prn_1 = PcvT.sat -1;
	double diff = timediff(PcvT.ts, PcvS[prn_1].ts);

	if(diff > 0.0)
	{
		memcpy(&PcvS[prn_1], &PcvT, sizeof(Pcv_t));
	}
	return;
}
/* read ngs antenna parameter file -------------------------------------------*/
int readngspcv(const char *file, Pcv_t *pcvs)
{
	FILE *fp;
	static const Pcv_t pcv0={0};
	Pcv_t pcv;
	double neu[3];
	int n=0;
	char buff[256];

	if (!(fp=fopen(file,"r"))) {
		return 0;
	}
	while (fgets(buff,sizeof(buff),fp)) {

		if (strlen(buff)>=62&&buff[61]=='|') continue;

		if (buff[0]!=' ') n=0; /* start line */
		if (++n==1) {
			pcv=pcv0;
			strncpy(pcv.type,buff,61); pcv.type[61]='\0';
		}
		else if (n==2) {
			if (decodef(buff,3,neu)<3) continue;
			pcv.off[0][0]=neu[1];
			pcv.off[0][1]=neu[0];
			pcv.off[0][2]=neu[2];
		}
		else if (n==3) decodef(buff,10,pcv.var[0]);
		else if (n==4) decodef(buff,9,pcv.var[0]+10);
		else if (n==5) {
			if (decodef(buff,3,neu)<3) continue;;
			pcv.off[1][0]=neu[1];
			pcv.off[1][1]=neu[0];
			pcv.off[1][2]=neu[2];
		}
		else if (n==6) decodef(buff,10,pcv.var[1]);
		else if (n==7) {
			decodef(buff,9,pcv.var[1]+10);
			addpcv(pcv,pcvs);
		}
	}
	fclose(fp);

	return 1;
}
/* read antex file ----------------------------------------------------------*/
int readantex(const char *file, Pcv_t *pcvs)
{
	FILE *fp;
	static const Pcv_t pcv0={0};
	Pcv_t pcv;
	double neu[3];
	int i,f,freq=0,state=0,freqs[]={1,2,5,6,7,8,0};
	char buff[256];

	if ((fp=fopen(file,"r")) == NULL) 
	{
		return 0;
	}
	while (fgets(buff,sizeof(buff),fp)) {

		if (strlen(buff)<60||strstr(buff+60,"COMMENT")) continue;

		if (strstr(buff+60,"START OF ANTENNA")) {
			pcv=pcv0;
			state=1;
		}
		if (strstr(buff+60,"END OF ANTENNA")) {
			addpcv(pcv,pcvs);
			state=0;
		}
		if (!state) continue;

		if (strstr(buff+60,"TYPE / SERIAL NO")) {
			strncpy(pcv.type,buff   ,20); pcv.type[20]='\0';
			strncpy(pcv.code,buff+20,20); pcv.code[20]='\0';
			if (!strncmp(pcv.code+3,"        ",8)) 
			{
				pcv.sat=Satid2No(pcv.code);
			}
		}
		else if (strstr(buff+60,"VALID FROM")) {
			if (!str2time(buff,0,43,&pcv.ts)) continue;
		}
		else if (strstr(buff+60,"VALID UNTIL")) {
			if (!str2time(buff,0,43,&pcv.te)) continue;
		}
		else if (strstr(buff+60,"START OF FREQUENCY")) {
			if (sscanf(buff+4,"%d",&f)<1) continue;
			for (i=0;i<NEFREQ;i++) if (freqs[i]==f) break;
			if (i<NEFREQ) freq=i+1;
		}
		else if (strstr(buff+60,"END OF FREQUENCY")) {
			freq=0;
		}
		else if (strstr(buff+60,"NORTH / EAST / UP")) {
			if (freq<1||NEFREQ<freq) continue;
			if (decodef(buff,3,neu)<3) continue;
			pcv.off[freq-1][0]=neu[pcv.sat?0:1]; /* x or e */
			pcv.off[freq-1][1]=neu[pcv.sat?1:0]; /* y or n */
			pcv.off[freq-1][2]=neu[2];           /* z or u */
		}
		else if (strstr(buff,"NOAZI")) {
			if (freq<1||NEFREQ<freq) continue;
			if ((i=decodef(buff+8,19,pcv.var[freq-1]))<=0) continue;
			for (;i<19;i++) pcv.var[freq-1][i]=pcv.var[freq-1][i-1];
		}
	}
	fclose(fp);

	return 1;
}
