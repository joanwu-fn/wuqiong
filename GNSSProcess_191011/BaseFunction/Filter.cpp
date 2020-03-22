#include "../Common/BaseFunction.h"

#include <stdarg.h>

/* fatal error ---------------------------------------------------------------*/
void fatalerr(const char *format, ...)
{
	va_list ap;
	va_start(ap,format); vfprintf(stderr,format,ap); va_end(ap);
	exit(-9);
}

/* copy matrix -----------------------------------------------------------------
* copy matrix
* args   : double *A        O   destination matrix A (n x m)
*          double *B        I   source matrix B (n x m)
*          int    n,m       I   number of rows and columns of matrix
* return : none
*-----------------------------------------------------------------------------*/
void matcpy(double *A, const double *B, int n, int m)
{
	memcpy(A,B,sizeof(double)*n*m);
}
/* new matrix ------------------------------------------------------------------
* allocate memory of matrix 
* args   : int    n,m       I   number of rows and columns of matrix
* return : matrix pointer (if n<=0 or m<=0, return NULL)
*-----------------------------------------------------------------------------*/
double *mat(int n, int m)
{
	double *p;

	if (n<=0||m<=0) return NULL;
	if (!(p=(double *)malloc(sizeof(double)*n*m))) {
		fatalerr("matrix memory allocation error: n=%d,m=%d\n",n,m);
	}
	return p;
}
/* new integer matrix ----------------------------------------------------------
* allocate memory of integer matrix 
* args   : int    n,m       I   number of rows and columns of matrix
* return : matrix pointer (if n<=0 or m<=0, return NULL)
*-----------------------------------------------------------------------------*/
int *imat(int n, int m)
{
	int *p;

	if (n<=0||m<=0) return NULL;
	if (!(p=(int *)malloc(sizeof(int)*n*m))) {
		fatalerr("integer matrix memory allocation error: n=%d,m=%d\n",n,m);
	}
	return p;
}
double *zeros(int n, int m)
{
	double *p;
	if(n<=0 || m<=0)
	{
		return NULL;
	}
	if((p = (double *)calloc(sizeof(double), n*m)) == NULL)
	{
		printf("matrix memory allocation error:n=%d,m=%d\n",n,m);
	}
	return p;
}
/* identity matrix -------------------------------------------------------------
* generate new identity matrix
* args   : int    n         I   number of rows and columns of matrix
* return : matrix pointer (if n<=0, return NULL)
*-----------------------------------------------------------------------------*/
double *eye(int n)
{
	double *p;
	int i;

	if ((p=zeros(n,n))) for (i=0;i<n;i++) p[i+i*n]=1.0;
	return p;
}
/* LU decomposition ----------------------------------------------------------*/
int ludcmp(double *A, int n, int *indx, double *d)
{
	double big,s,tmp,*vv=mat(n,1);
	int i,imax=0,j,k;

	*d=1.0;
	for (i=0;i<n;i++) {
		big=0.0; for (j=0;j<n;j++) if ((tmp=fabs(A[i+j*n]))>big) big=tmp;
		if (big>0.0) vv[i]=1.0/big; else {free(vv); return -1;}
	}
	for (j=0;j<n;j++) {
		for (i=0;i<j;i++) {
			s=A[i+j*n]; for (k=0;k<i;k++) s-=A[i+k*n]*A[k+j*n]; A[i+j*n]=s;
		}
		big=0.0;
		for (i=j;i<n;i++) {
			s=A[i+j*n]; for (k=0;k<j;k++) s-=A[i+k*n]*A[k+j*n]; A[i+j*n]=s;
			if ((tmp=vv[i]*fabs(s))>=big) {big=tmp; imax=i;}
		}
		if (j!=imax) {
			for (k=0;k<n;k++) {
				tmp=A[imax+k*n]; A[imax+k*n]=A[j+k*n]; A[j+k*n]=tmp;
			}
			*d=-(*d); vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (A[j+j*n]==0.0) {free(vv); return -1;}
		if (j!=n-1) {
			tmp=1.0/A[j+j*n]; for (i=j+1;i<n;i++) A[i+j*n]*=tmp;
		}
	}
	free(vv);
	return 0;
}
/* LU back-substitution ------------------------------------------------------*/
void lubksb(const double *A, int n, const int *indx, double *b)
{
	double s;
	int i,ii=-1,ip,j;

	for (i=0;i<n;i++) {
		ip=indx[i]; s=b[ip]; b[ip]=b[i];
		if (ii>=0) for (j=ii;j<i;j++) s-=A[i+j*n]*b[j]; else if (s) ii=i;
		b[i]=s;
	}
	for (i=n-1;i>=0;i--) {
		s=b[i]; for (j=i+1;j<n;j++) s-=A[i+j*n]*b[j]; b[i]=s/A[i+i*n];
	}
}
/* inverse of matrix ---------------------------------------------------------*/
int matinv(double *A, int n)
{
	double d,*B;
	int i,j,*indx;

	indx=imat(n,1); B=mat(n,n); matcpy(B,A,n,n);
	if (ludcmp(B,n,indx,&d)) {free(indx); free(B); return -1;}
	for (j=0;j<n;j++) {
		for (i=0;i<n;i++) A[i+j*n]=0.0; A[j+j*n]=1.0;
		lubksb(B,n,indx,A+j*n);
	}
	free(indx); free(B);
	return 0;
}
/* multiply matrix -----------------------------------------------------------*/
void matmul(const char *tr, int n, int k, int m, double alpha,
			const double *A, const double *B, double beta, double *C)
{
	double d;
	int i,j,x,f=tr[0]=='N'?(tr[1]=='N'?1:2):(tr[1]=='N'?3:4);

	for (i=0;i<n;i++) for (j=0;j<k;j++) {
		d=0.0;
		switch (f) {
			case 1: for (x=0;x<m;x++) d+=A[i+x*n]*B[x+j*m]; break;
			case 2: for (x=0;x<m;x++) d+=A[i+x*n]*B[j+x*k]; break;
			case 3: for (x=0;x<m;x++) d+=A[x+i*m]*B[x+j*m]; break;
			case 4: for (x=0;x<m;x++) d+=A[x+i*m]*B[j+x*k]; break;
		}
		if (beta==0.0) C[i+j*n]=alpha*d; else C[i+j*n]=alpha*d+beta*C[i+j*n];
	}
	i=0;
}
/* least square estimation -----------------------------------------------------
* least square estimation by solving normal equation (x=(A*A')^-1*A*y)
* args   : double *A        I   transpose of (weighted) design matrix (n x m)
*          double *y        I   (weighted) measurements (m x 1)
*          int    n,m       I   number of parameters and measurements (n<=m)
*          double *x        O   estmated parameters (n x 1)
*          double *Q        O   esimated parameters covariance matrix (n x n)
* return : status (0:ok,0>:error)
* notes  : for weighted least square, replace A and y by A*w and w*y (w=W^(1/2))
*          matirix stored by column-major order (fortran convention)
*-----------------------------------------------------------------------------*/
int lsq(const double *A, const double *y, int n, int m, double *x,
		double *Q)
{
	double *Ay;
	int info;

	if (m<n) return -1;
	Ay=zeros(n,1);
	matmul("NN",n,1,m,1.0,A,y,0.0,Ay); /* Ay=A*y */
	matmul("NT",n,n,m,1.0,A,A,0.0,Q);  /* Q=A*A' */
	if (!(info=matinv(Q,n))) 
	{
		matmul("NN",n,1,n,1.0,Q,Ay,0.0,x); /* x=Q^-1*Ay */
	}
	free(Ay);
	return info;
}
/* kalman filter ---------------------------------------------------------------
* kalman filter state update as follows:
*
*   K=P*H*(H'*P*H+R)^-1, xp=x+K*v, Pp=(I-K*H')*P
*
* args   : double *x        I   states vector (n x 1)
*          double *P        I   covariance matrix of states (n x n)
*          double *H        I   transpose of design matrix (n x m)
*          double *v        I   innovation (measurement - model) (m x 1)
*          double *R        I   covariance matrix of measurement error (m x m)
*          int    n,m       I   number of states and measurements
*          double *xp       O   states vector after update (n x 1)
*          double *Pp       O   covariance matrix of states after update (n x n)
* return : status (0:ok,<0:error)
* notes  : matirix stored by column-major order (fortran convention)
*          if state x[i]==0.0, not updates state x[i]/P[i+i*n]
*-----------------------------------------------------------------------------*/
int filter_(const double *x, const double *P, const double *H,
				   const double *v, const double *R, int n, int m,
				   double *xp, double *Pp)
{
	double *F=mat(n,m),*Q=mat(m,m),*K=mat(n,m),*I=eye(n);
	int info;
	double *Q1=mat(m,m),*l=mat(m,1);double *xv=mat(n,1);

	matcpy(Q,R,m,m);
	matcpy(xp,x,n,1);

	matcpy(xv,x,n,1);//////////////////////

	matmul("NN",n,m,n,1.0,P,H,0.0,F);       /* Q=H'*P*H+R */
	matmul("TN",m,m,n,1.0,H,F,1.0,Q);
	if (!(info=matinv(Q,m))) {
		matmul("NN",n,m,m,1.0,F,Q,0.0,K);   /* K=P*H*Q^-1 */
		matmul("NN",n,1,m,1.0,K,v,1.0,xp);  /* xp=x+K*v */
		matmul("NT",n,n,m,-1.0,K,H,1.0,I);  /* Pp=(I-K*H')*P */
		matmul("NN",n,n,n,1.0,I,P,0.0,Pp);
	}
	free(F); free(Q); free(K); free(I);

	free(Q1);free(l);free(xv);
	return info;
}

int printMatfp(const double *A,int m,int n,FILE *fp)
{
	int i,j;
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			fprintf(fp,"%8.3f",A[j+i*n]);
		}
		fprintf(fp,"\n");
	}
	return 1;
}
int filter(double *x, double *P, const double *H, const double *v,
		   const double *R, int n, int m)
{
	double *x_,*xp_,*P_,*Pp_,*H_;
	int i,j,k,info,*ix;
// 	FILE *fp=fopen("d:\\0.txt","w");
	ix=imat(n,1); for (i=k=0;i<n;i++) if (x[i]!=0.0&&P[i+i*n]>0.0) ix[k++]=i;
	x_=mat(k,1); xp_=mat(k,1); P_=mat(k,k); Pp_=mat(k,k); H_=mat(k,m);
	for (i=0;i<k;i++) 
	{
		x_[i]=x[ix[i]];
		for (j=0;j<k;j++) P_[i+j*k]=P[ix[i]+ix[j]*n];
		for (j=0;j<m;j++) H_[i+j*k]=H[ix[i]+j*n];
	}
// 	fprintf(fp,"   --------------------   H  ------------------------\n");
// 	printMatfp(H_,m,k,fp);
// 	fprintf(fp,"   --------------------   V  ------------------------\n");
// 	printMatfp(v,1,m,fp);
// 	fprintf(fp,"   --------------------   R  ------------------------\n");
// 	printMatfp(R,m,m,fp);
// 	fprintf(fp,"   --------------------   X  ------------------------\n");
// 	printMatfp(x_,1,k,fp);
// 	fprintf(fp,"   --------------------   before filter P  ------------------------\n");
// 	printMatfp(P_,k,k,fp);
	info=filter_(x_,P_,H_,v,R,k,m,xp_,Pp_);
	for (i=0;i<k;i++) {
		x[ix[i]]=xp_[i];
		for (j=0;j<k;j++) P[ix[i]+ix[j]*n]=Pp_[i+j*k];
	}
// 	fprintf(fp,"   ----------------       after filter  X    --------------------\n");
// 	printMatfp(xp_,1,k,fp);
// 	fprintf(fp,"   --------------------   after filter P  ------------------------\n");
// 	printMatfp(Pp_,k,k,fp);
	free(ix); free(x_); free(xp_); free(P_); free(Pp_); free(H_);
// 	fclose(fp);
	return info;
}