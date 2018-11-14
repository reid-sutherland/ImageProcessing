/* Mathematical routines */

#include "levelutil.h"


/* LU decomposition */
void ludcmp(double **a, int n, int *indx, double *d)
{
  int i,imax,j,k;
  double big,dum,sum,temp;
  double *vv;
  
  vv=(double*)malloc(n*sizeof(double));
  if(vv==NULL)
    mexErrMsgTxt("LU Decomp , Memory allocation error");

  *d=1.0;
  for (i=0;i<n;i++) {
    big=0.0;
    for (j=0;j<n;j++)
      if ((temp=fabs(a[i][j])) > big) 
	big=temp;
    if (big == 0.0)
      mexErrMsgTxt("LU Decomp , Singular matrix");
    vv[i]=1.0/big;
  }
  for (j=0;j<n;j++) {
    for (i=0;i<j;i++) {
      sum=a[i][j];
      for (k=0;k<i;k++) 
	sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for (i=j;i<n;i++) {
      sum=a[i][j];
      for (k=0;k<j;k++)
	sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {
	big=dum;
	imax=i;
      }
    }
    if (j != imax) {
      for (k=0;k<n;k++) {
	dum=a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) 
      a[j][j]=TINY;
    if (j != (n-1)) {
      dum=1.0/(a[j][j]);
      for (i=j+1;i<n;i++) 
	a[i][j] *= dum;
    }
  }
  free(vv);
}
/* LU backward substitution */
void lubksb(double **a, int n, int *indx, double b[])
{
  int i,ii=(-1),ip,j;
  double sum;
  
  for (i=0;i<n;i++) {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii>=0)
      for (j=ii;j<=i-1;j++) 
	sum -= a[i][j]*b[j];
    else 
      if (sum) 
	ii=i;
    b[i]=sum;
  }
  for (i=n-1;i>=0;i--) {
    sum=b[i];
    for (j=i+1;j<n;j++) 
      sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}
/* Determinant of inverse of a matrix */
double detinv(double **a, int n){
  
  double **y;
  double *col, d=1.0;
  int i,j, *indx;

  if((indx = (int *)malloc(n*sizeof(int)))==NULL)
    mexErrMsgTxt("Gauss Vect , Memory allocation error");
  if((col = (double *)malloc(n*sizeof(double)))==NULL)
    mexErrMsgTxt("Gauss Vect , Memory allocation error"); 
  if((y = (double **)malloc(n*sizeof(double*)))==NULL)
    mexErrMsgTxt("Gauss Vect , Memory allocation error"); 
  for(i=0;i<n;i++)
    if((y[i] = (double *)malloc(n*sizeof(double)))==NULL)
      mexErrMsgTxt("Gauss Vect , Memory allocation error"); 
  
  ludcmp(a,n,indx,&d);
  for(i=0;i<n;i++)
    d *= a[i][i];

  for(j=0;j<n;j++){
    for(i=0;i<n;i++)
      col[i] = 0.0;
    col[j] = 1.0;
    lubksb(a,n,indx,col);
    for(i=0;i<n;i++)
      y[i][j] = col[i];
  }
  
  for(i=0;i<n;i++){
    memcpy(a[i],y[i],n*sizeof(double));
    free(y[i]);
  }
  
  free(y);
  free(col);
  free(indx);
  
  return(d);

}

/* Function to compute uu and vv from vector of parameters in
   the vectorial case */
double cal_uu(double *I, double *param, int N){
  
  double uu, uu_sum, uu_m, uu_v, uu_sumi, uu_mi;
  double uu_d;
  int i, j;
  

  uu_sum = 0.0;
  for(i=0;i<N;i++){
    uu_m = param[N*N+N+1+i];
    uu_sumi = 0.0;
    for(j=0;j<N;j++){
      uu_mi = param[N*N+N+1+j];
      uu_v = param[N*N+2*N+1+i*N+j];
      uu_sumi += (I[j]-uu_mi)*uu_v;
    }
    uu_sum += (I[i]-uu_m)*uu_sumi;
  }

  uu_d =  param[2*N*N+2*N+1]; 
  uu = (log(uu_d) + uu_sum)/2.0; 
  
  return(uu);
}


/* Function to update the parameters (mean and covariance) */
void update_param(double *param,int N){
  
  int j,l;
  double **temp;
  
  if((temp=(double**)malloc(N*sizeof(double*)))==NULL)
    mexErrMsgTxt("Gauss Vect , Memory allocation error");
  
  /* means */
  for(j=0;j<N;j++){
    param[N*(N+1)+1 + j] = param[j+1]/param[0];
  }
  /* covariance = Int. II^T/A - mean*mean^T */
  for(j=0;j<N;j++){
    for(l=0;l<N;l++){
      param[N*(N+2)+1 + j*N+l] = param[(N+1) + j*N+l] / param[0] -
	param[N*(N+1)+1 + j] * param[N*(N+1)+1 + l];
    }
  } 
  
  /* Compute determinant */
  /* Set temporal matrix */
  for(j=0;j<N;j++)
    temp[j] = param + N*N + 2*N + 1 + N*j;
  
  /* Compute inverse covariance */
  param[2*N*N + 2*N + 1] = detinv(temp,N);
 
  free(temp); 
}


/* Function to add (or substract)  a point to a level set */

void add_point(double *I, double *param, int N, int sign){

  int l,j;
  
  param[0] += sign;
  for(j=0;j<N;j++){
    param[j+1] += sign*I[j];
    for(l=0;l<N;l++){
      param[j*N+l + (N+1)] += sign*I[j]*I[l];
    }
  }
}
