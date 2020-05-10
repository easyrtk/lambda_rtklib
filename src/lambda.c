/*------------------------------------------------------------------------------
* lambda.c : integer ambiguity resolution
*
*          Copyright (C) 2007-2008 by T.TAKASU, All rights reserved.
*
* reference :
*     [1] P.J.G.Teunissen, The least-square ambiguity decorrelation adjustment:
*         a method for fast GPS ambiguity estimation, J.Geodesy, Vol.70, 65-82,
*         1995
*     [2] X.-W.Chang, X.Yang, T.Zhou, MLAMBDA: A modified LAMBDA method for
*         integer least-squares estimation, J.Geodesy, Vol.79, 552-565, 2005
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
* history : 2007/01/13 1.0 new
*           2015/05/31 1.1 add api lambda_reduction(), lambda_search()
*-----------------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>

#include "lambda.h"

/* constants/macros ----------------------------------------------------------*/

#define LOOPMAX     10000           /* maximum count of search loop */

#define SGN(x)      ((x)<=0.0?-1.0:1.0)
#define ROUND(x)    (floor((x)+0.5))
#define SWAP(x,y)   do {double tmp_; tmp_=x; x=y; y=tmp_;} while (0)

/* fatal error ---------------------------------------------------------------*/
static void fatalerr(const char *format, ...)
{
    va_list ap;
    va_start(ap,format); vfprintf(stderr,format,ap); va_end(ap);
    exit(-9);
}

/* new matrix ------------------------------------------------------------------
* allocate memory of matrix 
* args   : int    n,m       I   number of rows and columns of matrix
* return : matrix pointer (if n<=0 or m<=0, return NULL)
*-----------------------------------------------------------------------------*/
static double *mat(int n, int m)
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
static int *imat(int n, int m)
{
    int *p;
    
    if (n<=0||m<=0) return NULL;
    if (!(p=(int *)malloc(sizeof(int)*n*m))) {
        fatalerr("integer matrix memory allocation error: n=%d,m=%d\n",n,m);
    }
    return p;
}
/* zero matrix -----------------------------------------------------------------
* generate new zero matrix
* args   : int    n,m       I   number of rows and columns of matrix
* return : matrix pointer (if n<=0 or m<=0, return NULL)
*-----------------------------------------------------------------------------*/
static double *zeros(int n, int m)
{
    double *p;
    
#if NOCALLOC
    if ((p=mat(n,m))) for (n=n*m-1;n>=0;n--) p[n]=0.0;
#else
    if (n<=0||m<=0) return NULL;
    if (!(p=(double *)calloc(sizeof(double),n*m))) {
        fatalerr("matrix memory allocation error: n=%d,m=%d\n",n,m);
    }
#endif
    return p;
}
/* identity matrix -------------------------------------------------------------
* generate new identity matrix
* args   : int    n         I   number of rows and columns of matrix
* return : matrix pointer (if n<=0, return NULL)
*-----------------------------------------------------------------------------*/
static double *eye(int n)
{
    double *p;
    int i;
    
    if ((p=zeros(n,n))) for (i=0;i<n;i++) p[i+i*n]=1.0;
    return p;
}

/* copy matrix -----------------------------------------------------------------
* copy matrix
* args   : double *A        O   destination matrix A (n x m)
*          double *B        I   source matrix B (n x m)
*          int    n,m       I   number of rows and columns of matrix
* return : none
*-----------------------------------------------------------------------------*/
static void matcpy(double *A, const double *B, int n, int m)
{
    memcpy(A,B,sizeof(double)*n*m);
}

/* multiply matrix -----------------------------------------------------------*/
static void matmul(const char *tr, int n, int k, int m, double alpha,
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
}
/* LU decomposition ----------------------------------------------------------*/
static int ludcmp(double *A, int n, int *indx, double *d)
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
static void lubksb(const double *A, int n, const int *indx, double *b)
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
static int matinv(double *A, int n)
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
/* solve linear equation -----------------------------------------------------*/
static int solve(const char *tr, const double *A, const double *Y, int n,
                 int m, double *X)
{
    double *B=mat(n,n);
    int info;
    
    matcpy(B,A,n,n);
    if (!(info=matinv(B,n))) matmul(tr[0]=='N'?"NN":"TN",n,m,n,1.0,B,Y,0.0,X);
    free(B);
    return info;
}

/* LD factorization (Q=L'*diag(D)*L) -----------------------------------------*/
static int LD(int n, const double *Q, double *L, double *D)
{
    int i,j,k,info=0;
    double a,*A=mat(n,n);
    
    memcpy(A,Q,sizeof(double)*n*n);
    for (i=n-1;i>=0;i--) {
        if ((D[i]=A[i+i*n])<=0.0) {info=-1; break;}
        a=sqrt(D[i]);
        for (j=0;j<=i;j++) L[i+j*n]=A[i+j*n]/a;
        for (j=0;j<=i-1;j++) for (k=0;k<=j;k++) A[j+k*n]-=L[i+k*n]*L[i+j*n];
        for (j=0;j<=i;j++) L[i+j*n]/=L[i+i*n];
    }
    free(A);
    if (info) fprintf(stderr,"%s : LD factorization error\n",__FILE__);
    return info;
}
/* integer gauss transformation ----------------------------------------------*/
static void gauss(int n, double *L, double *Z, int i, int j)
{
    int k,mu;
    
    if ((mu=(int)ROUND(L[i+j*n]))!=0) {
        for (k=i;k<n;k++) L[k+n*j]-=(double)mu*L[k+i*n];
        for (k=0;k<n;k++) Z[k+n*j]-=(double)mu*Z[k+i*n];
    }
}
/* permutations --------------------------------------------------------------*/
static void perm(int n, double *L, double *D, int j, double del, double *Z)
{
    int k;
    double eta,lam,a0,a1;
    
    eta=D[j]/del;
    lam=D[j+1]*L[j+1+j*n]/del;
    D[j]=eta*D[j+1]; D[j+1]=del;
    for (k=0;k<=j-1;k++) {
        a0=L[j+k*n]; a1=L[j+1+k*n];
        L[j+k*n]=-L[j+1+j*n]*a0+a1;
        L[j+1+k*n]=eta*a0+lam*a1;
    }
    L[j+1+j*n]=lam;
    for (k=j+2;k<n;k++) SWAP(L[k+j*n],L[k+(j+1)*n]);
    for (k=0;k<n;k++) SWAP(Z[k+j*n],Z[k+(j+1)*n]);
}
/* lambda reduction (z=Z'*a, Qz=Z'*Q*Z=L'*diag(D)*L) (ref.[1]) ---------------*/
static void reduction(int n, double *L, double *D, double *Z)
{
    int i,j,k;
    double del;
    
    j=n-2; k=n-2;
    while (j>=0) {
        if (j<=k) for (i=j+1;i<n;i++) gauss(n,L,Z,i,j);
        del=D[j]+L[j+1+j*n]*L[j+1+j*n]*D[j+1];
        if (del+1E-6<D[j+1]) { /* compared considering numerical error */
            perm(n,L,D,j,del,Z);
            k=j; j=n-2;
        }
        else j--;
    }
}
/* modified lambda (mlambda) search (ref. [2]) -------------------------------
* args   : n      I  number of float parameters
*          m      I  number of fixed solution
           L,D    I  transformed covariance matrix
           zs     I  transformed double-diff phase biases
           zn     O  fixed solutions
           s      O  sum of residuals for fixed solutions                    */
static int search(int n, int m, const double *L, const double *D, const double *zs, double *zn, double *s)
{
    int i,j,k,c,nn=0,imax=0;
    double newdist,maxdist=1E99,y;
    double *S=zeros(n,n),*dist=mat(n,1),*zb=mat(n,1),*z=mat(n,1),*step=mat(n,1);
    
    k=n-1; dist[k]=0.0;
    zb[k]=zs[k];
    z[k]=ROUND(zb[k]);
    y=zb[k]-z[k];
    step[k]=SGN(y);  /* step towards closest integer */
    for (c=0;c<LOOPMAX;c++) {
        newdist=dist[k]+y*y/D[k];  /* newdist=sum(((z(j)-zb(j))^2/d(j))) */
        if (newdist<maxdist) {
            /* Case 1: move down */
            if (k!=0) {
                dist[--k]=newdist;
                for (i=0;i<=k;i++)
                    S[k+i*n]=S[k+1+i*n]+(z[k+1]-zb[k+1])*L[k+1+i*n];
                zb[k]=zs[k]+S[k+k*n];
                z[k]=ROUND(zb[k]); /* next valid integer */
                y=zb[k]-z[k];
                step[k]=SGN(y);
            }
            /* Case 2: store the found candidate and try next valid integer */
            else {
                if (nn<m) {  /* store the first m initial points */
                    if (nn==0||newdist>s[imax]) imax=nn;
                    for (i=0;i<n;i++) zn[i+nn*n]=z[i];
                    s[nn++]=newdist;
                }
                else {
                    if (newdist<s[imax]) {
                        for (i=0;i<n;i++) zn[i+imax*n]=z[i];
                        s[imax]=newdist;
                        for (i=imax=0;i<m;i++) if (s[imax]<s[i]) imax=i;
                    }
                    maxdist=s[imax];
                }
                z[0]+=step[0]; /* next valid integer */
                y=zb[0]-z[0];
                step[0]=-step[0]-SGN(step[0]);
            }
        }
        /* Case 3: exit or move up */
        else {
            if (k==n-1) break;
            else {
                k++;  /* move up */
                z[k]+=step[k];  /* next valid integer */
                y=zb[k]-z[k];
                step[k]=-step[k]-SGN(step[k]);
            }
        }
    }
    for (i=0;i<m-1;i++) { /* sort by s */
        for (j=i+1;j<m;j++) {
            if (s[i]<s[j]) continue;
            SWAP(s[i],s[j]);
            for (k=0;k<n;k++) SWAP(zn[k+i*n],zn[k+j*n]);
        }
    }
    free(S); free(dist); free(zb); free(z); free(step);
    
    if (c>=LOOPMAX) {
        fprintf(stderr,"%s : search loop count overflow\n",__FILE__);
        return -2;
    }
    return 0;
}
/* lambda/mlambda integer least-square estimation ------------------------------
* integer least-square estimation. reduction is performed by lambda (ref.[1]),
* and search by mlambda (ref.[2]).
* args   : int    n      I  number of float parameters
*          int    m      I  number of fixed solutions
*          double *a     I  float parameters (n x 1) (double-diff phase biases)
*          double *Q     I  covariance matrix of float parameters (n x n)
*          double *F     O  fixed solutions (n x m)
*          double *s     O  sum of squared residulas of fixed solutions (1 x m)
* return : status (0:ok,other:error)
* notes  : matrix stored by column-major order (fortran convension)
*-----------------------------------------------------------------------------*/
extern int lambda(int n, int m, const double *a, const double *Q, double *F, double *s)
{
    int info;
    double *L,*D,*Z,*z,*E;
    
    if (n<=0||m<=0) return -1;
    L=zeros(n,n); D=mat(n,1); Z=eye(n); z=mat(n,1); E=mat(n,m);
    
    /* LD (lower diaganol) factorization (Q=L'*diag(D)*L) */
    if (!(info=LD(n,Q,L,D))) {
        
        /* lambda reduction (z=Z'*a, Qz=Z'*Q*Z=L'*diag(D)*L) */
        reduction(n,L,D,Z);
        matmul("TN",n,1,n,1.0,Z,a,0.0,z); /* z=Z'*a */
        
        /* mlambda search 
            z = transformed double-diff phase biases
            L,D = transformed covariance matrix */
        if (!(info=search(n,m,L,D,z,E,s))) {  /* returns 0 if no error */
            
            info=solve("T",Z,E,n,m,F); /* F=Z'\E */
        }
    }
    free(L); free(D); free(Z); free(z); free(E);
    return info;
}
/* lambda reduction ------------------------------------------------------------
* reduction by lambda (ref [1]) for integer least square
* args   : int    n      I  number of float parameters
*          double *Q     I  covariance matrix of float parameters (n x n)
*          double *Z     O  lambda reduction matrix (n x n)
* return : status (0:ok,other:error)
*-----------------------------------------------------------------------------*/
static int lambda_reduction(int n, const double *Q, double *Z)
{
    double *L,*D;
    int i,j,info;
    
    if (n<=0) return -1;
    
    L=zeros(n,n); D=mat(n,1);
    
    for (i=0;i<n;i++) for (j=0;j<n;j++) {
        Z[i+j*n]=i==j?1.0:0.0;
    }
    /* LD factorization */
    if ((info=LD(n,Q,L,D))) {
        free(L); free(D);
        return info;
    }
    /* lambda reduction */
    reduction(n,L,D,Z);
     
    free(L); free(D);
    return 0;
}
/* mlambda search --------------------------------------------------------------
* search by  mlambda (ref [2]) for integer least square
* args   : int    n      I  number of float parameters
*          int    m      I  number of fixed solutions
*          double *a     I  float parameters (n x 1)
*          double *Q     I  covariance matrix of float parameters (n x n)
*          double *F     O  fixed solutions (n x m)
*          double *s     O  sum of squared residulas of fixed solutions (1 x m)
* return : status (0:ok,other:error)
*-----------------------------------------------------------------------------*/
extern int lambda_search(int n, int m, const double *a, const double *Q, double *F, double *s)
{
    double *L,*D;
    int info;
    
    if (n<=0||m<=0) return -1;
    
    L=zeros(n,n); D=mat(n,1);
    
    /* LD factorization */
    if ((info=LD(n,Q,L,D))) {
        free(L); free(D);
        return info;
    }
    /* mlambda search */
    info=search(n,m,L,D,a,F,s);
    
    free(L); free(D);
    return info;
}
