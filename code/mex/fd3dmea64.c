#include <math.h>
#include "mex.h"
#include "matrix.h"
#include "ldl.h"

/* MEA=fd3dmea(x,y,z,El,C,C1,Potmap,p) */

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{

  double *x, *y, *z, *SIG, *EL, *phip, *Ax, *A1x, *Potx, *MEA, *rhoBg;
  mwSize nx, ny, nz, nn, nel, n;
  mwIndex *Ai, *Ap, *A1i, *A1p, *Poti, *Potp, *Lp, *Li;
  int ix, iy, iz, ie, je, ii, jj;
  double dx, dy, dz, ma, sdel;
  const double _PI=3.1415926535897;
  long *P, *Pinv, *Lnz, lnz, *Parent, *Pattern, *Flag;
  double *Y, *D, *Lx, *X, *B, *p;

  if (nrhs != 9) 
    mexErrMsgTxt("Nine inputs required.");
  if ((nlhs != 1)&(nlhs != 2))
    mexErrMsgTxt("One output required.");

  nn=0; 
  nx=mxGetM(prhs[nn])*mxGetN(prhs[nn]);
  x=mxGetPr(prhs[nn]);
  nn++;
  ny=mxGetM(prhs[nn])*mxGetN(prhs[nn]);
  y=mxGetPr(prhs[nn]);
  nn++;
  nz=mxGetM(prhs[nn])*mxGetN(prhs[nn]);
  z=mxGetPr(prhs[nn]);
  nn++;
  nel=mxGetM(prhs[nn]);
  EL=mxGetPr(prhs[nn]);
  nn++;
  n = mxGetM(prhs[nn]);
  Ai = mxGetIr(prhs[nn]);
  Ap = mxGetJc(prhs[nn]);
  Ax = mxGetPr(prhs[nn]);
  nn++;
  A1i = mxGetIr(prhs[nn]);
  A1p = mxGetJc(prhs[nn]);
  A1x = mxGetPr(prhs[nn]);
  nn++;
  Poti = mxGetIr(prhs[nn]);
  Potp = mxGetJc(prhs[nn]);
  Potx = mxGetPr(prhs[nn]);
  nn++;
  rhoBg = mxGetPr(prhs[nn]);
  nn++;
  p = mxGetPr(prhs[nn]);

  P = (long *) mxMalloc( n * sizeof( long ) );
  Pinv = (long *) mxMalloc( n * sizeof( long ) );
  for( ii=0 ; ii<n ; ii++) { P[ ii ] = p[ii]-1; }
  phip = (double *) mxMalloc( n * sizeof(double) );

  Lp = (long *) mxMalloc( (n+1)*sizeof(long) );
  nn = (n == 0) ? 1 : n;
  Parent = (long *) mxMalloc( nn*sizeof(long) );
  Y = (double *) mxMalloc ( nn*sizeof(double) );
  Flag = (long *) mxMalloc( nn*sizeof(int) );  
  Pattern = (long *) mxMalloc( nn*sizeof(long) );  
  Lnz = (long *) mxMalloc( nn*sizeof(long) );  
  LDL_symbolic (n, Ap, Ai, Lp, Parent, Lnz, Flag, P, Pinv);
  D = (double *) mxMalloc( nn*sizeof(double) );
   
  LDL_symbolic(n, Ap, Ai, Lp, Parent, Lnz, Flag, P, Pinv);
  lnz = Lp[n];
  Li = (long *) mxMalloc ( (lnz+1)*sizeof(long) );
  Lx = (double *) mxMalloc ( (lnz+1)*sizeof(double) );
  LDL_numeric(n, Ap, Ai, Ax, Lp, Parent, Lnz, Li, Lx, D, Y, Flag, Pattern, P, Pinv);
 
  plhs[0]=mxCreateDoubleMatrix(nel,nel,mxREAL);
  MEA = mxGetPr(plhs[0]);

  X = (double *) mxMalloc ( n*sizeof(double) );
  sdel = sqrt(pow(EL[1]-EL[2],2)+pow(EL[nel+1]-EL[nel+2],2));
  if(sdel==0) sdel=1;
  
  for( ie=0 ; ie<nel ; ie++ ) {
    for( je=0 ; je<nel ; je++) {
     if (ie!=je) MEA[ ie*nel+je ] = rhoBg[ie]*0.5/_PI/sqrt(pow(EL[ie]-EL[je],2)+pow(EL[nel+ie]-EL[nel+je],2)+1e-8);
    }
    ii=0;
    for( iz=0 ; iz<nz ; iz++ ) {
      dz = z[iz];
      for( iy=0 ; iy<ny ; iy++ ) {
        dy = y[iy] - EL[nel+ie];
        for( ix=0 ; ix<nx ; ix++ ) {
          dx = x[ix] - EL[ie];
          ma = dx*dx + dy*dy + dz*dz;
          if( ma==0 ) ma = sdel / 5;
          phip[ ii++ ] = 1.0/( 2 * _PI * sqrt( ma ) );
        }
      }
    }
    jj=0;
    for( ii=0 ; ii<n ; ii++) { X[ii]=0; }
    for( ii=0 ; ii<n ; ii++) {
      while( jj<Ap[ii+1] ) {
	    X[Ai[jj]] += ( A1x[jj] - rhoBg[ie] * Ax[jj] ) * phip[ii];
	    jj++;
      }
    }
    
    LDL_perm (n, Y, X, P);
    LDL_lsolve (n, Y, Lp, Li, Lx);
    LDL_dsolve (n, Y, D);
    LDL_ltsolve(n, Y, Lp, Li, Lx);
    LDL_permt (n, X, Y, P);
    jj=0;
    for( ii=0 ; ii<n ; ii++) {
      while( jj<Potp[ii+1] ) {
/*       	MEA[ ie*nel + Poti[jj] ] += Potx[jj] * ( phip[ii] * rhoBg[ie] + X[ii] ); */
       	MEA[ ie*nel + Poti[jj] ] += Potx[jj] * X[ii];
	    jj++;
      }
    }
    /* ma=0.0;
    for( ii=0 ; ii<nx*ny*nz ; ii++ ) { if( X[ii]<ma ) ma=X[ii]; } */
    MEA[ ie*nel + ie ] = 0;
  }
  if(nlhs>1) {
    plhs[1]=mxCreateDoubleMatrix(n,1,mxREAL);
    phip = mxGetPr(plhs[1]);
    for( ii=0 ; ii<n ; ii++) { phip[ii] = X[ii]; }
  }
  mxFree(Lp); mxFree(Li); mxFree(Lx); mxFree(D); mxFree(P); mxFree(Pinv);
  mxFree(Parent); mxFree(Y); mxFree(Flag); mxFree(Pattern); mxFree(Lnz);  
}
