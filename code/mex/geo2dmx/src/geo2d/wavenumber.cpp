/*
 * wavenumber.cpp
 */

#include "geo2d.h"

void gauleg(double x1, double x2, double *x, double *w, int n);
void gaulag(double *x, double *w, int n, double alf);
double gammln(double xx);

void geo2d::wavenumber()
{
	double k0;
	int i;

   i = config.NoLegendre + config.NoLaguerre;
   if (i != nw)
   {
      nw = i;
      if (wavnum != NULL) FREE(wavnum);
      if (weight != NULL) FREE(weight);
      wavnum = (double *) CALLOC(nw, sizeof(double));
      weight = (double *) CALLOC(nw, sizeof(double));
   }
	
  /* Minimal wavenumber depends on grid size */
	k0 = xcell[0];
	for (i=1; i<cx; i++){
		k0 = min(k0,xcell[i]);
	}
	for (i=0; i<cz; i++){
		k0 = min(k0,zcell[i]);
	}
   k0 = 0.5 / k0;

  /* Gauss-Legendre abscissas and weights for (0..1) interval,
  	 legendre points */
  gauleg(0.0, 1.0, wavnum, weight, config.NoLegendre);

  /* Scaling according to real wavenumber interval */
  for (i=0; i<config.NoLegendre; i++)
  {
    weight[i] = 2.0*k0/M_PI * weight[i]*wavnum[i];	/* weights for summation */
    wavnum[i] = k0*wavnum[i]*wavnum[i]; 				/* proper wavenumbers */
  }

  /* Gauss-Laguerre abscissas and weights for laguerre points */
  gaulag(wavnum+config.NoLegendre, weight+config.NoLegendre, config.NoLaguerre, 0.0);

  /* Rescale wavenumbers and weights */
  for (i=config.NoLegendre; i<nw; i++)
  {
    weight[i] = k0 * exp(wavnum[i])*weight[i]/M_PI;
    wavnum[i] = k0 * (wavnum[i]+1.0);
  }
}
