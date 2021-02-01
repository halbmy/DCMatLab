/*
 * superpos.h
 */
#include "geo2d.h"

void geo2d::superpos(double *pot)
/*
 * Performs superposition of anomalous potential calculated
 * before and normal potential calculated here.
 */
{
   int i, j, k;
   double x, z, h;
   
   k = 0;
   for (i=0; i<nsrc; i++)
   {
      h = zsrc[i] - znode[0];
      for (j=0; j<npot; j++)
      {
         x = xpot[j] - xsrc[i];
         z = zpot[j] - znode[0];
         pot[k++] = interpol(i,j) + potential(i,h,x,z);
      }
   }
}

double geo2d::interpol(int is, int ip)
/*
 * Perform interpolation to get anomalous potential for (arbitrarily
 * placed) potential electrode ip with respect to source electrode is.
 */
{
   double p, s, t;
   int i, j, k;
   
   /* First, locate cell within which potential is to be calculated */
   i = locate(xnode, nx, xpot[ip]);
   i = max(i,0);
   i = min(i,nx-2);
   j = locate(znode, nz, zpot[ip]);
   j = max(j,0);
   j = min(j,nz-2);
   k = (is*nx + i) * nz + j;
   
   /* Second, calculate bilinear interpolation */
   s = (xpot[ip] - xnode[i]) / (xnode[i+1] - xnode[i]);
   t = (zpot[ip] - znode[j]) / (znode[j+1] - znode[j]);
   p = (1.0 - s) * (1.0 - t) * phi[k] + 
      (1.0 - s) *        t  * phi[k+1] +
      s  * (1.0 - t) * phi[k+nz] +
      s  *        t  * phi[k+nz+1];
   
   return(p);
}

