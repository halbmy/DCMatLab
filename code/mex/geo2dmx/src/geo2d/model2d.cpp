/*
 * model2d.cpp 
 */

#include "model2d.h"

double pythag(double a, double b);

model2d::model2d()
/* Standard constructor */
{
   flag1d = false;
   chgrd = false; chsrc = false;
   nxgrd = 0; nzgrd = 0;      // Initialize size values of vectors and arrays
   nxcnd = 0; nzcnd = 0;
   nsrc = 0;  npot = 0;
   ndata = 0;
   xgrd = NULL; zgrd = NULL;
   xsrc = NULL; zsrc = NULL;
   xpot = NULL; zpot = NULL;
   cond = NULL; confac = NULL;
   is1 = NULL; is2 = NULL;
   ip1 = NULL; ip2 = NULL;
}

model2d::~model2d()
/* Destructor */
{
    if (xgrd != NULL) FREE(xgrd);   // Free allocated memory of vectors and arrays
    if (zgrd != NULL) FREE(zgrd);
    if (xsrc != NULL) FREE(xsrc);
    if (zsrc != NULL) FREE(zsrc);
    if (xpot != NULL) FREE(xpot);
    if (zpot != NULL) FREE(zpot);
    if (cond != NULL) FREE(cond);
    if (confac != NULL) FREE(confac);
    if (is1 != NULL) FREE(is1);
    if (is2 != NULL) FREE(is2);
    if (ip1 != NULL) FREE(ip1);
    if (ip2 != NULL) FREE(ip2);
}

void model2d::setmodel(int nx, double *x, int nz, double *z, double *c)
/* Set up grid and conductivities */
{
   int i, j, k, n;
   
   if (nxgrd != nx)
   {
      chgrd = true;
      nxgrd = nx;
      if (xgrd != NULL) FREE(xgrd);
      xgrd = (double *) CALLOC(nxgrd, sizeof(double));
   }
   if (nzgrd != nz)
   {
      chgrd = true;
      nzgrd = nz;
      if (zgrd != NULL) FREE(zgrd);
      zgrd = (double *) CALLOC(nzgrd, sizeof(double));
   }
   
   nxcnd = nxgrd - 1;
   nzcnd = nzgrd - 1;
   n = nxcnd*nzcnd;
   if (cond == NULL)
   {
      cond = (double *) CALLOC(nxcnd*nzcnd, sizeof(double));
   }
   else
   {
      if (chgrd)
      {
         FREE(cond);
         cond = (double *) CALLOC(nxcnd*nzcnd, sizeof(double));
      }
   }
   
   for (i = 0; i < nxgrd; i++)
   {
      xgrd[i] = x[i];
   }
   for (i = 0; i < nzgrd; i++)
   {
      zgrd[i] = z[i];
   }
   for (i = 0; i < n; i++)
   {
      cond[i] = c[i];
   }

   /* Check if model is 1D, use faster modelling routine then */
   flag1d = true;
   for (i = 1; i < nxcnd; i++)
   {
      for (j = 0, k = i*nzcnd; j < nzcnd; j++, k++)
      {
         if (cond[k-nzcnd] != cond[k])
         {
            flag1d = false;
            i = nxcnd; // exit outer loop as well
            break;
         }
      }
   }
}

void model2d::setsources(int ns, double *xs, double *zs,
                         int np, double *xp, double *zp)
/* Sets ns sources at locations given by the vectors xs and zs
    and  np receiver locations given by vectors xp and zp */
{
   double x, z, xmin, xmax, zmin, zmax;
   int i;
 
   if (nsrc != ns)
   {
      chsrc = true;
      nsrc = ns;
      if (xsrc != NULL) FREE(xsrc);
      if (zsrc != NULL) FREE(zsrc);
      xsrc = (double *) CALLOC(ns, sizeof(double));
      zsrc = (double *) CALLOC(ns, sizeof(double));
   }
   if (npot != np)
   {
      chsrc = true;
      npot = np;
      if (xpot != NULL) FREE(xpot);
      if (zpot != NULL) FREE(zpot);
      xpot = (double *) CALLOC(np, sizeof(double));
      zpot = (double *) CALLOC(np, sizeof(double));
   }   
   
   if (xgrd != NULL && zgrd != NULL)
   {
      xmin = min(xgrd[0],xgrd[nxgrd-1]);
      xmax = max(xgrd[0],xgrd[nxgrd-1]);
      zmin = min(zgrd[0],zgrd[nzgrd-1]);
      zmax = max(zgrd[0],zgrd[nzgrd-1]);
      for (i=0; i<ns; i++)
      {
         x = xs[i];
         z = zs[i];
         if (x >= xmin && x <= xmax && z >= zmin && z <= zmax)
         {
            xsrc[i] = x;
            zsrc[i] = z;
          }
         else
         {
            PRINTF("%d-th source electrode (%.1f,%.1f) is placed out of"
               " the model grid.\n", i+1, x, z);
            exit(0);
         }
      }
      for (i=0; i<np; i++)
      {
         x = xp[i];
         z = zp[i];
         if (x >= xmin && x <= xmax && z >= zmin && z <= zmax)
         {
            xpot[i] = x;
            zpot[i] = z;
         }
         else
         {
            PRINTF("%d-th potential electrode (%.1f,%.1f) is placed out of"
               " the model grid.\n", i+1, x, z);
            exit(0);
         }
      }
      
   }
   else
   {   
      for (i=0; i<ns; i++)
      {
         xsrc[i] = xs[i];
         zsrc[i] = zs[i];
      }
      for (i=0; i<np; i++)
      {
         xpot[i] = xp[i];
         zpot[i] = zp[i];
      }
   }
}

void model2d::setsources(int ns, double *xs, double *zs,
      int nd, int **id)
/* Sets ns source/receiver locations given by the vectors xs and zs
   and nd source/receiver combinations given in 2D-array id. Calculate
   configuration factors for each source/receiver combination */
{
   int i, is, ip;
   double t11, t12, t21, t22;

   setsources(ns, xs, zs, ns, xs, zs);

   if (ndata != nd)
   {
      ndata = nd;
      if (is1 != NULL) FREE(is1);
      if (is2 != NULL) FREE(is2);
      if (ip1 != NULL) FREE(ip1);
      if (ip2 != NULL) FREE(ip2);
      if (confac != NULL) FREE(confac);
      is1 = (int *) CALLOC(ndata, sizeof(int));
      is2 = (int *) CALLOC(ndata, sizeof(int));
      ip1 = (int *) CALLOC(ndata, sizeof(int));
      ip2 = (int *) CALLOC(ndata, sizeof(int));
      confac = (double *) CALLOC(ndata, sizeof(double));
   }

   for (i = 0; i < ndata; i++)
   {
      is1[i] = id[0][i];
      is2[i] = id[1][i];
      ip1[i] = id[2][i];
      ip2[i] = id[3][i];
   }

   for (i = 0; i < ndata; i++)
   {
      is = is1[i];
      if (is--)
      {
         ip = ip1[i];
         if (ip--)
         {
            t11 = 1/pythag(xsrc[is] - xpot[ip], zsrc[is] - zpot[ip]);
         }
         else
         {
            t11 = 0.0;
         }
         ip = ip2[i];
         if (ip--)
         {
            t12 = 1/pythag(xsrc[is] - xpot[ip], zsrc[is] - zpot[ip]);
         }
         else
         {
            t12 = 0.0;
         }
      }
      else
      {
         t11 = 0.0;
         t12 = 0.0;
      }
      is = is2[i];
      if (is--)
      {
         ip = ip1[i];
         if (ip--)
         {
            t21 = 1/pythag(xsrc[is] - xpot[ip], zsrc[is] - zpot[ip]);
         }
         else
         {
            t21 = 0.0;
         }
         ip = ip2[i];
         if (ip--)
         {
            t22 = 1/pythag(xsrc[is] - xpot[ip], zsrc[is] - zpot[ip]);
         }
         else
         {
            t22 = 0.0;
         }
      }
      else
      {
         t21 = 0.0;
         t22 = 0.0;
      }
      confac[i] = M_PI2 / (t11 - t12 - t21 + t22);
   }
}

