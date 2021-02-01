/*
 * model2d.h
 *
 * (C) 2003 C. Schwarzbach
 * schwarzb@geophysik.tu-freiberg.de
 * 04-07-2003
 *
 */

#ifndef MODEL2D_H
#define MODEL2D_H

#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define PRINTF printf
#define CALLOC calloc
#define FREE free

#ifndef M_PI
   #define M_PI 3.1415926535897931
#endif
#ifndef M_PI2
   #define M_PI2 6.283185307179586
#endif
#ifndef min
  #define min(a,b) (a < b) ? a : b
#endif
#ifndef max
  #define max(a,b) (a > b) ? a : b
#endif 

/*
 * model2d
 * Basis class for 2D resistivity modelling.
 * Provides space and routines to set up the model geometry,
 * source and potential electrode locations.
 */

class model2d 
{
private:
protected:
   bool flag1d;   /* signal to use faster 1D modelling routine */
   int nxgrd, nzgrd, nxcnd, nzcnd, nsrc, npot, ndata;
   double *xgrd, *zgrd, *cond, *xsrc, *zsrc, *xpot, *zpot, *confac;
   int *is1, *is2, *ip1, *ip2;
   bool chgrd, chsrc;
public:
   model2d();  /* Standard constructor */
   ~model2d(); /* Standard destructor */
   void setmodel(int nx, double *x, int nz, double *z, double *c);
   /* setup of model geometry
    * nx ... no. of grid nodes in horizontal direction
    * x[0,...,nx-1] ... coordinates of horizontal grid nodes
    * nz ... no. of grid nodes in vertical direction
    * z[0,...,nz-1] ... coordinates of vertical grid nodes
    * c[0,...,(nx-1)*(nz-1)-1] ... array containing conductivities of 
    *        grid cells in columnwise order; c[i*(nz-1)+j] is the
    *        conductivity of cell between (x[i],z[j]) and (x[i+1],z[j+1]) 
    */
   void setsources(int ns, double *xs, double *zs,
      int np, double *xp, double *zp);
   /* setup electrode positions
    * ns ... no. of current electrodes
    * xs[0,...,ns-1] ... x-coordinates of current electrodes
    * zs[0,...,ns-1] ... z-coordinates of current electrodes
    * np ... no. of potential electrodes
    * xp[0,...,np-1] ... x-coordinates of potential electrodes
    * zp[0,...,np-1] ... z-coordinates of potential electrodes
    */
   void setsources(int ns, double *xs, double *zs,
      int nd, int **id);
   /* setup electrode positions and measurement layout
    * ns ... no. of electrodes
    * xs[0,...,ns-1] ... x-coordinates of electrodes
    * zs[0,...,ns-1] ... z-coordinates of electrodes
    * nd ... no. of data points
    * id[0][0,...,ns-1] ... indices of 1st current electrodes
    * id[1][0,...,ns-1] ... indices of 2nd current electrodes
    * id[2][0,...,ns-1] ... indices of 1st potential electrodes
    * id[3][0,...,ns-1] ... indices of 2nd potential electrodes
    *        index no. should be in the range of [1,2,...,ns];
    *        index 0 indicates electrode at (x,z)->infinity
    * see schematic datafile below
    */
};

/* schematic datafile:

ns
xs[0]  zs[0]
xs[1]  zs[1]
 :      :
xs[ns] zs[ns]
nd
id[0][0]  id[1][0]  id[2][0]  id[3][0]  dataval[0]  datastd[0]
id[0][1]  id[1][1]  id[2][1]  id[3][1]  dataval[1]  datastd[1]
 :         :         :         :           :           :
id[0][nd] id[1][nd] id[2][nd] id[3][nd] dataval[nd] datastd[nd]

 * (dataval and datastd ... value and standard deviation of measured
 *                          apparent resistivity)
 */
#endif

