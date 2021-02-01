/* geo2d.h */

/***************************************************************************
 * GEO2D
 ***************************************************************************
 *
 * Calculate the potential of dc-sources on a 2D conductivity distribution in 
 * the earth.
 * Features:
 * - finite differences modelling code (Dey and Morrison, 1979)
 * - singularity removal technique (Lowry et al., 1989; Zhao and Yedlin, 1996)
 *   extended for sources placed on conductivity discontinuities
 * - direct linear equations solver (LAPACK)
 * - inverse Fourier transform by Gauss-quadrature
 *
 ***************************************************************************
 *
 * (C) 2001-2003 R.-U. Börner and C. Schwarzbach
 * rub@geophysik.tu-freiberg.de
 * schwarzb@geophysik.tu-freiberg.de
 * 04-07-2003
 *
 **************************************************************************/

#ifndef GEO2D_H
#define GEO2D_H

#include "model2d.h"

struct Geo2dConfig
{
   double FacProlX;
   double FacProlZ;
   int NoEquiX, NoProlX, NoRefineX;
   int NoEquiZ, NoProlZ, NoRefineZ;
   int NoLegendre;
   int NoLaguerre;
   char ModProlong;
   char Boundary;
   char Background;
};

/***************************************************************************
*                                                                          *
* Class definition of geo2d derived from class model2d. The latter         *
* provides basic geometry and conductivity data and methods.               *
* geo2d provides the appriate methods to perform a 2d forward modelling of *
* the potential of a point current source.                                 *
*                                                                          *
***************************************************************************/

class geo2d : public model2d
{
private:
   /* variable declarations */
   Geo2dConfig config;                     /* configuration parameters */
   int nx, nz, cx, cz;                     /* geometry, conductivity */
   int nw;                                 /* wavenumbers */
   double *xnode, *znode, *xcell, *zcell;  /* grid geometry */
   double *cnorm, *cdiff;                  /* conductivities */
   double *cne, *cns, *cna;                /* normal coupling coefficients */
   double *cae, *cas, *caa, *cap;          /* anomalous coupling coefficients */
   double *wavnum, *weight;                /* wavenumbers and quadrature weights */
   double *C, *b, *phi;                    /* arrays for linear equations */
   double **cback;                         /* arrays for source conductivities */
   int **iback;
   
   /* function declarations */
   void readconfiguration(char *configname);
   void grid();
   char sourcepos(int isrc, int *ix, int *iz);
   int locate(double *xx, int n, double x);
   void wavenumber();
   void coupnesa();
   void coupaesa(int isrc);
   void coupap(double lambda);
   void coupn(double lambda);
   void righthandside(double lambda);
   void righthandside(double lambda, int isrc);
   void solvelineq();
   void integration(int iw);
   void superpos(double *pot);
   double interpol(int is, int ip);
   double potential(int isrc, double h, double x, double z);
   double potential(int isrc, double h, double x, double z, double lambda);
   void potential1d(double *potential){ cerr << "function unimplemented" << endl; }
   double dpotdx(int isrc, double h, double x, double z, double lambda);
   double dpotdz(int isrc, double h, double x, double z, double lambda);
   
public:
   geo2d();  /* standard constructor */
   ~geo2d(); /* standard destructor */
   void computephi(double *potential);
   /* Compute potential at each potential electrode for each source electrode.
    * Arguments:
    * 1) Array of potentials
    * potential[0,...,nsrc*npot-1] ... contains potential values on return
    *          potential[i*npot+j] contains potential for (i+1)-th source
    *          electrode at (j+1)-th potential electrode
    *          (i = 0,...,nsrc-1; j = 0,...,npot-1)
    * 2) Control prolongation of grid to improve validity of boundary
    *    conditions
    * prolnum ... number of grid nodes/cells to add at each edge
    * prolfac ... factor to stretch cell size
    * prolmod ... add cells with conductivity chosen from
    *     'a' or 'A' ... average of model conductivities
    *     'b' or 'B' ... conductivity of old boundary cells
    * 3) Control refinement of the grid
    * refine ... number of additional grid nodes in horizontal/vertical
    *               direction within each cell
    * 4) Set number of wavenumbers for quadrature
    * legendre ... no. of gauss-legendre quadrature abscissas
    * laguerre ... no. of gauss-laguerre quadrature abscissas
    * 5)
    * bound ... type of boundary condition
    *     'd' or 'D' ... Dirichlet
    *     'm' or 'M' ... mixed a la Dey & Morrison (1979)
    * 6)
    * backgrnd ... choice of background conductivity for singularity removal
    *     'a' or 'A' ... take average conductivity of the non-prolonged grid
    *     'm' or 'M' ... take average of conductivity at each source location
    *     's' or 'S' ... take conductivity of each source location
    *                    - exact singularity removal
    */
   double computerhoa(double *rhoa);
   /* Compute apparent resistivities 
    * rhoa[0,...,nd-1] ... apparent resistivity values for data points i = 1,...,nd.
    * For the other arguments see 'computephi' above.
    * computephi returns maximum reciprocity error of voltage/rhoa from all data points.
    */
};

#endif

