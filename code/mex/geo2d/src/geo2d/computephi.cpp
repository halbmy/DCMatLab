/*
 * computephi.cpp 
 */

#include "geo2d.h"

void geo2d::computephi(double *potential)
/*
 * User interface to start forward-modelling.
 *
 * Arguments:
 * 1) potential ... pointer to a matrix of size nsrc x npot;
 *                  the potential of each source is written columnwise
 *                  to this array.
 * 2) Control prolongation of grid to improve validity of boundary
 *    conditions
 *    prolnum ... number of grid nodes/cells to add at each edge
 *    prolfac ... factor to stretch cell size
 *    prolmod ... add cells with conductivity chosen from
 *     'a' or 'A' ... average of model conductivities
 *     'b' or 'B' ... conductivity of old boundary cells
 * 3) Control refinement of the grid
 *    refine ... number of additional grid nodes in horizontal/vertical
 *               direction within each cell
 * 4) Set number of wavenumbers for quadrature
 *    legendre ... first, use gauss-legendre quadrature formulae
 *    laguerre ... second, use gauss-laguerre quadrature formulae
 * 5) bound ... type of boundary condition
 *     'd' or 'D' ... Dirichlet
 *     'm' or 'M' ... mixed a la Dey & Morrison (1979)
 * 6) backgrnd ... choice of background conductivity for singularity removal
 *     'a' or 'A' ... take average conductivity of the non-prolonged grid
 *     'm' or 'M' ... take average of conductivity at each source location
 *     's' or 'S' ... take conductivity of each source location
 *                    - exact singularity removal
 */

{
   if (flag1d)
   {
      potential1d(potential);
   }
   else
   {
	int n, iw, isrc;

  /* refine grid and add cells to the boundary */
  grid();

  /* compute wavenumbers and quadrature weights */
  wavenumber();

  /* allocate memory for matrix of coupling coefficients:
     spd banded matrix type for direct solver (LAPACK) */
   n = nx*nz;
   if (chgrd)
   {
      if(C != NULL) FREE(C);
      C = (double *) CALLOC((nz+1)*n, sizeof(double));
   }
   if (chgrd || chsrc)
   {
      /* matrix of source terms, will be overwritten by solution of C*phi=b */
      if (b != NULL) FREE(b);
      b = (double *) CALLOC(n*nsrc, sizeof(double));
      /* matrix of potentials, accumulates integration results at each wavenumber */
      if (phi != NULL) FREE(phi);
      phi = (double *) CALLOC(n*nsrc, sizeof(double));
      chsrc = false;
   }  
  /* compute normal coupling coefficients */
	coupnesa();

  /* in case of averaged background conductivity compute anomalous
     coupling coefficients */
	if (config.Background == 'a')
  {
    coupaesa(0);
    chgrd = false;
  }

	/* loop over wavenumbers */
	for (iw = nw-1; iw>=0; iw--)
  {
		switch (config.Background)
    {
		case 'a':
			/* compute anomalous self-coupling coefficient ca_p */
			coupap(wavnum[iw]);
			/* compute right hand side b */
			righthandside(wavnum[iw]);
			break;
    case 'm':
    case 's':
			for (isrc=0; isrc<nsrc; isrc++)
      {
				/* compute anomalous coupling coefficients ca_e, ca_s, ca_a */
				coupaesa(isrc);
            chgrd = false;
				/* compute anomalous self-coupling coefficient ca_p */
				coupap(wavnum[iw]);
				/* compute right hand side b */
				righthandside(wavnum[iw],isrc);
			}
			break;
		}
	
		/* assemble spd banded matrix C */
    coupn(wavnum[iw]);
	
		/* solve C*phi=b, b contains solution phi on exit */
		solvelineq();

		/* integration of anomalous potential (inverse Fourier transform) */
		integration(iw);
	}
  /* end of wavenumber loop */

	/* interpolate anomalous potential for positions of potential electrodes
     and add normal potential */
  superpos(potential);
   }
}

