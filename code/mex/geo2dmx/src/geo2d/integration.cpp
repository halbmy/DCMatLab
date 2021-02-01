/*
 * integration.cpp
 */

#include "geo2d.h"

void geo2d::integration(int iw)
/* Performs integration for inverse Fourier transform. */
{
	double w = weight[iw];
	int i, n = nsrc*nx*nz;

	if (iw == nw-1){	// initialize phi at first call to integration()
		for (i=0; i<n; i++){
			phi[i] = w*b[i];
		}
	}
	else {			// successive calls have to add
		for (i=0; i<n; i++){
			phi[i] = phi[i] + w*b[i];
		}
	}
}

