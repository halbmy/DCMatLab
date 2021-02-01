/*
 * coupling.cpp
 */

#include "geo2d.h"
#include "bessel.h"

void geo2d::coupnesa()
/*
 * allocates memory for vectors cne, cns, cna and
 * calculates coupling coefficients cne, cns, and cna (area
 * integrational part) for normal conductivity
 */
{
	int i, j, k, l, n;
	
	n = nx*nz;
   if (chgrd)
   {
      if (cne != NULL) FREE(cne);
      if (cns != NULL) FREE(cns);
      if (cna != NULL) FREE(cna);
	 	cne = (double *) CALLOC(n-nz, sizeof(double));
		cns = (double *) CALLOC(n-1, sizeof(double));
		cna = (double *) CALLOC(n, sizeof(double));
   }
   
 	// Upper left corner
	cne[0] = -(zcell[0]*cnorm[0])/(2.0*xcell[0]);
	cns[0] = -(xcell[0]*cnorm[0])/(2.0*zcell[0]);
	cna[0] = (xcell[0]*zcell[0]*cnorm[0])/4.0;
 	// Left edge
 	for (k=1; k<cz; k++){
 		l = k-1;
 		cne[k] = -(zcell[l]*cnorm[l] + zcell[k]*cnorm[k])
 			/(2.0*xcell[0]);
 		cns[k] = -(xcell[0]*cnorm[k])/(2.0*zcell[k]);
 		cna[k] = (xcell[0]*zcell[l]*cnorm[l] +
 			xcell[0]*zcell[k]*cnorm[k])/4.0;
 	}
 	// Lower left corner
	l = k-1;
	cne[k] = -(zcell[l]*cnorm[l])/(2.0*xcell[0]);
	cns[k] = 0.0;
	cna[k] = (xcell[0]*zcell[l]*cnorm[l])/4.0;
 	for (i=1; i<cx; i++){
 		j = i-1;
 		// Upper edge
		cne[i*nz] = -(zcell[0]*cnorm[i*cz])/(2.0*xcell[i]);
		cns[i*nz] = -(xcell[j]*cnorm[j*cz] + xcell[i]*cnorm[i*cz])
			/(2.0*zcell[0]);
		cna[i*nz] = (xcell[i]*zcell[0]*cnorm[i*cz] +
			xcell[j]*zcell[0]*cnorm[j*cz])/4.0;
		// Inner coefficients
 		for (k=1; k<cz; k++){
 			l = k-1;
 			cne[i*nz+k] = -(zcell[l]*cnorm[i*cz+l] + zcell[k]*cnorm[i*cz+k])
 				/(2.0*xcell[i]);
 			cns[i*nz+k] = -(xcell[j]*cnorm[j*cz+k] + xcell[i]*cnorm[i*cz+k])
 				/(2.0*zcell[k]);
 			cna[i*nz+k] = (xcell[j]*zcell[l]*cnorm[j*cz+l] +
 				xcell[i]*zcell[l]*cnorm[i*cz+l] +
 				xcell[i]*zcell[k]*cnorm[i*cz+k] +
 				xcell[j]*zcell[k]*cnorm[j*cz+k])/4.0;
 		}
 		// Lower edge
 		l = k-1;
 		cne[i*nz+k] = -(zcell[l]*cnorm[i*cz+l])/(2.0*xcell[i]);
 		cns[i*nz+k] = 0.0;
 		cna[i*nz+k] = (xcell[j]*zcell[l]*cnorm[j*cz+l] +
 			xcell[i]*zcell[l]*cnorm[i*cz+l])/4.0;
 	}
 	// Upper right corner
 	j = i-1;
	cns[i*nz] = -(xcell[j]*cnorm[j*cz])/(2.0*zcell[0]);
	cna[i*nz] = (xcell[j]*zcell[0]*cnorm[j*cz])/4.0;
 	// Right edge
	for (k=1; k<cz; k++){
		l = k-1;
		cns[i*nz+k] = -(xcell[j]*cnorm[j*cz+k])/(2.0*zcell[k]);
		cna[i*nz+k] = (xcell[j]*zcell[l]*cnorm[j*cz+l] +
			xcell[j]*zcell[k]*cnorm[j*cz+k])/4.0;
	}
 	// Lower right corner
	l = k-1;
	cna[i*nz+k] = (xcell[j]*zcell[l]*cnorm[j*cz+l])/4.0;
}

void geo2d::coupaesa(int isrc)
/*
 * allocates memory for vectors cae, cas, caa, cap and
 * calculates coupling coefficients cae, cas, and caa (area
 * integrational part) for differential conductivity
 */
{
   double c;
	int i, j, k, l, n;
	
  // Allocate memory for coupling coefficients
	n = nx*nz;
   if (chgrd)
   {
      if (cae != NULL) FREE(cae);
      if (cas != NULL) FREE(cas);
      if (caa != NULL) FREE(caa);
      if (cap != NULL) FREE(cap);
	 	cae = (double *) CALLOC(n-nz, sizeof(double));
 		cas = (double *) CALLOC(n-1, sizeof(double));
 		caa = (double *) CALLOC(n, sizeof(double));
 		cap = (double *) CALLOC(n, sizeof(double));
   }

  // Compute differential conductivity overwriting cdiff
  k = iback[isrc][0];
  l = iback[isrc][1];
  for (i=0; i<k; i++)
  {
    c = cback[isrc][0];
    for (j=0; j<l; j++)
    {
      cdiff[i*cz+j] = c - cnorm[i*cz+j];
    }
    c = cback[isrc][3];
    for (j=l; j<cz; j++)
    {
      cdiff[i*cz+j] = c - cnorm[i*cz+j];
    }
  }
  for (i=k; i<cx; i++)
  {
    c = cback[isrc][1];
    for (j=0; j<l; j++)
    {
      cdiff[i*cz+j] = c - cnorm[i*cz+j];
    }
    c = cback[isrc][2];
    for (j=l; j<cz; j++)
    {
      cdiff[i*cz+j] = c - cnorm[i*cz+j];
    }
  }
    
  // Upper left corner
  cae[0] = -(zcell[0]*cdiff[0])/(2.0*xcell[0]);
  cas[0] = -(xcell[0]*cdiff[0])/(2.0*zcell[0]);
	caa[0] = (xcell[0]*zcell[0]*cdiff[0])/4.0;
 	// Left edge
 	for (k=1; k<cz; k++){
 		l = k-1;
 		cae[k] = -(zcell[l]*cdiff[l] + zcell[k]*cdiff[k])
 			/(2.0*xcell[0]);
 		cas[k] = -(xcell[0]*cdiff[k])/(2.0*zcell[k]);
 		caa[k] = (xcell[0]*zcell[l]*cdiff[l] +
 			xcell[0]*zcell[k]*cdiff[k])/4.0;
 	}
 	// Lower left corner
	l = k-1;
	cae[k] = -(zcell[l]*cdiff[l])/(2.0*xcell[0]);
	cas[k] = 0.0;
	caa[k] = (xcell[0]*zcell[l]*cdiff[l])/4.0;
 	for (i=1; i<cx; i++){
 		j = i-1;
 		// Upper edge
		cae[i*nz] = -(zcell[0]*cdiff[i*cz])/(2.0*xcell[i]);
		cas[i*nz] = -(xcell[j]*cdiff[j*cz] + xcell[i]*cdiff[i*cz])
			/(2.0*zcell[0]);
		caa[i*nz] = (xcell[i]*zcell[0]*cdiff[i*cz] +
			xcell[j]*zcell[0]*cdiff[j*cz])/4.0;
		// Inner coefficients
 		for (k=1; k<cz; k++){
 			l = k-1;
 			cae[i*nz+k] = -(zcell[l]*cdiff[i*cz+l] + zcell[k]*cdiff[i*cz+k])
 				/(2.0*xcell[i]);
 			cas[i*nz+k] = -(xcell[j]*cdiff[j*cz+k] + xcell[i]*cdiff[i*cz+k])
 				/(2.0*zcell[k]);
 			caa[i*nz+k] = (xcell[j]*zcell[l]*cdiff[j*cz+l] +
 				xcell[i]*zcell[l]*cdiff[i*cz+l] +
 				xcell[i]*zcell[k]*cdiff[i*cz+k] +
 				xcell[j]*zcell[k]*cdiff[j*cz+k])/4.0;
 		}
 		// Lower edge
 		l = k-1;
 		cae[i*nz+k] = -(zcell[l]*cdiff[i*cz+l])/(2.0*xcell[i]);
 		cas[i*nz+k] = 0.0;
 		caa[i*nz+k] = (xcell[j]*zcell[l]*cdiff[j*cz+l] +
 			xcell[i]*zcell[l]*cdiff[i*cz+l])/4.0;
 	}
 	// Upper right corner
 	j = i-1;
	cas[i*nz] = -(xcell[j]*cdiff[j*cz])/(2.0*zcell[0]);
	caa[i*nz] = (xcell[j]*zcell[0]*cdiff[j*cz])/4.0;
 	// Right edge
	for (k=1; k<cz; k++){
		l = k-1;
		cas[i*nz+k] = -(xcell[j]*cdiff[j*cz+k])/(2.0*zcell[k]);
		caa[i*nz+k] = (xcell[j]*zcell[l]*cdiff[j*cz+l] +
			xcell[j]*zcell[k]*cdiff[j*cz+k])/4.0;
	}
 	// Lower right corner
	l = k-1;
	caa[i*nz+k] = (xcell[j]*zcell[l]*cdiff[j*cz+l])/4.0;
}

void geo2d::coupap(double lambda)
/*
 * Computes wavenumber dependent anomalous self-coupling coefficient
 * cap, memory was allocated at previous calls to coupaesa
 */
{
	
	double lambda2 = lambda*lambda;
	int i, j, k, l;

 	// Upper left corner
  cap[0] = -cae[0] - cas[0] + lambda2*caa[0];
 	// Left edge
 	for (k=1; k<cz; k++){
 		l = k-1;
 		cap[k] = -cae[k] - cas[k] - cas[l] + lambda2*caa[k];
 	}
 	// Lower left corner
 	l = k-1;
 	cap[k] = -cae[k] - cas[l] + lambda2*caa[k];
 	for (i=1; i<cx; i++){
 		j = i-1;
 		// Upper edge
 		cap[i*nz] = -cae[i*nz] - cas[i*nz] - cae[j*nz] + lambda2*caa[i*nz];
 		// Inner coefficients
 		for (k=1; k<cz; k++){
 			l = k-1;
 			cap[i*nz+k] = -cae[i*nz+k] - cas[i*nz+k]
 				- cae[j*nz+k] - cas[i*nz+l] + lambda2*caa[i*nz+k];
 		}
 		// Lower edge
 		l = k-1;
 		cap[i*nz+k] = -cae[i*nz+k] - cae[j*nz+k] - cas[i*nz+l]
 			+ lambda2*caa[i*nz+k];
 	}
 	j = i-1;
 	// Upper right corner
 	cap[i*nz] = -cae[j*nz] - cas[i*nz] + lambda2*caa[i*nz];
 	// Right edge
 	for (k=1; k<cz; k++){
 		l = k-1;
 		cap[i*nz+k] = -cae[j*nz+k] - cas[i*nz+k] - cas[i*nz+l]
 			+ lambda2*caa[i*nz+k];
 	}
 	// Lower right corner
 	l = k-1;
 	cap[i*nz+k] = -cae[j*nz+k] - cas[i*nz+l] + lambda2*caa[i*nz+k];
}

void geo2d::coupn(double lambda)
/*
 * Assembly of coupling matrix C, computes normal self-coupling
 * coefficient cnp in place.
 * Format of matrix C is LAPACK's banded spd format DPB:
 * m-th subdiagonal is stored in m-th row of C (m=0,...,nz).
 * Applies 'bound' boundary conditions to the anomalous potential:
 * 'd' or 'D' ... Dirichlet type
 * 'm' or 'M' ... mixed type a la Dey & Morrison 1979
 */
{
	double xc, zc, x, z, r, alf1, alf2, lambda2 = lambda*lambda;
	int i, j, k, l, m, n, band = nz+1;

  switch (config.Boundary)
  {
  /* Dirichlet boundary conditions */
  case 'd':
  	// Upper left corner
  	C[nz] = -2.0*cne[0] - cns[0] + lambda2*cna[0];
  	// Left edge
  	for (k=1; k<cz; k++){
  		l = k-1;
  		m = k*band+nz-k;
  		n = k*band+nz-1;
  		while (m<n){
  			C[m++]=0.0;
  		}
  		C[m++] = cns[l];
  		C[m++] = -2.0*cne[k] - cns[k] - cns[l] + lambda2*cna[k];
  	}
  	// Lower left corner
  	l = k-1;
  	m = k*band+nz-k;
  	n = k*band+nz-1;
 		while (m<n){
 			C[m++]=0.0;
 		}
 		C[m++] = cns[l];
  	C[m++] = -2.0*(cne[k] + cns[l]) + lambda2*cna[k];
  	for (i=1; i<cx; i++){
  		j = i-1;
  		// Upper edge
  		C[m++] = cne[j*nz];
  		n = i*nz*band+nz-1;
  		while (m<n){
  			C[m++] = 0.0;
  		}
  		C[m++] = 0.0;
  		C[m++] = -cne[i*nz] - cns[i*nz] - cne[j*nz] + lambda2*cna[i*nz];
  		// Inner coefficients
  		for (k=1; k<cz; k++){
  			l = k-1;
  			C[m++] = cne[j*nz+k];
  			n = (i*nz+k)*band+nz-1;
  			while (m<n){
  				C[m++] = 0.0;
  			}
  			C[m++] = cns[i*nz+l];
  			C[m++] = -cne[i*nz+k] - cns[i*nz+k]
  				- cne[j*nz+k] - cns[i*nz+l] + lambda2*cna[i*nz+k];
  		}
  		// Lower edge
  		l = k-1;
  		C[m++] = cne[j*nz+k];
  		n = (i*nz+k)*band+nz-1;
  		while (m<n){
  			C[m++] = 0.0;
  		}
  		C[m++] = cns[i*nz+l];
  		C[m++] = -cne[i*nz+k] - cne[j*nz+k] - 2.0*cns[i*nz+l]
  			+ lambda2*cna[i*nz+k];
  	}
  	j = i-1;
  	// Upper right corner
  	C[m++] = cne[j*nz];
  	n = i*nz*band+nz-1;
  	while (m<n){
  		C[m++] = 0.0;
  	}
  	C[m++] = 0.0;
  	C[m++] = -2.0*cne[j*nz] - cns[i*nz] + lambda2*cna[i*nz];
  	// Right edge
  	for (k=1; k<cz; k++){
  		l = k-1;
  		C[m++] = cne[j*nz+k];
  		n = (i*nz+k)*band+nz-1;
  		while (m<n){
  			C[m++] = 0.0;
  		}
  		C[m++] = cns[i*nz+l];
  		C[m++] = -2.0*cne[j*nz+k] - cns[i*nz+k] - cns[i*nz+l]
  			+ lambda2*cna[i*nz+k];
  	}
  	// Lower right corner
  	l = k-1;
  	C[m++] = cne[j*nz+k];
  	n = (i*nz+k)*band+nz-1;
  	while (m<n){
  		C[m++] = 0.0;
  	}
  	C[m++] = cns[i*nz+l];
  	C[m] = -2.0*(cne[j*nz+k] + cns[i*nz+l]) + lambda2*cna[i*nz+k];
  	break;
  /* mixed boundary conditions with respect to a source at the top-center */
  case 'm':
 		xc = (xnode[0]+xnode[nx-1])/2.0;
 		zc = znode[0];
 		// Upper left corner
 		x = fabs(xnode[0]-xc);
 		z = znode[0]-zc;
 		r = sqrt(x*x + z*z);
 		alf1 = lambda*r;
 		if (alf1 >= 670.0){
 			alf1 = lambda*xcell[0]*x/r;
 		}
 		else {
  		alf1 = lambda*xcell[0]*bessk1(alf1)/bessk0(alf1)*x/r;
  	}
 		C[nz] = -(1.0+alf1)*cne[0] - cns[0] + lambda2*cna[0];
  	// Left edge
  	for (k=1; k<cz; k++){
  		l = k-1;
  		m = k*band+nz-k;
  		n = k*band+nz-1;
  		while (m<n){
  			C[m++]=0.0;
  		}
  		C[m++] = cns[l];
 			z = znode[k]-zc;
 			r = sqrt(x*x + z*z);
 			alf1 = lambda*r;
 			if (alf1 >= 670.0){
 				alf1 = lambda*xcell[0]*x/r;
 			}
 			else {
 		 		alf1 = lambda*xcell[0]*bessk1(alf1)/bessk0(alf1)*x/r;
  		}
  		C[m++] = -(1.0+alf1)*cne[k] - cns[k] - cns[l] + lambda2*cna[k];
  	}
  	// Lower left corner
  	l = k-1;
  	m = k*band+nz-k;
  	n = k*band+nz-1;
 		while (m<n){
 			C[m++]=0.0;
 		}
 		C[m++] = cns[l];
 		z = fabs(znode[k]-zc);
 		r = sqrt(x*x + z*z);
 		alf1 = lambda*r;
 		alf2 = alf1;
 		if (alf1 >= 670.0){
 			alf1 = lambda*xcell[0]*x/r;
 			alf2 = lambda*zcell[l]*z/r;
 		}
 		else {
  		alf1 = lambda*xcell[0]*bessk1(alf1)/bessk0(alf1)*x/r;
  		alf2 = lambda*zcell[l]*bessk1(alf2)/bessk0(alf2)*z/r;
  	}
  	C[m++] = -(1.0+alf1)*cne[k] - (1.0+alf2)*cns[l] + lambda2*cna[k];
  	for (i=1; i<cx; i++){
  		j = i-1;
  		// Upper edge
  		C[m++] = cne[j*nz];
  		n = i*nz*band+nz-1;
  		while (m<n){
  			C[m++] = 0.0;
  		}
  		C[m++] = 0.0;
  		C[m++] = -cne[i*nz] - cns[i*nz] - cne[j*nz] + lambda2*cna[i*nz];
  		// Inner coefficients
  		for (k=1; k<cz; k++){
  			l = k-1;
  			C[m++] = cne[j*nz+k];
  			n = (i*nz+k)*band+nz-1;
  			while (m<n){
  				C[m++] = 0.0;
  			}
  			C[m++] = cns[i*nz+l];
  			C[m++] = -cne[i*nz+k] - cns[i*nz+k]
  				- cne[j*nz+k] - cns[i*nz+l] + lambda2*cna[i*nz+k];
  		}
  		// Lower edge
  		l = k-1;
  		C[m++] = cne[j*nz+k];
  		n = (i*nz+k)*band+nz-1;
  		while (m<n){
  			C[m++] = 0.0;
  		}
  		C[m++] = cns[i*nz+l];
 			x = xnode[i]-xc;
 			r = sqrt(x*x + z*z);
 			alf2 = lambda*r;
 			if (alf2 >= 670.0){
 				alf2 = lambda*znode[l]*z/r;
	 		}
 			else {
  			alf2 = lambda*znode[l]*bessk1(alf2)/bessk0(alf2)*z/r;
	  	}
  		C[m++] = -cne[i*nz+k] - cne[j*nz+k] - (1.0+alf2)*cns[i*nz+l]
  			+ lambda2*cna[i*nz+k];
  	}
  	j = i-1;
  	// Upper right corner
  	C[m++] = cne[j*nz];
  	n = i*nz*band+nz-1;
  	while (m<n){
  		C[m++] = 0.0;
  	}
  	C[m++] = 0.0;
  	x = fabs(xnode[i]-xc);
		z = znode[0]-zc;
		r = sqrt(x*x + z*z);
		alf1 = lambda*r;
		if (alf1 >= 670.0){
			alf1 = lambda*xcell[j]*x/r;
		}
		else {
 			alf1 = lambda*xcell[j]*bessk1(alf1)/bessk0(alf1)*x/r;
  	}
  	C[m++] = -(1.0+alf1)*cne[j*nz] - cns[i*nz] + lambda2*cna[i*nz];
  	// Right edge
  	for (k=1; k<cz; k++){
  		l = k-1;
  		C[m++] = cne[j*nz+k];
  		n = (i*nz+k)*band+nz-1;
  		while (m<n){
  			C[m++] = 0.0;
  		}
  		C[m++] = cns[i*nz+l];
			z = znode[k]-zc;
			r = sqrt(x*x + z*z);
			alf1 = lambda*r;
			if (alf1 >= 670.0){
				alf1 = lambda*xcell[j]*x/r;
			}
			else {
 				alf1 = lambda*xcell[j]*bessk1(alf1)/bessk0(alf1)*x/r;
  		}
  		C[m++] = -(1.0+alf1)*cne[j*nz+k] - cns[i*nz+k] - cns[i*nz+l]
  			+ lambda2*cna[i*nz+k];
  	}
  	// Lower right corner
  	l = k-1;
  	C[m++] = cne[j*nz+k];
  	n = (i*nz+k)*band+nz-1;
  	while (m<n){
  		C[m++] = 0.0;
  	}
  	C[m++] = cns[i*nz+l];
		z = fabs(znode[k]-zc);
		r = sqrt(x*x + z*z);
		alf1 = lambda*r;
		alf2 = alf1;
		if (alf1 >= 670.0){
			alf1 = lambda*xcell[j]*x/r;
			alf2 = lambda*zcell[l]*z/r;
		}
		else {
 			alf1 = lambda*xcell[j]*bessk1(alf1)/bessk0(alf1)*x/r;
 			alf2 = lambda*zcell[l]*bessk1(alf2)/bessk0(alf2)*z/r;
  	}
  	C[m] = -(1.0+alf1)*cne[j*nz+k] - (1.0+alf2)*cns[i*nz+l]
  		+ lambda2*cna[i*nz+k];
  	break;
  }
}

