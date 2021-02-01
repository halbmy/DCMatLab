/*
 * righthandside.cpp
 */

#include "geo2d.h"

void geo2d::righthandside(double lambda)
{
	double x, z, h;
	double *v, *vb, *bi, *cp, *vp, *ce, *ve, *cs, *vs, *cn, *vn, *cw, *vw, *vbb;
	int isrc, i, j, n = nx*nz, nb = 2*nz+nx;
	
  v = (double *) CALLOC(n, sizeof(double));
  vb = (double *) CALLOC(nb, sizeof(double));	

  bi=b;
	for (i=n*nsrc; i; i--){	// Set b=0 in case of cancellation of digits
		*bi=0.0;							// to get equally biased results
		bi++;
	}

	for (isrc=0; isrc<nsrc; isrc++)
  {
    h = zsrc[isrc] - znode[0];
		// First, compute normal potential in wavenumber domain
		// all grid nodes	
    for (i=0; i<nx; i++)
    {
      x = xnode[i] - xsrc[isrc];
		  for (j=0; j<nz; j++)
      {
        z = znode[j] - znode[0];
			  v[i*nz+j] = 2.0 * potential(isrc, h, x, z, lambda);
      }
    }
    
    // and normal derivative of potential at the boundary nodes
    // Left boundary
    x = xnode[0] - xsrc[isrc];
    for (j=0; j<nz; j++)
    {
      z = znode[j] - znode[0];
      vb[j] = -2.0 * xcell[0] * dpotdx(isrc, h, x, z, lambda);
    }
    // Lower boundary
    z = znode[cz] - znode[0];
    for (i=0; i<nx; i++)
    {
      x = xnode[i] - xsrc[isrc];
      vb[j++] = 2.0 * zcell[cz-1] * dpotdz(isrc, h, x, z, lambda);
    }
    // Right boundary
    x = xnode[cx] - xsrc[isrc];
    for (i=0; i<nz; i++)
    {
      z = znode[i] - znode[0];
      vb[j++] = 2.0 * xcell[cx-1] * dpotdx(isrc, h, x, z, lambda);
    }
    
    // Second, matrix Ca * vector v + boundary terms
    bi = b+isrc*n; vbb = vb;
    cp = cap; vp = v;
		ce = cae; ve = v+nz;
		cs = cas; vs = v+1;
		// Upper left corner
		*bi = (*cp)*(*vp)
			+ (*cs)*(*vs)
			+ (*ce)*((*ve)+(*vbb));
		bi++; vbb++;
		cp++; vp++;
		ce++; ve++;
		cs++; vs++;
		cn = cas; vn = v;
		// Left boundary
		for (j=1; j<cz; j++){
			*bi = (*cn)*(*vn)
				+ (*cp)*(*vp)
				+ (*cs)*(*vs)
				+ (*ce)*((*ve)+(*vbb));
			bi++; vbb++;
			cp++; vp++;
			ce++; ve++;
			cs++; vs++;
			cn++; vn++;
		}
		// Lower left corner
		*bi = (*cp)*(*vp)
			+ (*cn)*((*vn)+(*(vbb+1)))
			+ (*ce)*((*ve)+(*vbb));
		bi++; vbb+=2;
		cp++; vp++;
		ce++; ve++;
		cs++; vs++;
		cn++; vn++;
		cw = cae; vw = v;
		for (i=1; i<cx; i++){
			// Upper boundary
			*bi = (*ce)*(*ve)
				+ (*cs)*(*vs)
				+ (*cp)*(*vp)
				+ (*cw)*(*vw);
			bi++;
			cp++; vp++;
			ce++; ve++;
			cs++; vs++;
			cn++; vn++;
			cw++; vw++;
			// Inner grid nodes
			for (j=1; j<cz; j++){
				*bi = (*ce)*(*ve)
					+ (*cs)*(*vs)
					+ (*cp)*(*vp)
					+ (*cn)*(*vn)
					+ (*cw)*(*vw);
				bi++;
				cp++; vp++;
				ce++; ve++;
				cs++; vs++;
				cn++; vn++;
				cw++; vw++;
		}
			// Lower boundary
			*bi = (*ce)*(*ve)
				+ (*cp)*(*vp)
				+ (*cn)*((*vn)+(*vbb))
				+ (*cw)*(*vw);
			bi++; vbb++;
			cp++; vp++;
			if (i<cx-1) {ce++; ve++;}
			cs++; vs++;
			cn++; vn++;
			cw++; vw++;
	}
		vbb++; // vb[nx+nz-1] skipped for call at lower right corner
		// Upper right corner
		*bi = (*cp)*(*vp)
			+ (*cs)*(*vs)
			+ (*cw)*((*vw)+(*vbb));
		bi++; vbb++;
		cp++; vp++;
		cs++; vs++;
		cn++; vn++;
		cw++; vw++;
		// Right boundary
		for (j=1; j<cz; j++){
			*bi = (*cn)*(*vn)
				+ (*cp)*(*vp)
				+ (*cs)*(*vs)
				+ (*cw)*((*vw)+(*vbb));
			bi++; vbb++;
			cp++; vp++;
			if (j<cz-1) {cs++; vs++;}
			cn++; vn++;
			cw++; vw++;
		}
		// Lower right corner
		*bi = (*cp)*(*vp)
			+ (*cn)*((*vn)+(*vbb))
			+ (*cw)*((*vw)+(*(vbb-nz)));
	}
	FREE(v);
	FREE(vb);
}

void geo2d::righthandside(double lambda, int isrc)
{
	double x, z, h;
	double *v, *vb, *bi, *cp, *vp, *ce, *ve, *cs, *vs, *cn, *vn, *cw, *vw, *vbb;
	int i, j, n = nx*nz, nb = 2*nz+nx;
 
	v = (double *) CALLOC(n, sizeof(double));	
	vb = (double *) CALLOC(nb, sizeof(double));

	bi=b+isrc*n;
	for (i=n; i; i--){	// Set b=0 in case of cancellation of digits
    *bi=0.0;					// to get equally biased results
    bi++;
  }
  
  h = zsrc[isrc] - znode[0];
  // First, compute normal potential in wavenumber domain
	// all grid nodes	
  for (i=0; i<nx; i++)
  {
    x = xnode[i] - xsrc[isrc];
		for (j=0; j<nz; j++)
    {
      z = znode[j] - znode[0];
      v[i*nz+j] = 2.0 * potential(isrc, h, x, z, lambda);
    }
  }

  // and normal derivative of potential at the boundary nodes
	// Left boundary
  x = xnode[0] - xsrc[isrc];
  for (j=0; j<nz; j++)
  {
    z = znode[j] - znode[0];
    vb[j] = -2.0 * xcell[0] * dpotdx(isrc, h, x, z, lambda);
  }
  // Lower boundary
  z = znode[cz] - znode[0];
  for (i=0; i<nx; i++)
  {
    x = xnode[i] - xsrc[isrc];
    vb[j++] = 2.0 * zcell[cz-1] * dpotdz(isrc, h, x, z, lambda);
  }
  // Right boundary
  x = xnode[cx] - xsrc[isrc];
  for (i=0; i<nz; i++)
  {
    z = znode[i] - znode[0];
    vb[j++] = 2.0 * xcell[cx-1] * dpotdx(isrc, h, x, z, lambda);
  }
  
  // Second, matrix Ca * vector v + boundary terms
  bi = b+isrc*n; vbb = vb;
  cp = cap; vp = v;
  ce = cae; ve = v+nz;
	cs = cas; vs = v+1;
	// Upper left corner
	*bi = (*cp)*(*vp)
		+ (*cs)*(*vs)
		+ (*ce)*((*ve)+(*vbb));
	bi++; vbb++;
	cp++; vp++;
	ce++; ve++;
	cs++; vs++;
	cn = cas; vn = v;
	// Left boundary
	for (j=1; j<cz; j++){
		*bi = (*cn)*(*vn)
			+ (*cp)*(*vp)
			+ (*cs)*(*vs)
			+ (*ce)*((*ve)+(*vbb));
		bi++; vbb++;
		cp++; vp++;
		ce++; ve++;
		cs++; vs++;
		cn++; vn++;
	}
	// Lower left corner
	*bi = (*cp)*(*vp)
		+ (*cn)*((*vn)+(*(vbb+1)))
		+ (*ce)*((*ve)+(*vbb));
	bi++; vbb+=2;
	cp++; vp++;
	ce++; ve++;
	cs++; vs++;
	cn++; vn++;
	cw = cae; vw = v;
	for (i=1; i<cx; i++){
		// Upper boundary
		*bi = (*ce)*(*ve)
			+ (*cs)*(*vs)
			+ (*cp)*(*vp)
			+ (*cw)*(*vw);
		bi++;
		cp++; vp++;
		ce++; ve++;
		cs++; vs++;
		cn++; vn++;
		cw++; vw++;
		// Inner grid nodes
		for (j=1; j<cz; j++){
			*bi = (*ce)*(*ve)
				+ (*cs)*(*vs)
				+ (*cp)*(*vp)
				+ (*cn)*(*vn)
				+ (*cw)*(*vw);
			bi++;
			cp++; vp++;
			ce++; ve++;
			cs++; vs++;
			cn++; vn++;
			cw++; vw++;
		}
		// Lower boundary
		*bi = (*ce)*(*ve)
			+ (*cp)*(*vp)
			+ (*cn)*((*vn)+(*vbb))
			+ (*cw)*(*vw);
		bi++; vbb++;
		cp++; vp++;
		if (i<cx-1) {ce++; ve++;}
		cs++; vs++;
		cn++; vn++;
		cw++; vw++;
	}
	vbb++; // vb[nx+nz-1] skipped for call at lower right corner
	// Upper right corner
	*bi = (*cp)*(*vp)
		+ (*cs)*(*vs)
		+ (*cw)*((*vw)+(*vbb));
	bi++; vbb++;
	cp++; vp++;
	cs++; vs++;
	cn++; vn++;
	cw++; vw++;
	// Right boundary
	for (j=1; j<cz; j++){
		*bi = (*cn)*(*vn)
			+ (*cp)*(*vp)
			+ (*cs)*(*vs)
			+ (*cw)*((*vw)+(*vbb));
		bi++; vbb++;
		cp++; vp++;
		if (j<cz-1) {cs++; vs++;}
		cn++; vn++;
		cw++; vw++;
	}
	// Lower right corner
	*bi = (*cp)*(*vp)
		+ (*cn)*((*vn)+(*vbb))
		+ (*cw)*((*vw)+(*(vbb-nz)));
	FREE(v);
	FREE(vb);
}

