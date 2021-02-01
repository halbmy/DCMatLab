/* 
 * grid.cpp
 */

#include "geo2d.h"

void geo2d::grid()
/*
 * Refine model and prolong in +/- x and +z - direction.
 * nprol ... number of additional cells in each direction (global) 
 * fprol ... factor stretching boundary cell sizes
 * mprol ... additional cells get
 *  'a' ..... average conductivity from all area weighted model cells
 *  'b' ..... conductivities of old boundary cells
 * nref .... number of grid nodes which are added within a cell in 
 *           horizontal/vertical direction to refine the grid (global)
 * mback ... choice of background conductivity for singularity removal
 *  'a' ..... average conductivities of the whole grid
 *  'm' ..... average conductivity of each source location
 *  's' ..... conductivity of source location - exact singularity removal
 */
{
	double f, g;
	int i, j, k, l, m, n, offsetx, offsetz, nrx, nrz;
   static bool outputflag = true;
	
   offsetx = config.NoProlX + config.NoEquiX;
   offsetz = config.NoProlZ + config.NoEquiZ;
   nrx = config.NoRefineX + 1;
   nrz = config.NoRefineZ + 1;

	cx = nrx * nxcnd + 2 * (config.NoEquiX + config.NoProlX);  // Compute size of prolonged arrays
	cz = nrz * nzcnd + config.NoEquiZ + config.NoProlZ;
	nx = cx + 1;
	nz = cz + 1;
   if (outputflag)
   {
      //cout << "GEO2D: input grid size is " << nxcnd << " x " << nzcnd << " cells." << endl
      //  << "       size of refined and prolonged grid is " << cx << " x " << cz << " cells." << endl;
      outputflag = false;
   }
   if (chgrd)
   {
      if (xcell != NULL) FREE(xcell);
      if (zcell != NULL) FREE(zcell);
      if (xnode != NULL) FREE(xnode);
      if (znode != NULL) FREE(znode);
      if (cnorm != NULL) FREE(cnorm);
      if (cdiff != NULL) FREE(cdiff);
      if (cback != NULL)
      {
         for (i = 0; i < nsrc; i++)
         {
            FREE(cback[i]);
         }
         FREE(cback);
      }      
      if (iback != NULL)
      {
         for (i = 0; i < nsrc; i++)
         {
            FREE(iback[i]);
         }
         FREE(iback);
      }
      xcell = (double *) CALLOC(cx, sizeof(double));			// Allocate memory for arrays
      zcell = (double *) CALLOC(cz, sizeof(double));
      xnode = (double *) CALLOC(nx, sizeof(double));
      znode = (double *) CALLOC(nz, sizeof(double));
      cnorm = (double *) CALLOC(cx*cz, sizeof(double));
      cdiff = (double *) CALLOC(cx*cz, sizeof(double));
      cback = (double **) CALLOC(nsrc, sizeof(double*));
      for (i = 0; i < nsrc; i++)
      {
         cback[i] = (double *) CALLOC(4, sizeof(double));
      }
      iback = (int **) CALLOC(nsrc, sizeof(int *));
      for (i = 0; i < nsrc; i++)
      {
         iback[i] = (int *) CALLOC(2, sizeof(int));
      }
   }

	/* Copy model to right positions in vectors xnode, znode and array cnorm */
	for (i=0,k=0; i<nzcnd; i++)
   {
		znode[k++] = zgrd[i];
		for (j=1; j<nrz; j++)
      {
         znode[k++] = zgrd[i] + (zgrd[i+1]-zgrd[i])*double(j)/double(nrz);
		}
	}
	znode[k] = zgrd[i];

	for (i=0, k=offsetx; i<nxcnd; i++)
   {
		xnode[k++] = xgrd[i];
		for (j=1; j<nrx; j++)
      {
			xnode[k++] = xgrd[i] + (xgrd[i+1]-xgrd[i])*double(j)/double(nrx);
		}
	}
	xnode[k] = xgrd[i];

   for (i=0, k=offsetx; i < nxcnd; i++, k+=nrx)
   {
      for (j=0, l=0; j<nzcnd; j++, l+=nrz)
      {
         for (m = 0; m < nrx; m++)
         {
            for (n = 0; n < nrz; n++)
            {
               cnorm[(k+m)*cz+l+n] = cond[i*nzcnd+j];
            }
         }
      }
   }

   /* Calculate cell sizes and average conductivity */
   f = 0.0;
   for (i=0; i<cz-offsetz; i++)
   {
      zcell[i] = znode[i+1]-znode[i];
   }
   for (i=offsetx; i<cx-offsetx; i++)
   {
      xcell[i] = xnode[i+1]-xnode[i];
      for (j=0; j<cz-offsetz; j++)
      {
         f += xcell[i] * zcell[j] * cnorm[i*cz+j];
      }
   }
   f /= ( (xnode[cx-offsetx] - xnode[offsetx])
      * (znode[cz-offsetz] - znode[0]) );
   
   /* Add prolonging values */
   j=offsetx-1;
   k=cx-offsetx;
   for (i=0; i<config.NoEquiX; i++,j--,k++)
   {
      xcell[j] = xcell[j+1];
      xnode[j] = xnode[j+1] - xcell[j];
      xcell[k] = xcell[k-1];
      xnode[k+1] = xnode[k] + xcell[k];
   }
   for (i=0; i<config.NoProlX; i++,j--,k++)
   {
      xcell[j] = config.FacProlX * xcell[j+1];
      xnode[j] = xnode[j+1] - xcell[j];
      xcell[k] = config.FacProlX * xcell[k-1];
      xnode[k+1] = xnode[k] + xcell[k];
   }
   switch (config.ModProlong)
   {
   case 'a':
      k = cx-1;
      for (j=0; j<offsetx; j++,k--)
      {
         for (i=0; i<cz-offsetz; i++)
         {
            cnorm[j*cz+i] = f;
            cnorm[k*cz+i] = f;
         }
      }
      break;
   case 'b':
      for (i=0; i<cz-offsetz; i++)
      {
         f = cnorm[offsetx*cz+i];
         g = cnorm[(cx-offsetx-1)*cz+i];
         k = cx-1;
         for (j=0; j<offsetx; j++,k--)
         {
            cnorm[j*cz+i] = f;
            cnorm[k*cz+i] = g;
         }
      }
      break;
   }

   k=cz-offsetz;
   for (i=0; i<config.NoEquiZ; i++,k++)
   {
      zcell[k] = zcell[k-1];
      znode[k+1] = znode[k] + zcell[k];
   }
   for (i=0; i<config.NoProlZ; i++,k++)
   {
      zcell[k] = config.FacProlZ * zcell[k-1];
      znode[k+1] = znode[k] + zcell[k];
   }
   switch (config.ModProlong)
   {
   case 'a':
      for (j=0; j<cx; j++)
      {
         for (i=cz-offsetz; i<cz; i++)
         {
            cnorm[j*cz+i] = f;
         }
      }
      break;
   case 'b':
      for (j=0; j<cx; j++)
      {
         f = cnorm[(j+1)*cz-offsetz-1];
         for (i=cz-offsetz; i<cz; i++)
         {
            cnorm[j*cz+i] = f;
         }
      }
      break;
   }
   
  /* Get conductivities and indices of source locations */
  switch (config.Background)
  {
  case 'a':        /* Take average conductivity of the whole grid */
    for(i=0; i<nsrc; i++)
    {
      if (fabs(zsrc[i]-zgrd[0]) < 1.0e-15) // surface source
      {
        cback[i][0] = cback[i][1] = 0.0;
        cback[i][2] = cback[i][3] = f;
      }
      else        // buried source
      {
        cback[i][0] = cback[i][1] = cback[i][2] = cback[i][3] = f;
      }
      iback[i][0] = iback[i][0] = 0;
    }
    break;
  case 'm':        /* Take average conductivity of each source location */
    for(i=0; i<nsrc; i++)
    {
      switch (sourcepos(i, &k, &l))
      {
      case 'n':
        // source on an interiour grid node
        if (l > 0) // buried source
        {
          f = (cnorm[(k-1)*cz+l-1] + cnorm[k*cz+l-1]
            + cnorm[k*cz+l] + cnorm[(k-1)*cz+l]) / 4.0;
          cback[i][0] = cback[i][1] = cback[i][2] = cback[i][3] = f;
        }
        else       // surface source
        {
          f = (cnorm[k*cz+l] + cnorm[(k-1)*cz+l]) / 2.0;
          cback[i][2] = cback[i][3] = f;
          cback[i][0] = cback[i][1] = 0.0;
        }
        break;
      case 'h':
        // source on a horizontal edge of cells
        if (l > 0) // buried source
        {
          f = (cnorm[k*cz+l-1] + cnorm[k*cz+l]) / 2.0;
          cback[i][0] = cback[i][1] = cback[i][2] = cback[i][3] = f;
        }
        else       // surface source
        {
          cback[i][0] = cback[i][1] = 0.0;
          cback[i][2] = cback[i][3] = cnorm[k*cz+l];
        }
        break;
      case 'v':
        // source on a vertical edge of cells
        f = (cnorm[k*cz+l] + cnorm[(k-1)*cz+l]) / 2.0;
        cback[i][0] = cback[i][1] = cback[i][2] = cback[i][3] = f;
        break;
      case 'i':
        // source in the interiour of a cell
        cback[i][0] = cback[i][1] = cback[i][2] = cback[i][3] = cnorm[k*cz+l];
        break;
      }
      iback[i][0] = k;
      iback[i][1] = l;
    }
    break;
  case 's':        /* Take conductivity of each source location */
    for(i=0; i<nsrc; i++)
    {
      switch (sourcepos(i, &k, &l))
      {
      case 'n':
        // source on a grid node
        if (l > 0) // buried source
        {
          cback[i][0] = cnorm[(k-1)*cz+l-1];
          cback[i][1] = cnorm[k*cz+l-1];
        }
        else        // surface source
        {
          cback[i][0] = cback[i][1] = 0.0;
        }
        cback[i][2] = cnorm[k*cz+l];
        cback[i][3] = cnorm[(k-1)*cz+l];
        break;
      case 'h':
        // source on a horizontal edge of cells
        if (l > 0) // buried source
        {
          cback[i][0] = cback[i][1] = cnorm[k*cz+l-1];
        }
        else        // surface source
        {
          cback[i][0] = cback[i][1] = 0.0;
        }
        cback[i][2] = cback[i][3] = cnorm[k*cz+l];
        break;
      case 'v':
        // source on a vertical edge of cells
        cback[i][0] = cback[i][3] = cnorm[(k-1)*cz+l];
        cback[i][1] = cback[i][2] = cnorm[k*cz+l];
        break;
      case 'i':
        // source in the interiour of a cell
        cback[i][0] = cback[i][1] = cback[i][2] = cback[i][3] = cnorm[k*cz+l];
        break;
      }
      iback[i][0] = k;
      iback[i][1] = l;
    }
    break;
  }
}

