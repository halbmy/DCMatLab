/* evaluation for 2d inversion */

#include "dc2d.h"

/* void app_setup() */
int application::setup(char *dataname, char *modelname,
               double key1, int key2, double xmin, double xmax)
{
  int i, n, cx, cz;
  double t0, t1;
   int NoSources, NCond;
   double *XSources, *ZSources;
   int **ArSources;
  /* 
     set key for first objective function: 
     data misfit measured by l_key1-norm
     1 ... l1-norm
     2 ... l2-norm (euclidian norm)
     3 ... linf-norm (maximum-norm)
  */
  if (key1 <= 0.0)
  {
    obj1 = 2.0;
    fprintf(stderr, "Invalid key value for 1st objective function. "
      "Reset to %g.\n", obj1);
  }
  else
  {
    obj1 = key1;
  }

  /* 
     set key for second objective function
     1 ... standard deviation of model vector
     2 ... number of discontinuities
     3 ... 1-norm of 1st derivative
  */
  if (key2 < 1 || key2 > 5)
  {
     if (key2 == 2)
        fprintf(stderr, "Model restriction #2 currently unavailable.\n");
    obj2 = 3;
    fprintf(stderr, "Invalid key value for 2nd objective function. "
      "Reset to %d.\n", obj2);
  }
  else
  {
    obj2 = key2;
  }

  /* Range of searched conductivity values */
  if (xmin < xmax)
  {
     C_min = pow(10.0, xmin);
     C_max = pow(10.0, xmax);
  }
  else
  {
     fprintf(stderr, "[%g, %g] is an invalid variable range.\n", xmin, xmax);
     exit(0);
  }

  /*
     get measured data:
     NoSources, XSources, ZSources,
     NoData, ISource1, ISource2, IPotential1,
     IPotential2, DataVal, DataStd
  */
  read2ddata(dataname, &NoSources, &XSources, &ZSources,
             &NoData, &ArSources, &DataVal, &DataStd);
#ifdef LOGRHOA
  for (i = 0; i < NoData; i++)
  {
	  t0 = DataStd[i] / DataVal[i];
     DataVal[i] = log(DataVal[i]);
     DataStd[i] = log(1.0 + t0);
  }
#endif

  /* Create array for synthetic data */
  DataSyn = new double [NoData];
  if (DataSyn == NULL)
  {
    fprintf(stderr, "Error allocating memory for DataSyn.\n");
    exit(0);
  }

  /*
     get grid and conductivity blocks for inversion:
     NXGrid, NZGrid, XGrid, ZGrid, IConductivity
  */
  read2dmodel(modelname, &NXGrid, &XGrid, &NZGrid, &ZGrid, &IConductivity);
  cx = NXGrid-1;
  cz = NZGrid-1;
  n = cx*cz;
  Conductivity = new double [n];
  NCond = 0;
  for (i = 0; i < n; i++)
  {
     NCond = max(NCond,IConductivity[i]);
     IConductivity[i]--; // since values are used as array indices
  }

  switch (obj2)
  {
  case 1:
     CXM = new double [cx]; // C**[0] is never used
     CXP = new double [cx];
     CZM = new double [cz];
     CZP = new double [cz];
     t1 = 0.5 * (XGrid[cx] - XGrid[cx-2]);
     CXM[cx-1] = 1.0 / (t1 * t1);
     // CXP[cx-1] = 0.0;
     for (i = cx-2; i > 0; i--)
     {
        t0 = 0.5 * (XGrid[i+1] - XGrid[i-1]);
        CXM[i] = 2.0 / (t0 * (t0 + t1));
        CXP[i] = 2.0 / (t1 * (t0 + t1));
        t1 = t0;
     }
     //  CXM[i] = 0.0;
     CXP[i] = 1.0 / (t1 * t1);

     t1 = 0.5 * (ZGrid[cz] - ZGrid[cz-2]);
     CZM[cz-1] = 1.0 / (t1 * t1);
     // CZP[cz-1] = 0.0;
     for (i = cz-2; i > 0; i--)
     {
        t0 = 0.5 * (ZGrid[i+1] - ZGrid[i-1]);
        CZM[i] = 2.0 / (t0 * (t0 + t1));
        CZP[i] = 2.0 / (t1 * (t0 + t1));
        t1 = t0;
     }
     //  CZM[i] = 0.0;
     CZP[i] = 1.0 / (t1 * t1);
     break;
  case 2:
     /* topology of conductivity blocks */
     findblocks(cx,cz,IConductivity);
     break;
  }

  /* Initialize forward modelling object */
  geo2dobj = new geo2d;
  geo2dobj->setmodel(NXGrid, XGrid, NZGrid, ZGrid, Conductivity);
  geo2dobj->setsources(NoSources, XSources, ZSources, NoData, ArSources);

  /* Free temporary storage */
  delete [] ZSources;
  delete [] XSources;
  delete [] ArSources;

  /* Return number of conductivity values to be found */
  return(NCond);
}
