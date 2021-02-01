/* evaluation for 2d inversion */

#include "dc2d.h"

void application::report(double *x, int nx, double *f, int nf, double *c, int nc, int gen, FILE *file)
{
  int i;

  /* Print table of block resitivities and function values */
  /* header */
  if (gen >= 0)
  {
	if (gen)
	  fprintf(file, "#\n# Pareto front at generation no. %d\n", gen);
	else
	  fprintf(file, "#\n# Pareto front at initial generation\n");
    fprintf(file, "# "); /* for gnuplot-ability */
    if (nf > 0)
      fprintf(file, "objective functions (%d), ", nf);
    if (nc > 0) 
      fprintf(file, "constraints (%d), ", nc);
    fprintf(file, "block resistivities (%d)\n", nx);
  }
  /* body of table */
  for (i = 0; i < nf; i++)
  {
    fprintf(file, "%g\t", f[i]);
  }
  for (i = 0; i < nc; i++)
  {
    fprintf(file, "%g\t", c[i]);
  }
  for (i = 0; i < nx-1; i++)
  {
    fprintf(file, "%g\t", pow(10.0, -x[i]));
  }
  fprintf(file, "%g\n", pow(10.0, -x[nx-1]));

/* Print function values and matrix of conductivities */
/*  int i, j, cx, cz, n;
  cx = NXGrid-1;
  cz = NZGrid-1;
  n = cx * cz;
  for(i = 0; i < n; i++)
  {
    j = IConductivity[i];
    if(j < nx)
    {
      Conductivity[i] = pow(10.0, x[j]);
    }
    else
    {
      fprintf(stderr, "Index for conductivity exceeds number of parameters "
        "provided by the GA.\n");
      exit(0);
    }
  }

  if(nf != 2)
  {
    fprintf(stderr, "Number of functions in app_eval doesn't match number of "
      "functions nfunc considered by the GA.\n");
    exit(0);
  }

  fprintf(file, "%g %g\n", f[0], f[1]);
  for (i = 0; i < cz; i++)
  {
    for (j = 0; j < cx; j++)
    {
      fprintf(file, "%g ", Conductivity[j*cz+i]);
    }
    fprintf(file, "\n");
  }
  fprintf(file, "\n");*/
}

