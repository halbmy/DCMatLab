#include "nsga2.h"

/* Standard constructor:
   Read GA-parameters from stdin */
nsga2::nsga2()
{
  long iseed = -1;
  verb = 1;
  app = NULL;

  /* Read GA-parameters from stdin */
  input();
  
  /* Allocate memory for population structure */
  memory();

  /*Set up linear transformation constants for objective function values using default arguments*/
  transf();

  /*Initialize the random no generator*/
  ran1(&iseed);
  
  /* Some counters */
  nmut = 0;
  ncross = 0;
  neval = 0;
}

nsga2::nsga2(int verbose)
{
  long iseed = -1;
  verb = verbose;
  app = NULL;

  /* Read GA-parameters from stdin */
  input();
  
  /* Allocate memory for population structure */
  memory();

  /*Set transform matrix for objective function values using default arguments*/
  transf();

  /*Initialize the random no generator*/
  ran1(&iseed);
  
  /* Some counters */
  nmut = 0;
  ncross = 0;
  neval = 0;
}

/* Extended constructor:
   Get GA-parameters from arguments of constructor */
nsga2::nsga2(int nind, int ngen, int ox, double px, double dx,
             double pm, double dm, double r, int nreal, int nbin,
             int bits, double vmin, double vmax, int vfix, int nf,
             int nc, int nsf, int nsp, long is, application *appl)
{
  int i, m;
  verb = 0;
  app = NULL;

  /* Get GA-parameters from arguments of constructor */
  /* number of real/binary encoded variables */
  nvar = nreal;
  nchrom = nbin;
  /*Allocate memory depending on nvar and nchrom:
    number of bits of binary parameters*/
  vlen = (int *) malloc( nchrom * sizeof(int) );
  if (vlen == NULL)
  {
    fprintf(stderr, "Error allocating memory for 'vlen'.\n");
    exit(1);
  }
  /*boundaries of binary parameters*/
  lim_b = (double **) malloc( nchrom * sizeof(double *) );
  if (lim_b == NULL)
  {
    fprintf(stderr, "Error allocating memory for 'lim_b'.\n");
    exit(1);
  }
  for (i = 0; i < nchrom; i++)
  {
    lim_b[i] = (double *) malloc( 2 * sizeof(double) );
    if (lim_b[i] == NULL)
    {
      fprintf(stderr, "Error allocating memory for 'lim_b[%d]'.\n", i);
      exit(1);
    }
  }
  /*boundaries of real parameters*/
  lim_r = (double **) malloc( nvar * sizeof(double *) );
  if (vlen == NULL)
  {
    fprintf(stderr, "Error allocating memory for 'lim_r'.\n");
    exit(1);
  }
  for (i = 0; i < nvar; i++)
  {
    lim_r[i] = (double *) malloc( 2 * sizeof(double) );
    if (lim_r[i] == NULL)
    {
      fprintf(stderr, "Error allocating memory for 'lim_r[%d]'.\n", i);
      exit(1);
    }
  }
  
  /*Asks for number of the functions*/
  nfunc = nf;

  /*Asks for number of constraints*/
  ncons = nc;
  
  /*Asks for number of the individuals in the population*/
  /*Default value = 30 to 100*/
  popsize = nind;
  popsize2 = 2 * popsize;
  
  /*No. of generations for which the GA will let the population evolve
    Default value is 100
    Too large value will take very long time and very small value will
    not let the GA reach the global Pareto front to the problem*/
  generi = 0;
  gener = ngen;

	/*Ask for control parameter for elitism
	  permissible values are 0 (uncontrolled elitsm)
	  to <1 (no elitism)*/
  relite = r;
	if (relite < 0.0) relite = 0.0;
	if (relite >= 1.0) relite = 0.999;

	/*Ask for length of sequences for sharing in fitness and parameter
	  space, respectively: nsharefit and nsharepar generations will be
	  shared with respect to distances in fitness and parameter space
	  switching between both*/
  nsharefit = nsf;
  nsharepar = nsp;
	if (nsharefit < 0 || nsharepar < 0)
	{
		nsharefit = 1;
		nsharepar = 0;
		fprintf(stderr, "Invalid values for nsharefit and nsharepar. "
			              "Set to %d and %d.\n", nsharefit, nsharepar);
	}
	nshare = nsharefit + nsharepar;

  /*No. of generations for which the GA will let the population evolve
    Default value is 100
    Too large value will take very long time and very small value will
    not let the GA reach the global Pareto front to the problem*/
  pcross = px;
  
  if (nvar > 0) 
  {
    pmut_r = pm / (double) nvar;
      
    /*Asks for distribution index for crossover (default = 20)*/
    di = dx;
      
    /*Asks for distrubution index for mutation (default = 10 to 500)*/
    dim = dm;
      
    /*Specify the limits of real-coded variables*/
    for(i = 0; i < nvar; i++)
    {
      lim_r[i][0] = vmin;
      lim_r[i][1] = vmax;
    }
    ans = vfix;
  }
  else
  {
    pmut_r = 0.0;
    di = 0.0;
    dim = 0.0;
  }
  
  if (nchrom > 0)
  {
    /*Asking for Crossover Type*/
    optype = ox;
    if (optype == 2)
      di = dx;
      
    /*Specify the no of bits for each variable
      Total sum of the bit value must be equal to chromosome length*/
		m = 8 * sizeof(unsigned long);
		if (bits > m)
		{
		  bits = m;
			fprintf(stderr,
			  "No. of bits for binary variables exceeds size of\n"
				"'unsigned long int' = %d bits. Possible overflow in decode.\n"
				"Set to %d bits.\n", m, m);
		}
    chrom = 0;
    for (i = 0; i < nchrom; i++)
    {
      vlen[i] = bits;
      chrom += bits;
    
      /*Specify the limits of binary-coded variables*/
      lim_b[i][0] = vmin;
      lim_b[i][1] = vmax;
    }

    /* Ask for mutation probablity */
    pmut_b = pm / (double) chrom;
  }
  else
  {
    chrom = 0;
    pmut_b = 0.0;
  }

  /* Allocate memory for population structure */
  memory();

  /*Set transform matrix for objective function values using default arguments*/
  transf();

  /* Get seed for random no generator and initialize */
  seed = ( (is > 0) ? (-is) : (is) );
  ran1(&seed);
  
  /* Some counters */
  nmut = 0;
  ncross = 0;
  neval = 0;

  /* pointer to class application */
  app = appl;
}

void nsga2::memory()
{
  int i;
  int *iptr;
  int **jptr;
  double *dptr;
  individual *ind_ptr;
  
  pop_ptr = (population *) malloc( sizeof(population) );
  if (pop_ptr == NULL)
  {
    fprintf(stderr,"Error allocating memory for 'pop'.\n");
    exit(1);
  }

  iptr = (int *) malloc( popsize2 * sizeof(int) );
  if (iptr == NULL)
  {
    fprintf(stderr,"Error allocating memory for 'pop.rankno'.\n");
    exit(1);
  }
  else
  {
    pop_ptr->rankno = iptr;
  }

  jptr = (int **) malloc( popsize2 * sizeof(int *) );
  if (jptr == NULL)
  {
    fprintf(stderr,"Error allocating memory for 'pop.rankar'.\n");
    exit(1);
  }
  else
  {
    pop_ptr->rankar = jptr;
  }
  iptr = (int *) malloc( popsize2 * sizeof(int) );
  if (iptr == NULL)
  {
    fprintf(stderr,"Error allocating memory for 'pop.rankar'.\n");
    exit(1);
  }
  else
  {
		for (i = 0; i < popsize2; i++)
		{
			*jptr++ = iptr;
		}
  }

  ind_ptr = (individual *) malloc( popsize2 * sizeof(individual) );
  if (ind_ptr == NULL)
  {
    fprintf(stderr,"Error allocating memory for 'pop.ind'.\n");
    exit(1);
  }
  else
  {
    pop_ptr->ind = ind_ptr;
  }
  for (i = 0; i < popsize2; i++)
  {
    if (nchrom)
    {
      iptr = (int *) malloc( chrom * sizeof(int) );
      if (iptr == NULL)
      {
        fprintf(stderr, "Error allocating memory for 'pop.ind.genes'.\n");
        exit(1);
      }
      else
      {
        ind_ptr->genes = iptr;
      }
      dptr = (double *) malloc( nchrom * sizeof(double) );
      if (dptr == NULL)
      {
        fprintf(stderr, "Error allocating memory for 'pop.ind.xbin'.\n");
        exit(1);
      }
      else
      {
        ind_ptr->xbin = dptr;
      }
    }
    if (nvar)
    {
      dptr = (double *) malloc( nvar * sizeof(double) );
      if (dptr == NULL)
      {
        fprintf(stderr, "Error allocating memory for 'pop.ind.xreal'.\n");
        exit(1);
      }
      else
      {
        ind_ptr->xreal = dptr;
      }
    }
    dptr = (double *) malloc( nfunc * sizeof(double) );
    if (dptr == NULL)
    {
      fprintf(stderr, "Error allocating memory for 'pop.ind.object'.\n");
      exit(1);
    }
    else
    {
      ind_ptr->object = dptr;
    }
    dptr = (double *) malloc( nfunc * sizeof(double) );
    if (dptr == NULL)
    {
      fprintf(stderr, "Error allocating memory for 'pop.ind.fitness'.\n");
      exit(1);
    }
    else
    {
      ind_ptr->fitness = dptr;
    }
    if (ncons)
    {
      dptr = (double *) malloc( ncons * sizeof(double) );
      if (dptr == NULL)
      {
        fprintf(stderr, "Error allocating memory for 'pop.ind.constr'.\n");
        exit(1);
      }
      else
      {
        ind_ptr->constr = dptr;
      }
    }
    ind_ptr++;
  } /*end of population loop*/

  
  /*Allocate memory for flag array*/
  flag = (int *) malloc( popsize2 * sizeof(int) );
  if (flag == NULL)
  {
    fprintf(stderr, "Error allocating memory for 'flag'.\n");
    exit(1);
  }
  /*Allocate memory for index array*/
  index = (int *) malloc( popsize2 * sizeof(int) );
  if (index == NULL)
  {
    fprintf(stderr, "Error allocating memory for 'index'.\n");
    exit(1);
  }
  /*Allocate memory for fit array*/
  fit = (double *) malloc( popsize2 * sizeof(double) );
  if (fit == NULL)
  {
    fprintf(stderr, "Error allocating memory for 'fit'.\n");
    exit(1);
  }
  /*Allocate memory for x array*/
  x = (double *) malloc( (nchrom + nvar) * sizeof(double) );
  if (x == NULL)
  {
    fprintf(stderr, "Error allocating memory for 'x'.\n");
    exit(1);
  }
  /*Allocate memory for f array*/
  f = (double *) malloc( nfunc * sizeof(double) );
  if (f == NULL)
  {
    fprintf(stderr, "Error allocating memory for 'f'.\n");
    exit(1);
  }
  /*Allocate memory for array lim_f*/
  lim_f = (double **) malloc( nfunc * sizeof(double *) );
  if (lim_f == NULL)
  {
    fprintf(stderr, "Error allocating memory for 'lim_f'.\n");
    exit(1);
  }
  for (i = 0; i < nfunc; i++)
  {
     lim_f[i] = (double *) malloc( 2 * sizeof(double) );
     if (lim_f[i] == NULL)
     {
        fprintf(stderr, "Error allocating memory for 'lim_f[%d]'.\n",i);
        exit(1);
     }
  }
  /*Allocate memory for array T*/
  T = (double **) malloc( nfunc * sizeof(double *) );
  if (T == NULL)
  {
    fprintf(stderr, "Error allocating memory for 'T'.\n");
    exit(1);
  }
  for (i = 0; i < nfunc; i++)
  {
     T[i] = (double *) malloc( nfunc * sizeof(double) );
     if (T[i] == NULL)
     {
        fprintf(stderr, "Error allocating memory for 'T[%d]'.\n",i);
        exit(1);
     }
  }
  /*Allocate memory for cstr array*/
  if (ncons)
  {
    cstr = (double *) malloc( ncons * sizeof(double) );
    if (cstr == NULL)
    {
      fprintf(stderr, "Error allocating memory for 'cstr'.\n");
      exit(1);
    }
  }
  /*Allocate memory for oldfit array*/
  oldfit = (double **) malloc( nfunc * sizeof(double *) );
  if (oldfit == NULL)
  {
    fprintf(stderr,"Error allocating memory for 'pop.oldfit'.\n");
    exit(1);
  }
  else
  {
    for (i = 0; i < nfunc; i++)
    {
      oldfit[i] = (double *) malloc( popsize * sizeof(double) );
      if (oldfit[i] == NULL)
      {
        fprintf(stderr,"Error allocating memory for 'pop.oldfit'.\n");
        exit(1);
      }
    }
  }
}

void nsga2::transf(int k, int n)
/* Calculate the transform matrix for objective functions, i.e., restrict search
   to a part of the Pareto front.
   Input values are k, k = 1, 2, ..., or n the k-th part of the subdivision 
                and n the total number of parts the Pareto-front is subdivided into.
 */
{
   double t;
   int i, j;

   switch (nfunc)
   {
   case 1:
      T[0][0] = 1.0;
      break;
   case 2:
      t = M_PI * double(n - k) / double(2 * n);
      T[0][0] = cos(t);
      T[0][1] = sin(t);
      t = M_PI * double(k - 1) / double(2 * n);
      T[1][0] = sin(t);
      T[1][1] = cos(t);
      break;
   default:
      cerr << "Splitting the Pareto front for " << nfunc 
         << " objective functions is not implemented, yet." << endl;
      for (i = 0; i < nfunc; i++)
         for (j = 0; j < nfunc; j++)
            T[i][j] = (i == j ? 1.0 : 0.0);
      break;
   }

   if (app != NULL && n > 1)
   {
      app->setrange(lim_f, nfunc);
   }
   else
   {
      for (i = 0; i < nfunc; i++)
      {
         lim_f[i][0] = 0.0;
         lim_f[i][1] = 1.0;
     }
   }
}

nsga2::~nsga2()
{
  int i;
  individual *ind_ptr;

  /* free memory allocated in memory() */
  for (i = 0; i < nfunc; i++)
  {
    free(oldfit[i]);
  }
  free(oldfit);
  if(ncons)
  {
    free(cstr);
  }
  for (i = 0; i < nfunc; i++)
  {
     free(T[i]);
  }
  free(T);
  for (i = 0; i < nfunc; i++)
  {
     free(lim_f[i]);
  }
  free(lim_f);
  free(f);
  free(x);
  free(fit);
  free(index);
  free(flag);

  ind_ptr = pop_ptr->ind + (popsize2 - 1);
  for (i = 0; i < popsize2; i++)
  {
    if (ncons)
    {
      free(ind_ptr->constr);
    }
    free(ind_ptr->fitness);
    free(ind_ptr->object);
    if (nvar)
    {
      free(ind_ptr->xreal);
    }
    if (nchrom)
    {
      free(ind_ptr->xbin);
      free(ind_ptr->genes);
    }
    ind_ptr--;
  }
  free(pop_ptr->ind);
  free(*pop_ptr->rankar);
  free(pop_ptr->rankar);
  free(pop_ptr->rankno);
  free(pop_ptr);

  /* free memory allocated in input() */
  for (i = 0; i < nvar; i++)
  {
    free(lim_r[i]);
  }
  free(lim_r);
  for (i = 0; i < nchrom; i++)
  {
    free(lim_b[i]);
  }
  free(lim_b);
  free(vlen);
}

