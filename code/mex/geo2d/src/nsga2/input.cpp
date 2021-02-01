/*This is a file to get the input for the GA program*/
#include "nsga2.h"

void nsga2::input()
{
  int i, l, m;
  
  if (verb)
  printf("                                                                  \n"
         "                           NSGA-2                                 \n"
         "           by Prof. Kalyanmoy Deb and his students                \n"
         "               revised by Christoph Schwarzbach                   \n"
         "------------------------------------------------------------------\n"
         "This is a multi-objective GA program to solve constraint problems.\n"
         "------------------------------------------------------------------\n"
         "                                                                  \n"
         "1) Problem specification\n"
         "------------------------\n");
  
  /*Asks for number of the variables*/
  if (verb)
	printf("Give no. of real and binary-coded variables\n"
         " > ");
  scanf("%d %d", &nvar, &nchrom);

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
  if (verb)
  printf("Give no. of objective functions\n"
         " > ");
  scanf("%d", &nfunc);

  /*Asks for number of constraints*/
  if (verb)
  printf("Give no. of constraints\n"
         " > ");
  scanf("%d", &ncons);
  
  if (verb)
  printf("\n"
         "2) GA parameters\n"
         "----------------\n");
  
  /*Asks for number of the individuals in the population*/
  /*Default value = 30 to 100*/
  if (verb)
  printf("Give population size (an even no.)\n"
         " > ");
  scanf("%d", &popsize);
  popsize2 = 2 * popsize;
  
  /*No. of generations for which the GA will let the population evolve
    Default value is 100
    Too large value will take very long time and very small value will
    not let the GA reach the global Pareto front to the problem*/
  if (verb)
  printf("Give no. of generations\n"
         " > ");
  scanf("%d", &gener);

	/*Ask for control parameter for elitism
	  permissible values are 0 (uncontrolled elitsm)
	  to <1 (no elitism)*/
	if (verb)
	printf("Give control parameter for elitism (0.0 to 1.0)\n"
	       " > ");
	scanf("%lf", &relite);
	if (relite < 0.0) relite = 0.0;
	if (relite >= 1.0) relite = 0.999;

	/*Ask for length of sequences for sharing in fitness and parameter
	  space, respectively: nsharefit and nsharepar generations will be
	  shared with respect to distances in fitness and parameter space
	  switching between both*/
	if (verb)
	printf("Give length of sequence to perform sharing in fitness and "
	       "parameter space\n"
				 " > ");
	scanf("%d %d", &nsharefit, &nsharepar);
	if (nsharefit + nsharepar > gener)
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
  if (verb)
  printf("Give crossover probability (0.5 to 1.0)\n"
         " > ");
  scanf("%lf", &pcross);
  
  if (nvar > 0) 
  {
    if (verb)
    printf("Give mutation probability for real-coded vectors"
      " (between 0 and %f)\n"
      " > ", 1.0/( (double) nvar ));
    scanf("%lf", &pmut_r);
      
    /*Asks for distribution index for crossover (default = 20)*/
    if (verb)
    printf("Give distribution index for real-coded crossover"
      " (0.5 to 100)\n"
      " > ");
    scanf("%lf", &di);
      
    /*Asks for distrubution index for mutation (default = 10 to 500)*/
    if (verb)
    printf("Give distribution index for real-coded mutation"
      " (0.5 to 500)\n"
      " > ");
    scanf("%lf", &dim);
      
    /*Specify the limits of real-coded variables*/
    for(i = 0; i < nvar; i++)
    {
      if (verb)
      printf("Give lower & upper limits of real-coded variable no. %d\n"
        " > ", i+1);
      scanf("%lf %lf",&lim_r[i][0],&lim_r[i][1]);
    }
    if (verb)
    printf("Are the limits for real-coded variables rigid?\n"
      " 0 for no, 1 for yes\n"
      " > ");
      scanf("%d", &ans);
  }
  
  if (nchrom > 0)
  {
    /*Asking for Crossover Type*/
    if (verb)
    printf("Give binary crossover type\n"
           " 1 for simple crossover\n"
           " 2 for uniform crossover\n"
					 " 3 for multipoint crossover\n"
           " > ");
    scanf("%d", &optype);
    /*If uniform Crossover, ask for its parameter*/
    if (optype == 2)
    {
      if (verb)
        printf("Give parameter for uniform crossover (0...0.5)\n"
               " > ");
      scanf("%lf", &di);
    }
      
    /*Specify the no of bits for each variable
      Total sum of the bit value must be equal to chromosome length*/
    chrom = 0;
		m = 8 * sizeof(unsigned long);
    for (i = 0; i < nchrom; i++)
    {
      if (verb)
      printf("Give no. of bits assigned to binary-coded variable no. %d\n"
        " > ", i+1);
      scanf("%d", &l);
			if (l > m)
			{
				l = m;
				fprintf(stderr,
					"No. of bits for binary variable no. %d exceeds size of\n"
					"'unsigned long int' = %d bits. Possible overflow in decode.\n"
					"Set to %d bits.\n",
					i+1, m, m);
			}
      vlen[i] = l;
      chrom += l;
    
      /*Specify the limits of binary-coded variables*/
      if (verb)
      printf("Give lower & upper limits of binary-coded variable no. %d\n"
        " > ", i+1);
      scanf("%lf %lf",&lim_b[i][0],&lim_b[i][1]);
    }

    /* Ask for mutation probablity */
    if (verb)
    printf("Give mutation probability for binary strings"
           " (0.0 to 1.0)\n"
           " > ");
    scanf("%lf", &pmut_b);
		pmut_b /= (double) chrom;
  }

  if (verb)
  printf("\n"
         "3) Random number generator\n"
         "--------------------------\n");
  /*Give the initial seed*/
  if (verb)
  printf("Give random seed (0 to %ld)\n"
    " > ", (long) IM);
  scanf("%ld", &seed);

  if (verb)
  printf("\n");
}

