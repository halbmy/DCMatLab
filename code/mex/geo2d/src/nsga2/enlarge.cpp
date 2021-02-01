#include "nsga2.h"

void nsga2::enlarge(int newpopsize)
/*
 * Enlarge population to newpopsize and copy old popsize2 individuals
 * to new population; remaining (newpopsize2 - popsize2) slots remain
 * uninitialized
 */
{
  int i, newpopsize2 = 2*newpopsize;
  int *iptr;
  int **jptr;
  double *dptr;
  individual *ind_ptr, *ind_ptr1, *ind_ptr2;
  
  /* rankno: no. of individuals per rank */
  iptr = (int *) malloc( newpopsize2 * sizeof(int) );
  if (iptr == NULL)
  {
    fprintf(stderr,"Error allocating memory for 'pop.rankno'.\n");
    exit(1);
  }
  free(pop_ptr->rankno);
  pop_ptr->rankno = iptr;

  /* rankar: array of pointers to array containing locations
     of individuals */
  jptr = (int **) malloc( newpopsize2 * sizeof(int *) );
  iptr = (int *) malloc( newpopsize2 * sizeof(int) );
  if (jptr == NULL || iptr == NULL)
  {
    fprintf(stderr,"Error allocating memory for 'pop.rankar'.\n");
    exit(1);
  }
  free(*pop_ptr->rankar);
  free(pop_ptr->rankar);
  pop_ptr->rankar = jptr;
  for (i = 0; i < newpopsize2; i++) *jptr++ = iptr;

  /* individuals */
  ind_ptr = (individual *) malloc( newpopsize2 * sizeof(individual) );
  if (ind_ptr == NULL)
  {
    fprintf(stderr,"Error allocating memory for 'pop.ind'.\n");
    exit(1);
  }
  ind_ptr1 = ind_ptr;
  ind_ptr2 = pop_ptr->ind;
  /* keep old individuals ... */
  for (i = 0; i < popsize2; i++)
  {
    ind_ptr1->genes = ind_ptr2->genes;
    ind_ptr1->xbin = ind_ptr2->xbin;
    ind_ptr1->xreal = ind_ptr2->xreal;
    ind_ptr1->object = ind_ptr2->object;
    ind_ptr1->fitness = ind_ptr2->fitness;
    ind_ptr1->constr = ind_ptr2->constr;
    ind_ptr1->error = ind_ptr2->error;
    ind_ptr1->change = ind_ptr2->change;
    ind_ptr1++; ind_ptr2++;
  }
  free(pop_ptr->ind);
  pop_ptr->ind = ind_ptr;
  /* ... and allocate memory for additional individuals */
  for (; i < newpopsize2; i++)
  {
    if (nchrom)
    {
      iptr = (int *) malloc( chrom * sizeof(int) );
      if (iptr == NULL)
      {
        fprintf(stderr, "Error allocating memory for 'pop.ind.genes'.\n");
        exit(1);
      }
      ind_ptr1->genes = iptr;

      dptr = (double *) malloc( nchrom * sizeof(double) );
      if (dptr == NULL)
      {
        fprintf(stderr, "Error allocating memory for 'pop.ind.xbin'.\n");
        exit(1);
      }
      ind_ptr1->xbin = dptr;
    }
    if (nvar)
    {
      dptr = (double *) malloc( nvar * sizeof(double) );
      if (dptr == NULL)
      {
        fprintf(stderr, "Error allocating memory for 'pop.ind.xreal'.\n");
        exit(1);
      }
      ind_ptr1->xreal = dptr;
    }
    dptr = (double *) malloc( nfunc * sizeof(double) );
    if (dptr == NULL)
    {
      fprintf(stderr, "Error allocating memory for 'pop.ind.fitness'.\n");
      exit(1);
    }
    ind_ptr1->object = dptr;
    dptr = (double *) malloc( nfunc * sizeof(double) );
    if (dptr == NULL)
    {
      fprintf(stderr, "Error allocating memory for 'pop.ind.fitness'.\n");
      exit(1);
    }
    ind_ptr1->fitness = dptr;
    if (ncons)
    {
      dptr = (double *) malloc( ncons * sizeof(double) );
      if (dptr == NULL)
      {
        fprintf(stderr, "Error allocating memory for 'pop.ind.constr'.\n");
        exit(1);
      }
      ind_ptr1->constr = dptr;
    }
    ind_ptr1++;
  }

  
  /*Allocate memory for flag array*/
  free(flag);
  flag = (int *) malloc( newpopsize2 * sizeof(int) );
  if (flag == NULL)
  {
    fprintf(stderr, "Error allocating memory for 'flag'.\n");
    exit(1);
  }
  /*Allocate memory for index array*/
  free(index);
  index = (int *) malloc( newpopsize2 * sizeof(int) );
  if (index == NULL)
  {
    fprintf(stderr, "Error allocating memory for 'index'.\n");
    exit(1);
  }
  /*Allocate memory for fit array*/
  free(fit);
  fit = (double *) malloc( newpopsize2 * sizeof(double) );
  if (fit == NULL)
  {
    fprintf(stderr, "Error allocating memory for 'fit'.\n");
    exit(1);
  }
  /*Allocate memory for oldfit array*/
  for (i = 0; i < nfunc; i++)
  {
    free(oldfit[i]);
    oldfit[i] = (double *) malloc( newpopsize * sizeof(double) );
    if (oldfit[i] == NULL)
    {
      fprintf(stderr,"Error allocating memory for 'pop.oldfit'.\n");
      exit(1);
    }
  }
  noldfit = 0;
  conv = 0;

  popsize = newpopsize;
  popsize2 = newpopsize2;
}

