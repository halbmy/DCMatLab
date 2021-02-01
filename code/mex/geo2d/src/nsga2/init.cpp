/*This is the file which initializes the population*/
#include "nsga2.h"

void nsga2::init(individual *ind_ptr)
{
  int i,j;
  int *gene;
  double d;
  double *val;
  
  /*Loop Over the population size*/
  for (i = 0; i < popsize; i++)
  { 
    /*Binary variables*/
    gene = ind_ptr->genes;
    for (j = 0; j < chrom; j++)
    {
      /*Generate a Random No. if it is less than 0.5 it 
        generates a 0 in the string otherwise 1*/
      d = ran1(&seed);
      if(d >= 0.5)
      {
        *gene++ = 1;
      }
      else
      {
        *gene++ = 0;
      } 
    }

    /*Real variables*/
    val = ind_ptr->xreal;
    for (j = 0; j < nvar; j++)
    {
      d = ran1(&seed);
      /*if limits are not fixed then generate any number between 
        -infinity and infinity*/
      if(ans != 1)
      {
        *val++ = 1.0 / (2.0 * d - 1.0) ;
      }
      /*if limits are fixed generate a value in the 
        range of minimum and maximum value of the variable*/
      else
      {
        *val++ = lim_r[j][0] + d * (lim_r[j][1] - lim_r[j][0]);
      }
    }

    /* Set cublen = 0 - looks better in output for 1st generation*/
    ind_ptr->cub_len = 0.0;

    /*Evaluation flag*/
    ind_ptr->change = 1;

    ind_ptr++;
  } /* end of population loop */

  /* Don't forget evaluation flags for second half of population
     where children are copied during selection */
  while (i < popsize2)
  {
    ind_ptr->change = 0;
    ind_ptr++;
    i++;
  }

  noldfit = 0;
  conv = 0;
}

void nsga2::init(individual *ind_ptr, FILE *stream)
{
  int i = 0, j, n, items;
  int *gene;
  double d;
  double *val;

//  printf("init: ");
  if (stream != NULL)
  {
//	  printf("reading population from file ... ");
	  items = fread(&generi, sizeof(generi), 1, stream);
	  if (items != 1 || generi < 0)
	  {
		  printf(" Input error during initialization using dumpfile:\n" 
			  " Generating initial population at random.\n");
		  generi = 0;
		  init(pop_ptr->ind);
		  return;
	  }
	  items = fread(&n, sizeof(popsize), 1, stream);
	  if (items != 1)
	  {
		  printf(" Input error during initialization using dumpfile:\n" 
			  " Generating initial population at random.\n");
		  generi = 0;
		  init(pop_ptr->ind);
		  return;
	  }
	  items = fread(&i, sizeof(chrom), 1, stream);
	  if (items != 1)
	  {
		  printf(" Input error during initialization using dumpfile:\n" 
			  " Generating initial population at random.\n");
		  generi = 0;
		  init(pop_ptr->ind);
		  return;
	  }
	  items = fread(&j, sizeof(nvar), 1, stream);
	  if (items != 1)
	  {
		  printf(" Input error during initialization using dumpfile:\n" 
			  " Generating initial population at random.\n");
		  generi = 0;
		  init(pop_ptr->ind);
		  return;
	  }
	  if (i != chrom)
	  {
		  printf(" Error during initialization using dumpfile:\n"
			  " Mismatch in length of binary chromosome.\n"
			  " Generating initial population at random.\n");
		  generi = 0;
		  n = 0;
	  }
	  if (j != nvar)
	  {
		  printf(" Error during initialization using dumpfile:\n"
			  " Mismatch in length of real chromosome.\n"
			  " Generating initial population at random.\n");
		  generi = 0;
		  n = 0;
	  }
	  if (n > popsize)
	  {
		  printf(" Warning during initialization using dumpfile:\n"
			  " Number of records exceeds population size.\n"
			  " Skipping last %d records.\n", popsize-n);
		  n = popsize;
	  }
	  
	  /* Read initial values for n individuals from dumpfile */
	  for (i = 0; i < n; i++)
	  { 
		  /*Binary variables*/
		  items = fread(ind_ptr->genes, sizeof(*(ind_ptr->genes)), chrom, stream);
		  if (items != chrom)
		  {
			  printf(" Input error during initialization using dumpfile:\n" 
				  " Generating initial population at random.\n");
			  generi = 0;
			  init(pop_ptr->ind);
			  return;
		  }
		  
		  /*Real variables*/
		  items = fread(ind_ptr->xreal, sizeof(*(ind_ptr->xreal)), nvar, stream);
		  if (items != nvar)
		  {
			  printf(" Input error during initialization using dumpfile:\n" 
				  " Generating initial population at random.\n");
			  generi = 0;
			  init(pop_ptr->ind);
			  return;
		  }
		  
		  /* Set cublen = 0 - looks better in output for 1st generation*/
		  ind_ptr->cub_len = 0.0;
		  
		  /*Evaluation flag*/
		  ind_ptr->change = 1;
		  
		  ind_ptr++;
	  }
//	  printf("done.\n");
  }

  /* Generate remaining individuals at random*/
  for (; i < popsize; i++)
  { 
    /*Binary variables*/
    gene = ind_ptr->genes;
    for (j = 0; j < chrom; j++)
    {
      /*Generate a Random No. if it is less than 0.5 it 
        generates a 0 in the string otherwise 1*/
      d = ran1(&seed);
      if(d >= 0.5)
      {
        *gene++ = 1;
      }
      else
      {
        *gene++ = 0;
      } 
    }

    /*Real variables*/
    val = ind_ptr->xreal;
    for (j = 0; j < nvar; j++)
    {
      d = ran1(&seed);
      /*if limits are not fixed then generate any number between 
        -infinity and infinity*/
      if(ans != 1)
      {
        *val++ = 1.0 / (2.0 * d - 1.0) ;
      }
      /*if limits are fixed generate a value in the 
        range of minimum and maximum value of the variable*/
      else
      {
        *val++ = lim_r[j][0] + d * (lim_r[j][1] - lim_r[j][0]);
      }
    }

    /* Set cublen = 0 - looks better in output for 1st generation*/
    ind_ptr->cub_len = 0.0;

    /*Evaluation flag*/
    ind_ptr->change = 1;

    ind_ptr++;
  } /* end of population loop */


  /* Don't forget evaluation flags for second half of population
     where children are copied during selection */
  while (i < popsize2)
  {
    ind_ptr->change = 0;
    ind_ptr++;
    i++;
  }

  noldfit = 0;
  conv = 0;
}

