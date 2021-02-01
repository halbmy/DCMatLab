/* This is a Multi-Objective GA program.
**********************************************************************
*  This program is the implementation of the NSGA-2 proposed by      *
*                                                                    *
*  Prof. Kalyanmoy Deb and his students .                            *
*                                                                    *
*  copyright Kalyanmoy Deb                                           *
*                                                                    *
*  revised by Christoph Schwarzbach                                  *
**********************************************************************
The user have to give the input manualy or through a data file.

The user needs to enter objective functions in func-con.h
The code can also take care of the constraints. Enter the constraints
in the space provided in the func-con.h file.
Constraints must be of the following type:
g(x) >= 0.0
Also normalize all constraints (see the example problem in func-con.h)

The program generates few output files. These are described as
 1.report.out
*           This file has the detailed record for all the variables,
*           the fitness values, constraint values, overall constraint
*           violation (penalty) and their ranks for all the members
*           of the initial and subsequent populations.

 2.fitness.out
*                 This file has the fitness value of all feasible and
*                 non-dominated individuals at the final generation

 3.parameters.out
*                 This file has the all the variables of the feasible
*                 and non-dominated individuals at the final generation.
*                 The i-th solutions here corresponds to the i-th solution
*                 in the final_fitness.out file. 

Compilation procedure:  gcc -o nsga2 nsga2.c -lm
Run nsga2 -[q|v] with or without an input file.
Optional command line parameters are
 -q   suppress detailed output of all generations in report.out
      and progress on screen
 -v   (default) write complete report on report.out and print progress on
      screen

Input data files: Three files are included, but at one time one is needed
depending on the type of variables used:
inp-r (template file input-real)  : All variables are real-coded
inp-b (template file input-binary): All variables are binary-coded
inp-rb(template file input-rl+bin): Some variables are real and some are binary  
*/

#include "nsga2.h"
#include <time.h>
#include <string.h>


void nsga2::run(int verbose)
{
  int i;
  individual *ind;
  time_t t1, t2;

  verb = verbose;

  /*  */
  printf("NSGA-2 is processing ...\n");
	i = sizeof(population) + 2*popsize * sizeof(individual) +
		(nchrom + 2*popsize * (4 + chrom)) * sizeof(int) +
		2*popsize * sizeof(int *) + 
		(3*(nchrom+nvar) + nfunc + ncons + 
		 2*popsize * (1 + (nchrom+nvar) + nfunc + ncons)) * sizeof(double) + 
		(nchrom+nvar) * sizeof(double *);
	if(i > 1024)
	{
		i /= 1024;
	  if(i > 1024)
		{
	    printf(" ... requires about %0.1f MBytes of memory on heap.\n", 
				(double) i / 1024.0);
		}
		else
		{
	    printf(" ... requires about %d kBytes of memory on heap.\n", i);
		}
	}
	else
	{
	  printf(" ... requires about %d Bytes of memory on heap.\n", i);
	}

  time(&t1);

  /*Pointer to first block of individuals*/
  ind = pop_ptr->ind;

  /*Initialization of population*/
  init(ind);  

  /*Decode binary strings*/
  if (nchrom)
  {
    decode(ind); 
  }
  
  /*Function Calculation*/
  func(ind);

  /*Ranking*/
  if(ncons == 0)
  {
    rank(popsize);
  }
  else
  {
    rankc(popsize);
  }

  /*Report of initial population*/
  report(0);
  
  /*Pointer to second block of individuals*/
  ind += popsize;

  /********************************************************************/
  /*----------------------GENERATION STARTS HERE----------------------*/
  for (i = 0; i < gener; i++)
  {
		if (verb)
		{
			printf("generation no. %d\n",i+1);
		}

    /*--------SELECT----------------*/
    selection();
    
    /*CROSSOVER----------------------------*/      
    if (nchrom > 0) 
    {
      if(optype == 1) /*Binary single-point crossover*/
      {
        simplecross(ind);
      }
      if(optype == 2) /*Binary uniform crossover*/
      {
        unicross(ind);
      }
			if(optype == 3) /*Binary multi-point crossover*/
			{
				multicross(ind);
			}
    }
    if (nvar > 0)     
    {
      sbcross(ind); /*Real SBX crossover*/
    }
      
    /*------MUTATION-------------------*/
    if (nchrom > 0) /*Binary Mutation */
    {
      mutate(ind);
    }
    if (nvar > 0)   /*Real Mutation*/
    {
      realmutate(ind);
    }
      
    /*-------DECODING----------*/
    if(nchrom > 0)
    {
      decode(ind);
    }
      
    /*----------FUNCTION EVALUATION-----------*/
    func(ind);

    /*-------------------SELECTION KEEPING FRONTS ALIVE--------------*/
      
    /*Elitism and sharing implemented,
      Ranking of whole population to form new subpopulation for next 
      generation*/
    keepalive(popsize2);

    /*------------------CHECK CONVERGENCE--------------*/
    converge(i+1);

    /*------------------REPORT PRINTING--------------------------------*/
    report(i+1);

  }  /* end of i */
  /*                   Generation Loop Ends                                */
  /************************************************************************/

  time(&t2);
  printf("Finished in ");
  printtime(difftime(t2,t1));
  printf(".\n"
    "%d of %d function evaluations (%.1f%%).\n",
    neval, (gener+1)*popsize, 
    100.0 * ((double) neval) / ((double) ((gener+1)*popsize)) );
}

void nsga2::printtime(double secs)
{
  double mins, hrs, days;
  
  if(secs > 60.0)
  {
    mins = floor(secs/60.0);
    secs -= 60.0*mins;
    if(mins > 60.0)
    {
      hrs = floor(mins/60.0);
      mins -= 60.0*hrs;
      if(hrs > 24.0)
      {
        days = floor(hrs/24.0);
        hrs -= 24.0*hrs;
        printf("%.0f d, ", days);
      }
      printf("%.0f h, ", hrs);
    }
    printf("%.0f min, ", mins);
  }
  printf("%.0f s", secs);
}


/* Start evolution (initialization) */
void nsga2::start(int verbose, FILE *stream)
{
  individual *ind;

  verb = 0;

  if (verbose)
  {
    int i = sizeof(population) + 2*popsize * sizeof(individual) +
      (nchrom + 2*popsize * (4 + chrom)) * sizeof(int) +
      2*popsize * sizeof(int *) + 
      (3*(nchrom+nvar) + nfunc + ncons + 
      2*popsize * (1 + (nchrom+nvar) + nfunc + ncons)) * sizeof(double) + 
      (nchrom+nvar) * sizeof(double *);
    if(i > 1024)
    {
      i /= 1024;
      if(i > 1024)
      {
        printf("NSGA2 with population size %d "
          "requires about %0.1f MBytes of memory on heap.\n", 
          popsize, (double) i / 1024.0);
      }
      else
      {
        printf("NSGA2 with population size %d "
          "requires about %d kBytes of memory on heap.\n",
          popsize, i);
      }
    }
    else
    {
      printf("NSGA2 with population size %d "
        "requires about %d Bytes of memory on heap.\n",
        popsize, i);
    }
  }

  /*Pointer to first block of individuals*/
  ind = pop_ptr->ind;

  /*Initialization of population*/
  init(ind, stream);  

  /*Decode binary strings*/
  if (nchrom)
  {
    decode(ind); 
  }
  
  /*Function Calculation*/
  func(ind);

  /*Ranking*/
  if(ncons == 0)
  {
    rank(popsize);
  }
  else
  {
    rankc(popsize);
  }
}

/* Perform 1 generation */
void nsga2::step(int ngen)
{
  individual *ind;
  
  ind = pop_ptr->ind + popsize;
  
  while (ngen > 0)
  {
    /*--------SELECT----------------*/
    selection();
  
    /*CROSSOVER----------------------------*/      
    if (nchrom > 0) 
    {
      if(optype == 1) /*Binary single-point crossover*/
      {
        simplecross(ind);
      }
      if(optype == 2) /*Binary uniform crossover*/
      {
        unicross(ind);
      }
      if(optype == 3) /*Binary multi-point crossover*/
      {
        multicross(ind);
      }
    }
    if (nvar > 0)     
    {
      sbcross(ind); /*Real SBX crossover*/
    }
  
    /*------MUTATION-------------------*/
    if (nchrom > 0) /*Binary Mutation */
    {
      mutate(ind);
    }
    if (nvar > 0)   /*Real Mutation*/
    {
      realmutate(ind);
    }
  
    /*-------DECODING----------*/
    if(nchrom > 0)
    {
      decode(ind);
    }
  
    /*----------FUNCTION EVALUATION-----------*/
    func(ind);
  
    /*-------------------SELECTION KEEPING FRONTS ALIVE--------------*/
  
    /*Elitism and sharing implemented,
    Ranking of whole population to form new subpopulation for next 
    generation*/
    keepalive(popsize2);

    /*Keep track of convergency in conv*/
    converge();

    ngen--;
  }
}

/*
  statistics()
    Calculates some statistics to describe current population
*/
void nsga2::statistics(double *mean)
{
  int i, j, n;
  individual *ind_ptr = pop_ptr->ind;

  n = 0;
  for (j = 0; j < nfunc; j++)
  {
    mean[j] = 0.0;
  }

  for (i = 0; i < popsize; i++)
  {
    if (ind_ptr->error <= 0.0)
    {
      n++;
      for (j = 0; j < nfunc; j++)
      {
        mean[j] += ind_ptr->fitness[j];
      }
    }
  }

  if (n > 0)
  {
    for (j = 0; j < nfunc; j++)
    {
      mean[j] /= (double) n;
    }
  }
}

