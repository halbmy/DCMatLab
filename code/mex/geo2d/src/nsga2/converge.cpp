#include "nsga2.h"

void nsga2::converge(int gen)
{
  static FILE *stream;

  if(gen == 1)
  {
    stream = fopen("convergence.out", "w");
    if (stream == NULL)
    {
      fprintf(stderr, "Error opening file 'convergence.out'.\n");
      exit(1);
    }
    fprintf(stream, "# Generation\tno. of solutions in pareto front\n");
  }

  fprintf(stream, "%d\t%d\n", gen, pop_ptr->rankno[0]);

  if(gen == gener)
  {
    fclose(stream);    
	}
}

/* Return measure for rate of convergence: number of generations where no
   individual in the child population is found which dominates the parent
   population.
   class variables:
     int noldfit, conv;
     double **oldfit; */
void nsga2::converge()
{
  int i, j, k, n, dom;
  individual *ind_ptr;

  n = pop_ptr->rankno[0];
  k = nfunc;
  if (noldfit > 0)
  {
    dom = 0;
    for (i = 0; i < n; i++)
    {
      ind_ptr = pop_ptr->ind + pop_ptr->rankar[0][i];
      for (j = 0; j < noldfit; j++)
      {
        for (k = 0; k < nfunc; k++)
        {
          if (oldfit[k][j] > ind_ptr->fitness[k])
          {
            break;
          }
        }
        if (k < nfunc) 
        {
          break;
        }
      }
      if (k < nfunc)
      {
        dom++;
      }
    }
  }
  else
  {
    dom = n;
  }
  noldfit = n;
  for (i = 0; i < noldfit; i++)
  {
    ind_ptr = pop_ptr->ind + pop_ptr->rankar[0][i];
    for (j = 0; j < nfunc; j++)
    {
      oldfit[j][i] = ind_ptr->fitness[j];
    }
  }
  if (dom == 0)
  {
    conv++;
  }
  else
  {
    if (conv > 0)
    {
      conv--;
    }
    else
    {
      conv = 0;
    }
  }
}

