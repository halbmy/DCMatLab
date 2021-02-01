/* This is the module used to formulate the mutation routine*/
#include "nsga2.h"

void nsga2::mutate(individual *ind_ptr)
{
  int i,j,flag=0;
  int *gene;
  double rnd;

  for(j = 0; j < popsize; j++)
  {
    gene = ind_ptr->genes;
    for (i = 0; i < chrom; i++)
    {
      rnd = ran1(&seed); /*Check whether to do mutation or not*/
      if(rnd <= pmut_b)
      {
        *gene = !(*gene);
        flag++;
      }
      gene++;
    }
    if (flag)
    {
      ind_ptr->change = 1;
      nmut += flag;
      flag = 0;
    }
    ind_ptr++;     
  }
}

