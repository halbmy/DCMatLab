#include "nsga2.h"

void nsga2::selection()
{
  int i,j,rnk1,rnk2;
  double len1,len2;
  individual *ind1,*ind2;
  int *bit1,*bit2;
  double *ptr1,*ptr2;

  for (i = 0; i < popsize; i++)
  {
    /* pick first candidate for tournament */
    j = (int) ((double) popsize * ran1(&seed));
    if (j == popsize)
    {
      j--;
    }
    ind1 = &(pop_ptr->ind[j]);
    rnk1 = ind1->rank;
    len1 = ind1->cub_len;
    
    /* pick second candidate for tournament */
    j = (int) ((double) popsize * ran1(&seed));
    if (j == popsize)
    {
      j--;
    }
    ind2 = &(pop_ptr->ind[j]);
    rnk2 = ind2->rank;
    len2 = ind2->cub_len;

    /* compare candidates, ind1 holds winner */
    if (rnk1 > rnk2)
    {
      ind1 = ind2;
    }
    else
    {
      if (rnk1 == rnk2 && len1 < len2)
      {
        ind1 = ind2;
      }
    }

    /* copy winner to new population */
    ind2 = &(pop_ptr->ind[popsize + i]);
    bit1 = ind1->genes;
    bit2 = ind2->genes;
    for (j = 0; j < chrom; j++)
    {
      *bit2++ = *bit1++;
    }
    ptr1 = ind1->xbin;
    ptr2 = ind2->xbin;
    for (j = 0; j < nchrom; j++)
    {
      *ptr2++ = *ptr1++;
    }
    ptr1 = ind1->xreal;
    ptr2 = ind2->xreal;
    for (j = 0; j < nvar; j++)
    {
      *ptr2++ = *ptr1++;
    }
    ptr1 = ind1->object;
    ptr2 = ind2->object;
    for (j = 0; j < nfunc; j++)
    {
      *ptr2++ = *ptr1++;
    }
    ptr1 = ind1->fitness;
    ptr2 = ind2->fitness;
    for (j = 0; j < nfunc; j++)
    {
      *ptr2++ = *ptr1++;
    }
    ptr1 = ind1->constr;
    ptr2 = ind2->constr;
    for (j = 0; j < ncons; j++)
    {
      *ptr2++ = *ptr1++;
    }
    ind2->error = ind1->error;
/*    ind2->change = ind1->change;
      Should be 0 if everything is correct! */
  }
}

