/*This is the file for formulating the crossover process*/
#include "nsga2.h"

void nsga2::simplecross(individual *ind_ptr)
{
  int i,j,mating_site,gene;
  int *gene1,*gene2;
  double rnd;

  for (i = 0; i < popsize/2; i++)
  {
    rnd = ran1(&seed);
    if (rnd < pcross)
    {
      ncross++;
      rnd = ran1(&seed);
      mating_site = (int) (rnd * (double) chrom);
      if(mating_site >= chrom)
      {
        mating_site = mating_site/2;
      }
      gene1 = ind_ptr->genes + mating_site;
      ind_ptr++;
      gene2 = ind_ptr->genes + mating_site;
      ind_ptr++;
      for(j = mating_site; j < chrom; j++)
      {
        gene = *gene1;
        *gene1++ = *gene2;
        *gene2++ = gene;
      }
    }
    else 
    {
      ind_ptr += 2;
    }
  }
}

/***********************************************************************/

void nsga2::multicross(individual *ind_ptr)
{
  int i,j,k,mating_site,gene;
  int *gene1,*gene2;
  double rnd;

  for (i = 0; i < popsize/2; i++)
  {
    rnd = ran1(&seed);
    if (rnd < pcross)
    {
      ncross++;
      gene1 = ind_ptr->genes;
      ind_ptr++;
      gene2 = ind_ptr->genes;
      ind_ptr++;
			for (j = 0; j < nchrom; j++)
			{
        rnd = ran1(&seed);
        mating_site = (int) (rnd * (double) vlen[j]);
				if(mating_site >= vlen[j])
				{
          mating_site = mating_site/2;
				}
        gene1 += mating_site;
				gene2 += mating_site;
        for(k = mating_site; k < vlen[j]; k++)
				{
          gene = *gene1;
          *gene1++ = *gene2;
          *gene2++ = gene;
				}
			}
    }
    else 
    {
      ind_ptr += 2;
    }
  }
}

/***********************************************************************/

void nsga2::unicross(individual *ind_ptr)
{
  int i, j, gene, flag=0;
  int *gene1,*gene2;
  double rnd;
  individual *ptr;

  for(i = 0; i < popsize/2; i++)
  {
    rnd = ran1(&seed);
    if (rnd < pcross)
    {
      gene1 = ind_ptr->genes;
      ptr = ind_ptr++;
      gene2 = ind_ptr->genes;
      
      for(j = 0;j < chrom;j++)
      {
        rnd = ran1(&seed);
        if(rnd <= di && *gene1 != *gene2)
        {
          gene = *gene1;
          *gene1 = *gene2;
          *gene2 = gene;
          flag++;
        }
        gene1++;
        gene2++;
      }
      
      if(flag)
      {
        ptr->change = 1;
        ind_ptr->change = 1;
        ncross += flag;
        flag = 0;
      }
      ind_ptr++;
    }
    else
    {
      ind_ptr += 2;
    }
  }
}

