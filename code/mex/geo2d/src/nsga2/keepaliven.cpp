/*This is the routine to keep the Pareto-fronts alive 
  controlled elitsm added by C. Schwarzbach 19/08/02*/
#include "nsga2.h"

void nsga2::keepalive(int size)
{
  int i,j,k,m,n,p,q,pool,res;
  int *i1,*i2;
  double tmp;
  double *d1,*d2;
  individual *ind1_ptr,*ind2_ptr;

  /*Finding the global ranks */
  if(ncons == 0)
  {
    rank(size);
  }
  else
  {
    rankc(size);
  }

  m = pop_ptr->maxrank;

  /* Sharing the fitness to get the dummy fitness */
  if (nshare > nsharepar)
  {
    share(size);
    nshare--;
    if(nshare == 0) nshare = nsharefit; /*nsharepar == 0*/
  }
  else
  {
    share2(size);
    nshare--;
    if(nshare == 0) nshare = nsharefit + nsharepar;
  }
  
  /*Sort indices of rankar according to descending crowding distance cub_len*/
  for(i = 0; i < m; i++)
  {
     i1 = pop_ptr->rankar[i];
     k = pop_ptr->rankno[i];
     for(j = 0; j < k; j++)
     {
        fit[j] = -pop_ptr->ind[i1[j]].cub_len; /* '-' since heapsort */
     }                                         /* sorts to ascending order */
     if (k > 1)
     {
        heapsort2(k, fit, i1);
        if (fit[0] == fit[1] &&               /* prefer minimum fitness[0] */
           pop_ptr->ind[i1[0]].fitness[0] > pop_ptr->ind[i1[1]].fitness[0])
        {
           j = i1[1];
           i1[1] = i1[0];
           i1[0] = j;
        }
     }
  }

  /* Initializing the flags of population to zero */
  for(i = 0; i < size ;i++)
  {
    flag[i] = 0;
    index[i] = 0;
  }
  /* Decide which solutions belong to new population,
     Controlled elitism applied here */
  pool = popsize;
  res = 0;
  /* First, take predefined number of individuals from different
     Pareto-fronts */
  tmp = popsize * (1.0 - relite)/(1.0 - pow(relite, (double) m));
  for(i = 0; i < m; i++)
  {
    i1 = pop_ptr->rankar[i];
    p = pop_ptr->rankno[i];
    /* take numbers according to geometric distribution */
    q = iround(tmp) + res;
    tmp *= relite;
    pool -= q;
    /* Check, if q slots are still free (differences
       may occur due to rounding errors when calculating
       int-type q) */
    if (pool < 0)
    {
      q += pool;
      pool = 0;
    }
    /* Break if predefined numbers approach zero or pool is filled */
    if (q == 0)
    {
      break;
    }
    /* if there are more individuals in current rank, take first q
       with largest crowding distances */
    if (p > q)
    {
      res = 0;
      for(k = 0; k < q; k++)
      {
        flag[*i1++] = 1;
      }
      index[i] = q; /* keep record of number of individuals per rank */
    }
    /* else, take them all and calculate number of slots which
       remain unfilled*/
    else
    {
      res = q - p;
      pool += res; /* correct poolsize counter */
      for(k = 0; k < p; k++)
      {
        flag[*i1++] = 1;
      }
      index[i] = p; /* keep record of number of individuals per rank */
    }
  }
  /* Second, fill remaining slots with individuals left over from
     previous sweep */
  for(j = 0; j < m && pool > 0; j++)
  {
    p = pop_ptr->rankno[j];
		q = p - index[j];
		if(q > 0)
		{
			i1 = pop_ptr->rankar[j] + index[j];
			if(q > pool) q = pool;
			for (k = 0; k < q; k++)
			{
        flag[*i1++] = 1;
      }
			pool -= q;
  		index[j] += q; /* keep record of number of individuals per rank */
    }
  }

  /*Copy selected individuals to first half-part of population*/
  pop_ptr->maxrank = 0;
  for (i = 0; i < size; i++) pop_ptr->rankno[i] = 0;
  n = 0; /* counter for all individuals copied */
  ind1_ptr = pop_ptr->ind;
  ind2_ptr = pop_ptr->ind;
  for(i = 0; i < size; i++)
  {
    if(flag[i] == 1)
    {
      if(nchrom > 0)
      {
        i1 = ind1_ptr->genes;

        i2 = ind2_ptr->genes;
        for(j = 0; j < chrom; j++)
        { 
          *i1++ = *i2++;
        }
        d1 = ind1_ptr->xbin;
        d2 = ind2_ptr->xbin;
        for(j = 0; j < nchrom; j++)
        {
          *d1++ = *d2++;
        }
      }
      if (nvar > 0)
      {
        d1 = ind1_ptr->xreal;
        d2 = ind2_ptr->xreal;
        for(j = 0; j < nvar; j++)
        {
          *d1++ = *d2++;
        } 
      }
      d1 = ind1_ptr->object;
      d2 = ind2_ptr->object;
      for(j = 0;j < nfunc;j++)
      {
        *d1++ = *d2++;
      }
      d1 = ind1_ptr->fitness;
      d2 = ind2_ptr->fitness;
      for(j = 0;j < nfunc;j++)
      {
        *d1++ = *d2++;
      }
      ind1_ptr->cub_len = ind2_ptr->cub_len;
      if(ncons)
      {
        ind1_ptr->error = ind2_ptr->error;
      }
      d1 = ind1_ptr->constr;
      d2 = ind2_ptr->constr;
      for (j = 0; j < ncons; j++)
      {
        *d1++ = *d2++;
      }
	  k = ind2_ptr->rank;
      ind1_ptr->rank = k;
/*    ind1_ptr->change = ind2_ptr->change;
      all individuals should be evaluated by now! */
      /* Update rankar, rankno, and maxrank for new generation */
      pop_ptr->maxrank = max(k, pop_ptr->maxrank);
	  k--; /* indices are zero-based, ranks not */
	  pop_ptr->rankar[k][pop_ptr->rankno[k]++] = n;
      n++;
      ind1_ptr++;
    }
    ind2_ptr++;
  } /*end of population loop for copying*/
}

/* sharing in fitness space */
void nsga2::share(int size)
{
   int i,j,k,m,n,p,rnk,maxrnk;
   int *iptr;
   double maxfit,val;
   double *ptr1,*ptr2;
   individual *ind_ptr;
   
   maxrnk = pop_ptr->maxrank;
   ind_ptr = pop_ptr->ind;
   for(i = 0; i < size; i++)
   {
      ind_ptr->cub_len = 0.0;
      ind_ptr++;
   }
   
   for(rnk = 0; rnk < maxrnk; rnk++)
   {
      iptr = pop_ptr->rankar[rnk];
      m = pop_ptr->rankno[rnk];
      n = m - 1;
      if (m > 2) /* make sure that there are at least three individuals
                 in current rank. Otherwise, cub_len has already been
                 set to appropriate value of 1.0 before */
      {
         for(i = 0; i < nfunc; i++)
         {
            for(j = 0; j < m; j++)
            {
               k = iptr[j];
               fit[j] = pop_ptr->ind[k].fitness[i];
            }
            
            /*Sort the arrays in ascending order of the fitness*/
            heapsort2(m,fit,iptr);
            
            /* Identify duplicate fitness values and vectors */
            flag[0] = 0;
            maxfit = fit[0];
            k = iptr[0];
            ptr1 = pop_ptr->ind[k].fitness;
            p = 1;
            for (j = 1; j < m; j++)
            {
               k = iptr[j];
               if (fit[j] == maxfit)
               {
                  flag[j] = p; // same fitness vector as fit[k]
                  ptr2 = pop_ptr->ind[k].fitness;
                  for (k = 0; k < nfunc; k++)
                  {
                     if (ptr1[k] != ptr2[k])
                     {
                        flag[j] = -p; // only partially identical fitness vectors
                        break;
                     }
                  }
                  flag[j-1] = flag[j];
               }
               else
               {
                  flag[j] = 0; // fitness value is different from previous item
                  maxfit = fit[j];
                  ptr1 = pop_ptr->ind[k].fitness;
                  p++;
               }
            }
            
            /* Assign crowding distances for inner points of
            current Pareto-front; extremal points are preserved
            by not diminuishing their crowding distance of 1.0 */
            
            maxfit = fit[n] - fit[0];
            for (j = 0; j < m; )
            {
               if (flag[j] == 0)
               {
                  k = iptr[j];
                  val = ( (j == 0 || j == n) ? (2.0) : ((fit[j+1] - fit[j-1]) / maxfit) );
                  pop_ptr->ind[k].cub_len += val;
                  j++;
               }
               else
               {
                  for (p = j + 1; p < m; p++)
                  {
                     if (flag[j] != flag[p])
                        break;
                  }
                  val = ( (j == 0 || p == m) ? (2.0) : ((fit[p] - fit[j-1]) / maxfit) );
                  if (flag[j] > 0)
                  {
                     if (i == 0)
                     {
                        k = iptr[j];
                        pop_ptr->ind[k].cub_len = val;
                        for (j++; j < p; j++)
                        {
                           k = iptr[j];
                           pop_ptr->ind[k].cub_len = val - 3.0;
                        }
                     }
                     else
                     {
                        for (; j < p; j++)
                        {
                           k = iptr[j];
                           if (pop_ptr->ind[k].cub_len > 0.0)
                              pop_ptr->ind[k].cub_len += val;
                           else
                              pop_ptr->ind[k].cub_len -= (2.0 - val);
                        }
                     }
                  }
                  else
                  {
                     for (; j < p; j++)
                     {
                        k = iptr[j];
                        pop_ptr->ind[k].cub_len += val;
                     }
                  }
               }
            }
         }
      }
   }
}



/* sharing in parameter space */
void nsga2::share2(int size)
{
  int i,j,k,m,n,rnk,maxrnk;
  double d1,d2;
  double *ptr1, *ptr2;
  individual *ind_ptr1, *ind_ptr2;

  maxrnk = pop_ptr->maxrank;
  ind_ptr1 = pop_ptr->ind;
  for(i = 0; i < size; i++)
  {
    ind_ptr1->cub_len = 0.0;
    ind_ptr1++;
  }

  for(rnk = 0; rnk < maxrnk; rnk++)
  {
    m = pop_ptr->rankno[rnk];
    n = m - 1;

    for (i = 0; i < n; i++)
    {
      k = pop_ptr->rankar[rnk][i];
      ind_ptr1 = pop_ptr->ind + k;
      for (j = i+1; j < m; j++)
      {
        k = pop_ptr->rankar[rnk][j];
        ind_ptr2 = pop_ptr->ind + k;
        d2 = 0.0;
        ptr1 = ind_ptr1->xreal;
        ptr2 = ind_ptr2->xreal;
        for (k = 0; k < nvar; k++)
        {
          d1 = *ptr1 - *ptr2;
          d2 += d1*d1;
          ptr1++;
          ptr2++;
        }
        ptr1 = ind_ptr1->xbin;
        ptr2 = ind_ptr2->xbin;
        for (k = 0; k < nchrom; k++)
        {
          d1 = *ptr1 - *ptr2;
          d2 += d1*d1;
          ptr1++;
          ptr2++;
        }
        ind_ptr1->cub_len += d2;
        ind_ptr2->cub_len += d2;
      }
    }
  }
  ind_ptr1 = pop_ptr->ind;
  for (i = 0; i < size; i++)
  {
    d1 = ind_ptr1->cub_len;
    ind_ptr1->cub_len = sqrt(d1);
    ind_ptr1++;
  }
}

int nsga2::iround(double x)
{
  if(x > 0.0)
    if(x - floor(x) < 0.5)
      return((int) x);
    else
      return(((int) x) + 1);
  else
    if(x - ceil(x) > -0.5)
      return((int) x);
    else
      return(((int) x) - 1);
}

void nsga2::heapsort(unsigned long n, double *ra)
/* Sorts an array ra[0..n-1] into ascending numerical order using the
   Heapsort algorithm. n is input; ra is replaced on output by its 
	 sorted rearrangement.*/
{
  unsigned long i,ir,j,l;
  double rra;
  if (n < 2) return;
  l = n >> 1;
  ir = n - 1;
  /* The index l will be decremented from its initial value down to 0
	   during the "hiring" (heap creation) phase. Once it reaches 0,
		 the index ir will be decremented from its initial value down to 0
		 during the "retirement-and-promotion" (heap selection) phase. */
  for (;;)
	{
    if (l > 0) /* Still in hiring phase. */
		{
      rra = ra[--l];
		}
		else       /* In retirement-and-promotion phase. */
		{
      rra = ra[ir];   /* Clear a space at end of array. */
      ra[ir] = ra[0]; /* Retire the top of the heap into it. */
      if (--ir == 0)  /* Done with the last promotion. */
      {
				ra[0] = rra;  /* The least competent worker of all! */
        break;
			}
		}
    i = l;         /* Whether in the hiring phase or promotion phase, we */
    j = l + l + 1; /* here set up to sift down element rra to its proper */
                   /* level. */
    
    while (j <= ir)
		{
      if (j < ir && ra[j] < ra[j+1]) 
			{
			  j++; /*Compare to the better underling. */
			}
      if (rra < ra[j]) /* Demote rra. */
			{
        ra[i] = ra[j];
        i = j;
        j = j + j + 1;
			}
			else
			{
				break; /* Found rra's level. Terminate the sift-down. */
			}
		}
    ra[i] = rra; /* Put rra into its slot. */
	}
}

void nsga2::heapsort2(unsigned long n, double *ra, int *ia)
/* Sorts an arrays ra[0..n-1] and ia[0..n-1] into ascending numerical
   order with respect to array ra using the Heapsort algorithm.
	 n is input; ra and indx are replaced on output by their sorted 
	 rearrangement.*/
{
  unsigned long i,ir,j,l,iia;
  double rra;
  if (n < 2) return;
  l = n >> 1;
  ir = n - 1;
  /* The index l will be decremented from its initial value down to 0
	   during the "hiring" (heap creation) phase. Once it reaches 0,
		 the index ir will be decremented from its initial value down to 0
		 during the "retirement-and-promotion" (heap selection) phase. */
  for (;;)
	{
    if (l > 0) /* Still in hiring phase. */
		{
      rra = ra[--l];
			iia = ia[l];
		}
		else       /* In retirement-and-promotion phase. */
		{
      rra = ra[ir];   /* Clear a space at end of array. */
			iia = ia[ir];
      ra[ir] = ra[0]; /* Retire the top of the heap into it. */
			ia[ir] = ia[0];
      if (--ir == 0)  /* Done with the last promotion. */
      {
				ra[0] = rra;  /* The least competent worker of all! */
				ia[0] = iia;
        break;
			}
		}
    i = l;         /* Whether in the hiring phase or promotion phase, we */
    j = l + l + 1; /* here set up to sift down element rra to its proper */
                   /* level. */
    
    while (j <= ir)
		{
      if (j < ir && ra[j] < ra[j+1]) 
			{
			  j++; /*Compare to the better underling. */
			}
      if (rra < ra[j]) /* Demote rra. */
			{
        ra[i] = ra[j];
				ia[i] = ia[j];
        i = j;
        j = j + j + 1;
			}
			else
			{
				break; /* Found rra's level. Terminate the sift-down. */
			}
		}
    ra[i] = rra; /* Put rra into its slot. */
		ia[i] = iia;
	}
}

