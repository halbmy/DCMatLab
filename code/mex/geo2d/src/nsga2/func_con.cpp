/********************************************************************
*                    Objective functions                            *
********************************************************************/
#include "nsga2.h"

void nsga2::func(individual *ind_ptr)
{
  int i,j,k; 
  double t1, t2;
  double *ptr1, *ptr2;

  for(i = 0;i < popsize;i++)
  {
    if (ind_ptr->change)
    {
      ind_ptr->change = 0; /* Unset evaluation flag */
      neval++;

      ptr1 = ind_ptr->xreal; /* Real-coded variables */
      for(j = 0; j < nvar; j++)
      {
        x[j] = *ptr1++;
      }
      ptr1 = ind_ptr->xbin; /* Binary-coded variables */
      k = nvar;
      for(j = 0; j < nchrom; j++)
      {
        x[k++] = *ptr1++;
      }

  	  /*
	     * Evaluate fitness functions and constraints with arguments
       * x[0...(nvar+nchrom-1)] and put results in vectors
       * f[0...(nfunc-1)] and cstr[0...(ncons-1)] 
	     */
	   app->eval(x, nvar+nchrom, f, nfunc, cstr, ncons);

      /* Copy temporary arrays to individuals*/
      ptr1 = ind_ptr->object;
      ptr2 = ind_ptr->fitness;
      for(j = 0; j < nfunc; j++)
      {
         ptr1[j] = f[j];
         ptr2[j] = 0.0;
         for (k = 0; k < nfunc; k++)
            ptr2[j] += T[j][k] * (f[k] - lim_f[k][0]) / (lim_f[k][1] - lim_f[k][0]);
      }
      ptr1 = ind_ptr->constr;
      t2  = 0.0;
      for (j = 0; j < ncons; j++)
      {
        t1 = cstr[j];
        *ptr1++ = t1;
        if (t1 < 0.0) t2 -= t1;
      }
      ind_ptr->error = t2;
    }
    ind_ptr++;
  }
}

