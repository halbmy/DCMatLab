#include "nsga2.h"

void nsga2::rank(int size)
{
  int i,j,k,rnk,val,q;
  double *ptr1,*ptr2;

  rnk = 0;
  for(i = 0; i < size; i++)
  { 
    flag[i] = 2;
  }

  for(k = 0; k < size; k++)
  {
    q = 0;
    rnk++;
    for(j = 0; j < size; j++)
    {
      if(flag[j] == 0)
      {
        flag[j] = 2;
      }
    }
    
    for(i = 0; i < size; i++)
    {
      if(flag[i] != 1 && flag[i] != 0) 
      {
        ptr1 = &(pop_ptr->ind[i].fitness[0]);
        for(j = 0;j < size ; j++)
        {
          if( i!= j)
          {
            if(flag[j] != 1)
            {
              ptr2 = &(pop_ptr->ind[j].fitness[0]);
              val = indcmp(ptr1,ptr2);
              if(val == 2)
              { /* individual 1 is dominated */
                flag[i] = 0;
                break;
              }
              if(val == 1)
              { /* individual 2 is dominated */
                flag[j] = 0;
              }
              if(val == 3)
              { /* individual 1 & 2 are non dominated */
                if(flag[j] != 0)
                {
                  flag[j] = 3;
                }
              }
            }
          }
        }
        if(j == size)
        {
          pop_ptr->ind[i].rank = rnk;
          flag[i] = 1;
          pop_ptr->rankar[rnk-1][q] = i;
          q++;
        }
      }
    }
    pop_ptr->rankno[rnk-1] = q;

    for(j = 0; j < size; j++)
    {
      if (flag[j] != 1) break;
    }
    if(j == size) break;

		pop_ptr->rankar[rnk] = pop_ptr->rankar[rnk-1] + q;
  } 
  pop_ptr->maxrank = rnk;
}

void nsga2::rankc(int size)
{
  int i,j,k,rnk,val,q;
  double *ptr1,*ptr2;
  double err1,err2;

  rnk = 0;
  for(i=0; i<size; i++)
  { 
    flag[i] = 2;
  }

  for(k = 0; k < size; k++)
  {
    q = 0;
    rnk++;
    for(j = 0; j < size; j++)
    {
      if(flag[j] == 0)
      {
        flag[j] = 2;
      }
    }
    
    for(i = 0; i < size; i++)
    {
      if(flag[i] != 1 && flag[i] != 0) 
      {
        ptr1 = &(pop_ptr->ind[i].fitness[0]);
        err1 = pop_ptr->ind[i].error;
        for(j = 0;j < size ; j++)
        {
          if(i != j)
          {
            if(flag[j] != 1)
            {
              ptr2 = &(pop_ptr->ind[j].fitness[0]);
              err2 = pop_ptr->ind[j].error;

              if(err1 < 1.0e-6 && err2 > 1.0e-6)
              { /* first individual is feasible, second infeasible */
                flag[j] = 0;
              }
              else
              {
                if(err1 > 1.0e-6 && err2 < 1.0e-6)
                { /* first individual is infeasible, second feasible */
                  flag[i] = 0;
                  break;
                }
                else
                { /* both individuals are feasible or both infeasible */
                  if(err1 > err2)
                  { /* first individual is more infeasible */
                    flag[i] = 0;
                    break;
                  }
                  else
                  {
                    if(err1 < err2)
                    { /* second individual is more infeasible */
                      flag[j] = 0;
                    }
                    else
                    {
                      val = indcmp(ptr1,ptr2);
                      if(val == 2)
                      { /* individual 1 is dominated */
                        flag[i] = 0;
                        break;
                      }
                      if(val == 1)
                      { /* individual 2 is dominated */
                        flag[j] = 0;
                      }
                      if(val == 3)
                      { /* individual 1 & 2 are non dominated */
                        if(flag[j] != 0)
                        {
                          flag[j] = 3;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
        if(j == size)
        {
          pop_ptr->ind[i].rank = rnk;
          flag[i] = 1;
          pop_ptr->rankar[rnk-1][q] = i;
          q++;
        }
      }
    }
    pop_ptr->rankno[rnk-1] = q;

    for(j = 0; j < size; j++)
    {
      if (flag[j] != 1) break;
    }
    if(j == size) break;
 
		pop_ptr->rankar[rnk] = pop_ptr->rankar[rnk-1] + q;
 } 
  pop_ptr->maxrank = rnk;
}



int nsga2::indcmp(double *fit1, double *fit2)
{
  int m,n,value;

  m = 0; n=0;
  while(m < nfunc && fit1[m] <= fit2[m]) 
  {
    if((fit2[m] - fit1[m]) < 1e-7)
    {
      n++;
    }
    m++;
  }
  if(m == nfunc) 
  {
    if(n != nfunc)
    {
      value = 1;     /* 1 dominates 2 */
    }
    else
    {
      value = 3;     /* incomparable */
    }
  }
  else 
  {
    m = 0;
    n = 0;
    while(m < nfunc && fit1[m] >= fit2[m]) 
    {
      if((fit1[m] - fit2[m]) < 1e-7)
      {
        n++;
      }
      m++;
    }
    if(m == nfunc)
    {
      if(n != nfunc)
      {
        value = 2;   /* 2 dominates 1 */
      }
      else
      {
        value = 3;   /* incomparable */
      }
    }
    else
    {
      value = 3;     /* incomparable */
    }
  }
  return(value);
}

