/********************************************************************
*            Objective functions and constraints                    *
********************************************************************/
#include "dc2d.h"

void application::eval(double *x, int nx, double *f, int nf, double *c, int nc)
{
  int i, j, n, cx, cz;
  double err;

  cx = NXGrid-1;
  cz = NZGrid-1;
  n = cx * cz;

  /* powerscale and sort arguments */
  for (i = 0; i < nx; i++)
  {
    x[i] = pow(10.0, x[i]);
  }
  for (i = 0; i < n; i++)
  {
    j = IConductivity[i];
    if(j < nx)
    {
      Conductivity[i] = x[j];
    }
    else
    {
      fprintf(stderr, "Index for conductivity exceeds number of parameters "
        "provided by the GA.\n");
      exit(0);
    }
  }

  if(nf != 2)
  {
    fprintf(stderr, "Number of functions in app_eval doesn't match number of "
      "functions nfunc considered by the GA.\n");
    exit(0);
  }

  datamisfit(f,&err);

  /* Constrain erroneous modelling */
  if(nc == 1)
  {
	   *c = err;
  }
  else if (nc != 0)
  {
    fprintf(stderr, "Number of constrains in app_eval doesn't match number of "
      "constrains ncons considered by the GA.\n");
    exit(0);
  }

  modelrestriction(f+1);
}

#ifndef LOGRHOA

void application::datamisfit(double *f, double *e)
/* 1st objective function: data misfit */
{
   double t1, t2, val;
   int i;
   /* calculate synthetic data set DataSyn */
   geo2dobj->setmodel(NXGrid, XGrid, NZGrid, ZGrid, Conductivity);
   *e = geo2dobj->computerhoa(DataSyn);
   
   t2 = 0.0;
   if (obj1 == 1.0) /* l1-norm */
   {
      for (i = 0; i < NoData; i++)
      {
         t1 = (DataSyn[i] - DataVal[i]) / DataStd[i];
         t2 += fabs(t1);
      }
      val = (NoData > 1) ? (t2 / sqrt(double(NoData-1))) : HUGE_VAL;
   }
   else if (obj1 == 2.0) /* l2-norm */
   {
      for (i = 0; i < NoData; i++)
      {
         t1 = (DataSyn[i] - DataVal[i]) / DataStd[i];
         t2 += t1 * t1;
      }
      val = (NoData > 1) ? (sqrt(t2 / double(NoData-1))) : HUGE_VAL;
   }
   else if (obj1 == HUGE_VAL) /* maximum-norm */
   {
      for (i = 0; i < NoData; i++)
      {
         t1 = (DataSyn[i] - DataVal[i]) / DataStd[i];
         t2 = max(fabs(t1), t2);
      }
      val = (NoData > 1) ? (t2 / sqrt(double(NoData-1))) : HUGE_VAL;
   }
   else   /* lp-norm, p = obj1 != 1, 2, Inf */
   {
      for (i = 0; i < NoData; i++)
      {
         t1 = (DataSyn[i] - DataVal[i]) / DataStd[i];
         t2 += pow(fabs(t1), obj1);
      }
      val = (NoData > 1) ? (pow(t2,1.0/obj1) / sqrt(double(NoData-1))) : HUGE_VAL;
   }
   *f = val;
}

#else

void application::datamisfit(double *f, double *e)
/* 1st objective function: data misfit */
{
   double t0, t1, t2, val, err;
   int i, n;
   /* calculate synthetic data set DataSyn */
   geo2dobj->setmodel(NXGrid, XGrid, NZGrid, ZGrid, Conductivity);
   err = geo2dobj->computerhoa(DataSyn);
   
   t2 = 0.0;
   n = 0;
   if (obj1 == 1.0) /* l1-norm */
   {
      for (i = 0; i < NoData; i++)
      {
         t0 = DataSyn[i];
         if (t0 > 0.0)
         {
            t1 = (log(DataSyn[i]) - DataVal[i]) / DataStd[i];
            t2 += fabs(t1);
            n++;
         }
      }
      val = (n > 1) ? (t2 / sqrt(double(n-1))) : HUGE_VAL;
   }
   else if (obj1 == 2.0) /* l2-norm */
   {
      for (i = 0; i < NoData; i++)
      {
         t0 = DataSyn[i];
         if (t0 > 0.0)
         {
            t1 = (log(DataSyn[i]) - DataVal[i]) / DataStd[i];
            t2 += t1 * t1;
            n++;
         }
      }
      val = (n > 1) ? (sqrt(t2 / double(n-1))) : HUGE_VAL;
   }
   else if (obj1 == HUGE_VAL) /* maximum-norm */
   {
      for (i = 0; i < NoData; i++)
      {
         t0 = DataSyn[i];
         if (t0 > 0.0)
         {
            t1 = (log(DataSyn[i]) - DataVal[i]) / DataStd[i];
            t2 = max(fabs(t1), t2);
            n++;
         }
      }
      val = (n > 1) ? (t2 / sqrt(double(n-1))) : HUGE_VAL;
   }
   else   /* lp-norm, p = obj1 != 1, 2, Inf */
   {
      for (i = 0; i < NoData; i++)
      {
         t0 = DataSyn[i];
         if (t0 > 0.0)
         {
            t1 = (log(DataSyn[i]) - DataVal[i]) / DataStd[i];
            t2 += pow(fabs(t1), obj1);
            n++;
         }
      }
      val = (n > 1) ? (pow(t2,1.0/obj1) / sqrt(double(n-1))) : HUGE_VAL;
   }
   *f = val;
   
  	n -= NoData;
   if (err >= 1.0) n--;
   *e = (n ? double(n) : err);
}

#endif

void application::modelrestriction(double *f)
{
   double t0, t1, t2;
   int cx, cz, i, j, k, n;

  cx = NXGrid-1;
  cz = NZGrid-1;

  /* 2nd objective function: constrain models */
  switch (obj2)
  {
  case 1: /* smoothness: 2nd FD derivative of conductivity array */
    n = cx * cz;
    for (i = 0; i < n; i++)
    {
       Conductivity[i] = log(Conductivity[i]);
    }
    // lower right corner
    i = cx - 1;
    j = cz - 1;
    k = i * cz;
    t1 = CXM[i] * Conductivity[k+j-cz]
       + CZM[j] * Conductivity[k+j-1]
       - (CXM[i] + CZM[j]) * Conductivity[k+j];
    t2 = t1 * t1;
    // right edge
    for (j--; j > 0; j--)
    {
       t1 = CXM[i] * Conductivity[k+j-cz]
          + CZM[j] * Conductivity[k+j-1] + CZP[j] * Conductivity[k+j+1]
          - (CXM[i] + CZM[j] + CZP[j]) * Conductivity[k+j];
       t2 += t1 * t1;
    }
    // upper right corner
    t1 = CXM[i] * Conductivity[k+j-cz]
       + CZP[j] * Conductivity[k+j+1]
       - (CXM[i] + CZP[j]) * Conductivity[k+j];
    t2 += t1 * t1;
    for (i = cx-2; i > 0; i--)
    {
       // lower edge
       j = cz - 1;
       k = i * cz;
       t1 = CXM[i] * Conductivity[k+j-cz] + CXP[i] * Conductivity[k+j+cz]
          + CZM[j] * Conductivity[k+j-1]
          - (CXM[i] + CXP[i] + CZM[j]) * Conductivity[k+j];
       t2 += t1 * t1;
       // interiour cells
       for (j--; j > 0; j--)
       {
          t1 = CXM[i] * Conductivity[k+j-cz] + CXP[i] * Conductivity[k+j+cz]
             + CZM[j] * Conductivity[k+j-1] + CZP[j] * Conductivity[k+j+1]
             - (CXM[i] + CXP[i] + CZM[j] + CZP[j]) * Conductivity[k+j];
          t2 += t1 * t1;
       }
       // upper edge
       t1 = CXM[i] * Conductivity[k-cz] + CXP[i] * Conductivity[k+cz]
          + CZP[j] * Conductivity[k+j+1]
          - (CXM[i] + CXP[i] + CZP[j]) * Conductivity[k];
       t2 += t1 * t1;
    }
    // lower left corner
    j = cz - 1;
    t1 = CXP[i] * Conductivity[j+cz]
       + CZM[j] * Conductivity[j-1]
       - (CXP[i] + CZM[j]) * Conductivity[j];
    t2 += t1 * t1;
    // left edge
    for (j--; j > 0; j--)
    {
       t1 = CXP[i] * Conductivity[j+cz]
          + CZM[j] * Conductivity[j-1] + CZP[j] * Conductivity[j+1]
          - (CXP[i] + CZM[j] + CZP[j]) * Conductivity[j];
       t2 += t1 * t1;
    }
    // upper left corner
    t1 = CXP[i] * Conductivity[j+cz]
       + CZP[j] * Conductivity[j+1]
       - (CXP[i] + CZP[j]) * Conductivity[j];
    t2 += t1 * t1;
    *f = sqrt(t2 / double(n-1));
    break;
  case 2: /* number of distinct non-connected conductivity blocks */
/* NEEDS TO BE REWRITTEN FOR OPERATION ON ARRAY CONDUCTIVITY */
    k = 0;
/*    for (i = 0; i < block.n; i++)
    {
      block.flag[i] = 0;
    }
    for (i = 0; i < block.n; i++)
    {
      if (block.flag[i] == 0)
      {
        k++;
        block.flag[i] = 1;
        is1 = 0;
        is2 = 0;
        block.list[is2++] = i;
        do
        {
          m = block.list[is1++];
          for (j = 0; j < block.nnb[m]; j++)
          {
            n = block.inb[m][j];
            t1 = x[block.cnd[m]] / x[block.cnd[n]];
            if (0.99 < t1 && t1 < 1.01 && block.flag[n] == 0)
            { 
              block.flag[n] = 1;
              block.list[is2++] = n;
            }
          }
        }
        while(is1 < is2);
      }
    }*/
    *f = double(k);
    break;
  case 3: /* l1-norm of 1st derivative of conductivity array */
    /* using logarithm of conductivity */
    for (i = 0; i < cx*cz; i++) 
    {
       Conductivity[i] = log(Conductivity[i]);
    }
    t1 = 0.0;
    for (i = 0; i < cx; i++)
    {
      t0 = 0.0;
      for (j = i*cz; j < (i+1)*cz-1; j++)
      {
        t0 += fabs(Conductivity[j+1] - Conductivity[j]);
      }
      t1 += t0 * (XGrid[i+1] - XGrid[i]);
    }
    t2 = 0.0;
    for (i = 0; i < cz; i++)
    {
      t0 = 0.0;
      for (j = i; j < (cx-1)*cz; j+=cz)
      {
        t0 += fabs(Conductivity[j+cz] - Conductivity[j]);
      }
      t2 += t0 * (ZGrid[i+1] - ZGrid[i]);
    }
    t0 = (cz-1) * (XGrid[cx]-XGrid[0]) + (cx-1) * (ZGrid[cz]-ZGrid[0]);
//    *f = (t1 * (cz-1) * (XGrid[cx]-XGrid[0]) + t2 * (cx-1) * (ZGrid[cz]-ZGrid[0])) / t0;
    *f = sqrt(t1 + t2);
    break;
  case 4: /* l2-norm of discrete first derivative using normalized distances (==1) */
    n = cx * cz;
    for (i = 0; i < n; i++)
    {
       Conductivity[i] = log(Conductivity[i]);
    }
    t1 = 0.0; // differences in x-direction
    for (i = 1; i < cx; i++)
    {
       for (j = 0; j < cz; j++)
       {
           t0 = Conductivity[i*cz+j] - Conductivity[(i-1)*cz+j];
           t1 += t0 * t0;
       }
    }
    t2 = 0.0; // differences in z-direction
    for (i = 0; i < cx; i++)
    {
       for (j = 1; j < cz; j++)
       {
           t0 = Conductivity[i*cz+j] - Conductivity[i*cz+j-1];
           t2 += t0 * t0;
       }
    }
//    *f = sqrt( t1 / double((cx-1) * cz) + t2 / double(cx * (cz - 1)) );
    *f = sqrt(t1 + t2);
    break;
  case 5: /* l2-norm of discrete second derivative using normalized distances (==1) */
    n = cx * cz;
    for (i = 0; i < n; i++)
    {
       Conductivity[i] = log(Conductivity[i]);
    }
    t1 = 0.0; // differences in x-direction
    for (j = 0; j < cz; j++)
    {
       t0 = Conductivity[j] - Conductivity[cz+j];
       t1 += t0 * t0;
    }
    for (i = 1; i < cx-1; i++)
    {
       for (j = 0; j < cz; j++)
       {
           t0 = 2.0 * Conductivity[i*cz+j] - Conductivity[(i-1)*cz+j] - Conductivity[(i+1)*cz+j];
           t1 += t0 * t0;
       }
    }
    for (j = 0; j < cz; j++)
    {
       t0 = Conductivity[i*cz+j] - Conductivity[(i-1)*cz+j];
       t1 += t0 * t0;
    }
    t2 = 0.0; // differences in z-direction
    for (i = 0; i < cx; i++)
    {
       t0 = Conductivity[i*cz] - Conductivity[i*cz+1];
       t2 += t0 * t0;
       for (j = 1; j < cz-1; j++)
       {
          t0 = 2.0 * Conductivity[i*cz+j] - Conductivity[i*cz+j-1] - Conductivity[i*cz+j+1];
          t2 += t0 * t0;
       }
       t0 = Conductivity[i*cz+j] - Conductivity[i*cz+j-1];
       t2 += t0 * t0;
    }
//    *f = sqrt( t1 / double(n) + t2 / double(n) );
    *f = sqrt(t1 + t2);
    break;
  default:
     *f = 0.0;
  }
}