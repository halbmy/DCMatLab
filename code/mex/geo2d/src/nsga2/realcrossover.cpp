/*This is the file used for crossover for Real Coded GA*/
#include "nsga2.h"

void nsga2::sbcross(individual *ind_ptr)
{
  bool flag = false;
  int i,j;
  double mean,diff,rnd,alpha1,alpha2,beta1,beta2,xl,xu,x1,x2,y1,y2,expp = di + 1.0;
  double *ptr1,*ptr2;
  individual *iptr1 = ind_ptr,
    *iptr2 = ind_ptr + 1;

  for(i = 0; i < popsize/2; i++)
  {
    ptr1 = iptr1->xreal;
    ptr2 = iptr2->xreal;
 
    /*Check whether the cross-over to be performed*/
    rnd = ran1(&seed);
    if(rnd <= pcross)
    {
	  ncross++;
      for(j = 0; j < nvar; j++) /*Loop over no of variables*/
      { 
		x1 = *ptr1;
		x2 = *ptr2;
		xl = lim_r[j][0];
		xu = lim_r[j][1];
		
		mean = 0.5 * (x1 + x2);
		diff = 0.5 * fabs(x1 - x2);
		
		if (diff > 1.0e-8)
		{
			if (ans)
			{
				beta1 = (mean - xl) / diff;
				beta2 = (xu - mean) / diff;
				
				alpha1 = 2.0 - 1.0 / pow(beta1,expp);
				alpha2 = 2.0 - 1.0 / pow(beta2,expp);
			}
			else
			{
				alpha1 = 2.0;
				alpha2 = 2.0;
			}

			rnd = ran1(&seed);
			beta1 = ((rnd <= 1.0 / alpha1) ?
				pow(rnd * alpha1,               1.0 / expp) : 
			    pow(1.0 / (2.0 - rnd * alpha1), 1.0 / expp));
			beta2 = ((rnd <= 1.0 / alpha2) ?
				pow(rnd * alpha2,               1.0 / expp) : 
			    pow(1.0 / (2.0 - rnd * alpha2), 1.0 / expp));
			
		    y1 = mean - beta1 * diff;
			y2 = mean + beta2 * diff;
			if (ans && y1 < xl) y1 = xl;
			if (ans && y2 < xl) y2 = xl;
			if (ans && y1 > xu) y1 = xu;
			if (ans && y2 > xu) y2 = xu;

			if (x1 < x2)
			{
				*ptr1 = y1;
				*ptr2 = y2;
			}
			else
			{
				*ptr1 = y2;
				*ptr2 = y1;
			}
			flag = true;
		}
        ptr1++;
        ptr2++;
      } /* end of variables loop */
    } /* endif crossover no. i*/
    if (flag)
    {
        iptr1->change = 1;
        iptr2->change = 1;
		flag = false;
    }
    iptr1 += 2;
    iptr2 += 2;
  } /* end of population loop */
}

