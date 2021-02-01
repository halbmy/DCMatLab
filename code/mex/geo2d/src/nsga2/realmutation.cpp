/* This is the module used to formulate the real mutation routine */
#include "nsga2.h"

void nsga2::realmutate(individual *ind_ptr)
{
  bool flag = false;
  int i,j;
  double rnd,delta,y,yl,yu,g1,g2,expp = dim + 1.0;
  double *ptr;

  for(j = 0;j < popsize;j++)
  {
    ptr = ind_ptr->xreal;

    for (i = 0;i < nvar; i++)
    {
      /*For each variable find whether to do mutation or not*/
      rnd = ran1(&seed);
      if(rnd <= pmut_r)
      {
        y = *ptr;
		yl = lim_r[i][0];
		yu = lim_r[i][1];

		if (ans)
		{
			g1 = 1.0 / pow(1.0 + y - yl, expp);
			g2 = 1.0 / pow(1.0 - y + yu, expp);
		}
		else
		{
			g1 = 0.0;
			g2 = 0.0;
		}

		rnd = ran1(&seed);
		if (rnd <= (1.0 - g1) / (2.0 - g1 - g2))
		{
			delta = 1.0 - 1.0 / pow((2.0 - g1) * (1.0 - rnd) + g2 * rnd, 1.0 / expp);
		}
		else
		{
			delta = 1.0 / pow(g1 * (1.0 - rnd) + (2.0 - g2) * rnd, 1.0 / expp) - 1.0;
		}
		if (!ans)
		{
			delta *= yu - yl;
		}

		y += delta;
		if (ans && y < yl) y = yl;
		if (ans && y > yu) y = yu;
		*ptr = y;
		flag = true;
      } /*endif mutation of i-th variable*/
      ptr++; /*next variable*/
    } /*end of variables loop*/
    if (flag)
    {
        ind_ptr->change = 1;
        nmut++;
		flag = false;
    }
    ind_ptr++;
  } /*end of population loop*/
}

