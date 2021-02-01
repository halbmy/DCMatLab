#include "dc2d.h"

void application::setrange(double **range, int nf)
{
   double s, t;
   int cx, cz, i, j, k, n;

   if (nf == 2)
   {
      cx = NXGrid - 1;
      cz = NZGrid - 1;
      n = cx * cz;
      /* First, calculate maximum datamisfit and minimum model restriction using
         optimum homogeneous halfspace model */
      s = 0.0;
      t = 0.0;
      for (i = 0; i < NoData; i++)
      {
         s += 1.0 / DataStd[i];
         t += DataVal[i] / DataStd[i];
      }
#ifndef LOGRHOA
      s /= t; // conductivity of best fitting homogeneous halfspace
#else
      s = exp(-t/s);
#endif
      
      for (i = 0; i < n; i++) // set conductivity array
         Conductivity[i] = s;
      
      datamisfit(&s, &t);   // calculate data misfit and reciprocity error
      modelrestriction(&t); // calculate model restriction
      
      range[0][1] = s;
      range[1][0] = t;
      
      /* Second, estimate maximum model restriction using 
         chequerboard pattern of minimum and maximum conductivity values.
         This is an extreme worst case estimation not corresponding with
         minimum data misfit at all. */
      for (i = 0, k = 0; i < cx; i++)
      {
         for (j = 0; j < cz; j++, k++)
         {
            Conductivity[k] = ((i+j) % 2 ? C_max : C_min);
         }
      }
      
      modelrestriction(&t);
      range[1][1] = t;

      /* Finally, a trivial estimate for minimum data misfit -- we are using
         norms to measure data misfit */
      range[0][0] = 0.0;
   }
}