/*
 * computerhoa.cpp
 */

#include "geo2d.h"

double geo2d::computerhoa(double *rhoa)
{
   double p11, p12, p21, p22, voltg, voltr, recp = 0.0, rhog, rhor;
   int i, is, ip;
   double *pot;

   pot = (double *) CALLOC(nsrc * npot, sizeof(double));

   computephi(pot);

   for (i = 0; i < ndata; i++)
   {
      is = is1[i];
      if (is--)
      {
         is *= npot;
         ip = ip1[i];
         if (ip--)
         {
            p11 = pot[is+ip];
         }
         else
         {
            p11 = 0.0;
         }
         ip = ip2[i];
         if (ip--)
         {
            p12 = pot[is+ip];
         }
         else
         {
            p12 = 0.0;
         }
      }
      else
      {
         p11 = 0.0;
         p12 = 0.0;
      }
      is = is2[i];
      if (is--)
      {
         is *= npot;
         ip = ip1[i];
         if (ip--)
         {
            p21 = pot[is+ip];
         }
         else
         {
            p21 = 0.0;
         }
         ip = ip2[i];
         if (ip--)
         {
            p22 = pot[is+ip];
         }
         else
         {
            p22 = 0.0;
         }
      }
      else
      {
         p21 = 0.0;
         p22 = 0.0;
      }
      voltg = (p11 - p12 - p21 + p22);
      rhog = confac[i] * voltg;
      /* Check accuracy by reciprocity test */
      is = ip1[i];
      if (is--)
      {
         is *= npot;
         ip = is1[i];
         if (ip--)
         {
            p11 = pot[is+ip];
         }
         else
         {
            p11 = 0.0;
         }
         ip = is2[i];
         if (ip--)
         {
            p12 = pot[is+ip];
         }
         else
         {
            p12 = 0.0;
         }
      }
      else
      {
         p11 = 0.0;
         p12 = 0.0;
      }
      is = ip2[i];
      if (is--)
      {
         is *= npot;
         ip = is1[i];
         if (ip--)
         {
            p21 = pot[is+ip];
         }
         else
         {
            p21 = 0.0;
         }
         ip = is2[i];
         if (ip--)
         {
            p22 = pot[is+ip];
         }
         else
         {
            p22 = 0.0;
         }
      }
      else
      {
         p21 = 0.0;
         p22 = 0.0;
      }
      voltr = p11 - p12 - p21 + p22;
	  rhor = confac[i] * voltr;
      recp = max((rhog != 0.0 ? fabs(rhor/rhog-1.0) : HUGE_VAL), recp);
	  if (rhor >= 0.0) /* geometric mean might be more accurate */
	  {
		  if (rhog >= 0.0)         rhoa[i] = sqrt(rhor * rhog);
		  else if (rhog >= -rhor)  rhoa[i] = rhog + sqrt(-rhor * rhog);
		  else                     rhoa[i] = rhor - sqrt(-rhor * rhog);
	  }
	  else
  	  {
		  if (rhog <= 0.0)         rhoa[i] = -sqrt(rhor * rhog);
		  else if (rhog <= -rhor)  rhoa[i] = rhog - sqrt(-rhor * rhog);
		  else                     rhoa[i] = rhor + sqrt(-rhor * rhog);
	  }
   }

   FREE(pot);

   return(recp);
}

