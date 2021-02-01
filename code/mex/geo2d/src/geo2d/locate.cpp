/*
 * locate.cpp
 */

#include "geo2d.h"

int geo2d::locate(double *xx, int n, double x)
/* 
 * Given an array xx[0..n-1], and given a value x, returns a value j
 * such that x is between xx[j] and xx[j+1]. xx must be monotonic,
 * either increasing or decreasing. j=-1 or j=n-1 is returned to indicate
 * that x is out of range.
 */
{
  int ju,jm,jl,ascnd;
  xx--;     /* C-arrays are zero based */
  jl = 0;   /* Initialize lower */
  ju = n+1; /* and upper limits.*/
  ascnd = (xx[n] >= xx[1]);
  while (ju - jl > 1)   /* If we are not yet done, */
  {
    jm = (ju+jl) >> 1;  /* compute a midpoint, */
    if (x >= xx[jm] == ascnd)
    {
      jl = jm;          /* and replace either the lower limit */
    }
    else
    {
      ju = jm;          /* or the upper limit, as appropriate. */
    }
  }         /* Repeat until the test condition is satisfied. */
  if (x == xx[1]) /* Then set the output and return */
  {
    jl = 1;
  }
  if(x == xx[n])
  {
    jl = n-1;
  }
  return(jl-1);    /* C-arrays are zero based */
}

