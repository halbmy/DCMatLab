#include <math.h>

double pythag(double a, double b)
/* Computes sqrt(a * a + b * b) without destructive underflow or overflow */
{
  double absa, absb;

  absa = fabs(a);
  absb = fabs(b);
  if (absa > absb)
  {
    absb /= absa;
    return( absa * sqrt(1.0 + absb*absb) );
  }
  else
  {
    if (absb == 0.0)
    {
      return( 0.0 );
    }
    else
    {
      absa /= absb;
      return( absb * sqrt(1.0 + absa*absa) );
    }
  }
}

