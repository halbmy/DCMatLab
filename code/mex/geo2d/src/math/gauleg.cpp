/* For a description of the following routines see
   Numerical Recipes, Press et al.
*/
#include <math.h>

void gauleg(double x1, double x2, double *x, double *w, int n)
/* Given the lower and upper limits of integration x1 and x2 and given n,
   this routine returns arrays x[0..n-1] and w[0..n-1] of length n,
   containing the abscissas and weights of the Gauss-Legendre
   n-point quadrature formula.
*/
{
	double z1, z, xm, xl, pp, p3, p2, p1;
	const double eps = 3.0e-11;
   const double pi = 3.1415926535897931;
	int m, j, i;
	
	m=(n+1)/2;
	xm=0.5*(x2+x1);
	xl=0.5*(x2-x1);
	for(i=1;i<=m;i++)
	{
		z=cos(pi*(i-0.25)/(n+0.5));
/* Starting with the above approximation to the ith root, we enter
   the main loop of refinements by Newton's method.
*/
		do
		{
			p1 = 1.0; p2 = 0.0;
			for(j=1;j<=n;j++)
			{
				p3 = p2; p2 = p1;
				p1 = ((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
			}
/* p1 is now the desired Legendre polynomial. We next compute pp,
   its derivative, by a standard relation involving also p2, the
   polynomial of one lower order.
*/
			pp = n*(z*p1-p2)/(z*z-1.0);
			z1 = z;
			z = z1-p1/pp; /* Newtons method */
		} while (fabs(z-z1) > eps);
/* Scale the root to the desired interval, and put in its
   symmetric counterpart.
*/
		x[i-1] = xm-xl*z;
		x[n-i] = xm+xl*z;
/* Compute the weight and its symmetric counterpart.
*/
		w[i-1] = 2.0*xl/((1.0-z*z)*pp*pp);
		w[n-i] = w[i-1];
	}
}
