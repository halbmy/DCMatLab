/* For a description of the following routines see
   Numerical Recipes, Press et al.
*/
#include <math.h>

double gammln(double xx);


void gaulag(double *x, double *w, int n, double alf)
/* Given alf, the parameter alpha of the Laguerre polynomials, this routine
   returns arrays x[0..n-1] and w[0..n-1] containing the abscissas and weights
   of the n-point Gauss-Laguerre quadrature formula. The smallest abscissa
   is returned in x[0], the largest in x[n-1].
*/
{
	double ai;
	double p1, p2, p3, pp, z, z1;
	const double eps = 3.0e-11;
	int	i, its, j;
	const int maxit = 10;

	z = 0.0;
	for(i=1;i<=n;i++) /* Loop over desired roots */
	{
		if(i==1)
		{
			z = (1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
		}
		else if(i == 2)
		{
			z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
		}
		else
		{
			ai = i-2;
			z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
				(1.0+3.5*ai))*(z-x[i-3])/(1.0+0.3*alf);
		}
		for(its=1;its<=maxit;its++)
		{
			p1 = 1.0; p2 = 0.0;
			for(j=1;j<=n;j++)
			{
				p3 = p2; p2 = p1;
				p1 = ((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
			}
			pp = (n*p1-(n+alf)*p2)/z;
			z1 = z;
			z  = z1-p1/pp;
			if(fabs(z-z1)<=eps) break;
		}
		x[i-1] = z;
		w[i-1] = -exp(gammln(alf+n)-gammln(double(n)))/(pp*n*p2);
	}
}
