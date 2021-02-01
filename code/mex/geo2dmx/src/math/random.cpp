#include "random.h"

double ran1(long *idum) 
/* ``Minimal'' random number generator of Park and Miller with Bays-Durham shuffle and added 
safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint 
values). Call with idum a negative integer to initialize; thereafter, do not alter idum between 
successive deviates in a sequence. RNMX should approximate the largest floating value that is 
less than 1. */
{
	int j; 
	long k; 
  static long iy=0; 
  static long iv[NTAB]; 
	double temp; 
	if (*idum <= 0 || !iy)			// Initialize. 
	{
		if (-(*idum) < 1) *idum=1;	// Be sure to prevent idum = 0. 
		else *idum = -(*idum); 
		for (j=NTAB+7;j>=0;j--)		// Load the shuffle table (after 8 warm-ups).
		{
			k=(*idum)/IQ; 
			*idum=IA*(*idum-k*IQ)-IR*k; 
			if (*idum < 0) *idum += IM; 
			if (j < NTAB) iv[j] = *idum; 
		} 
		iy=iv[0]; 
	} 
	k=(*idum)/IQ;					// Start here when not initializing. 
	*idum=IA*(*idum-k*IQ)-IR*k;		// Compute idum=(IA*idum) % IM without over- 
	if (*idum < 0) *idum += IM;		//  flows by Schrage's method. 
	j=iy/NDIV;						// Will be in the range 0..NTAB-1. 
	iy=iv[j];						// Output previously stored value and refill the 
	iv[j] = *idum;					// 	shuffle table.
	if ((temp=AM*iy) > RNMX) return RNMX;	// Because users don't expect endpoint values. 
	else return temp; 
}

double gasdev(long *idum) 
{
	/*Returns a normally distributed deviate with zero mean and unit variance,
	using ran1(idum) as the source of uniform deviates. */
	static int iset=0; 
	static double gset; 
	double fac,rsq,v1,v2; 
	if (*idum < 0) iset=0;	// Reinitialize. 

	if (iset == 0)			// We don't have an extra deviate handy, so 
	{
		do
		{
			v1=2.0*ran1(idum)-1.0;	// pick two uniform numbers in the square ex- 
			v2=2.0*ran1(idum)-1.0;	//  tending from -1 to +1 in each direction,
			rsq=v1*v1+v2*v2;		// see if they are in the unit circle,
		}
		while (rsq >= 1.0 || rsq == 0.0);	// and if they are not, try again. 
		fac=sqrt(-2.0*log(rsq)/rsq); 
		/* Now make the Box-Muller transformation to get two normal deviates.
		Return one and 	save the other for next time. */
		gset=v1*fac; 
		iset=1;						// Set flag. 
		return v2*fac; 
	}
	else
	{						// We have an extra deviate handy, 
		iset=0;				// so unset the flag, 
		return gset;		// and return it.
	}
}

