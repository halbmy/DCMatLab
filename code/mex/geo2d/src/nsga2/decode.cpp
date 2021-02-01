/*This is the program to decode the chromosome to get real values*/
#include "nsga2.h"

void nsga2::decode(individual *ind_ptr)
{
	int bit, i, j, k;
	unsigned long bin;
	int *gene_ptr;
	double *real_ptr;
	
	for(i = 0; i < popsize; i++, ind_ptr++)
	{
		if (ind_ptr->change)
		{
			gene_ptr = ind_ptr->genes;
			real_ptr = ind_ptr->xbin;
			
			for(j = 0; j < nchrom; j++, real_ptr++)
			{
				bin = 0UL;
				bit = 0;
				for (k = vlen[j]-1; k >= 0; k--, gene_ptr++)
				{
					bit ^= *gene_ptr;           // Gray -> standard binary
					if (bit) bin += (1 << k);   // standard binary -> unsigned long
				}
				
				*real_ptr = lim_b[j][0] + 
					double(bin) / double((2<<(vlen[j]-1))-1) *
					(lim_b[j][1] - lim_b[j][0]);
			}
		}
	}
}

