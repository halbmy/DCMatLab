/*
 * solvelineq.cpp
 */

#include "geo2d.h"
#include "lapack.h"

void geo2d::solvelineq()
/*
 * Solve system of linear equations C*phi = b using LAPACK routine
 * dpbsv (Cholesky factorization, and backward-/forward substitution)
 */
{
  int n, kd, nrhs, ldab, ldb, info;
  char uplo;

  uplo = 'U';
	n = nx*nz;
	kd = nz;
	nrhs = nsrc;
	ldab = nz+1;
	ldb = n;	
  info=-32768;

	F77_NAME(dpbsv)(&uplo, &n, &kd, &nrhs, C, &ldab, b, &ldb, &info);
			
	if (info<0){
		PRINTF("Error while solving linear equations.\n"
			"%d-th argument of dpbsv had an illegal value.\n",
			-1*info);
		exit(0);
	}
	if (info>0){
		PRINTF("Error while solving linear equations.\n"
			"The leading minor of order %d"
			" of matrix C is not\n"
			"positive definite, so the factorization could not be\n"
			"completed, and the solution has not been computed.\n",
      info);
		exit(0);
	}
}

