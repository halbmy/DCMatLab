/*
 * lapack.h
 */

#ifdef WIN32
#define F77_NAME(x) x
#else
#define F77_NAME(x) x##_
#endif

extern "C"
{
	/* dpbsv: Solve system of linear equations by Cholesky decomposition
						for a banded, symmetric positive definite matrix */
	void F77_NAME(dpbsv)(char *uplo, int *n, int *kd, int *nrhs, double *ab,
		int *ldab, double *b, int *ldb, int *info);
}

