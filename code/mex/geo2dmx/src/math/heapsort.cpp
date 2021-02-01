void heapsort(unsigned long n, double *ra)
/* Sorts an array ra[0..n-1] into ascending numerical order using the
   Heapsort algorithm. n is input; ra is replaced on output by its 
	 sorted rearrangement.*/
{
  unsigned long i,ir,j,l;
  double rra;
  if (n < 2) return;
  l = n >> 1;
  ir = n - 1;
  /* The index l will be decremented from its initial value down to 0
	   during the "hiring" (heap creation) phase. Once it reaches 0,
		 the index ir will be decremented from its initial value down to 0
		 during the "retirement-and-promotion" (heap selection) phase. */
  for (;;)
	{
    if (l > 0) /* Still in hiring phase. */
		{
      rra = ra[--l];
		}
		else       /* In retirement-and-promotion phase. */
		{
      rra = ra[ir];   /* Clear a space at end of array. */
      ra[ir] = ra[0]; /* Retire the top of the heap into it. */
      if (--ir == 0)  /* Done with the last promotion. */
      {
				ra[0] = rra;  /* The least competent worker of all! */
        break;
			}
		}
    i = l;         /* Whether in the hiring phase or promotion phase, we */
    j = l + l + 1; /* here set up to sift down element rra to its proper */
                   /* level. */
    
    while (j <= ir)
		{
      if (j < ir && ra[j] < ra[j+1]) 
			{
			  j++; /*Compare to the better underling. */
			}
      if (rra < ra[j]) /* Demote rra. */
			{
        ra[i] = ra[j];
        i = j;
        j = j + j + 1;
			}
			else
			{
				break; /* Found rra's level. Terminate the sift-down. */
			}
		}
    ra[i] = rra; /* Put rra into its slot. */
	}
}

