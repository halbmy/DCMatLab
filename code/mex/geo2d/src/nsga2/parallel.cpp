#include "nsga2.h"
#include <mpi.h>
#include <iostream.h>
#include <fstream.h>

void nsga2::runs(int igap, int gap, const bool split)
{
	int i, j, k, mypop, prevpop, nextpop, npops, nhypercube;
	int bufsize, bufpos, flag;
   double *ptr1, *ptr2;
	int *inbuf, *outbuf, *neighbours;
	char *processor_name, *pareto_name, *popfile_name;
	MPI_Status status;
	FILE *pareto, *popfile;
	
	processor_name = new char[MPI_MAX_PROCESSOR_NAME];
	pareto_name = new char [128];
	popfile_name = new char[128];
	MPI_Comm_rank(MPI_COMM_WORLD, &mypop);
	MPI_Comm_size(MPI_COMM_WORLD, &npops);
	MPI_Get_processor_name(processor_name, &flag);
   
   if (split)
   {
      /* Set transform matrix for objective functions such as to split Pareto front between 
         different populations */
      transf(mypop + 1, npops);
   }

	/* Print start message */
	if (mypop > 0)
	{
		MPI_Recv(&flag, 0, MPI_INT, mypop-1, 0, MPI_COMM_WORLD, &status);
	}
	cout << "Population " << mypop + 1 << " is evolving on " << processor_name 
		<< endl << flush;
	delete [] processor_name;
	if (mypop < npops - 1)
	{
		MPI_Send(&flag, 0, MPI_INT, mypop+1, 0, MPI_COMM_WORLD);
	}
	else if (npops > 1)
	{
		MPI_Send(&flag, 0, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
	if (mypop == 0 && npops > 1)
		MPI_Recv(&flag, 0, MPI_INT, npops - 1, 0, MPI_COMM_WORLD, &status);
		
	i = npops;
	nhypercube = 0;
	while (i%2 == 0)
	{
		i >>= 1;
		nhypercube++;
	}
	if (i > 1)
		/* Ring topology if npops != 2^n */
	{
		nhypercube = 0;
		prevpop = 0;
		nextpop = 1;
		neighbours = new int[2];
		neighbours[prevpop] = (mypop == 0) ? (npops - 1) : (mypop - 1);
		neighbours[nextpop] = (mypop == npops - 1) ? 0 : (mypop + 1);
	}
	else
		/* Hypercube topology if npops == 2^n */
	{
		prevpop = 0;
		nextpop = 0;
		neighbours = new int[nhypercube];
		for (i = 0; i < nhypercube; i++)
		{
			neighbours[i] = 0;
		}
		flag = 1 << (nhypercube - 1);
		k = mypop;
		for (i = nhypercube - 1; i >= 0; i--)
		{
			if (k / flag)
			{
				for (j = nhypercube - 1; j >= 0; j--)
				{
					if (i != j) neighbours[j] += flag;
				}
				k -= flag;
			}
			else
			{
				neighbours[i] += flag;
			}
			flag >>= 1;
		}
	}
	
	/* Create input and output buffers for Isend/Ireceive*/
	MPI_Pack_size(1 + chrom * popsize, MPI_INT, MPI_COMM_WORLD, &bufsize);
	MPI_Pack_size((nchrom + nvar + nfunc + ncons + 1) * popsize, 
		MPI_DOUBLE, MPI_COMM_WORLD, &flag);
	bufsize += flag;
	inbuf = (int *) malloc(bufsize);
	outbuf = (int *) malloc(bufsize);
	if (inbuf == NULL || outbuf == NULL)
	{
		cerr << "Error allocating buffers for send/receive." << endl;
		MPI_Finalize();
		exit(0);
	}
	
	/* Initialize population, try using binary file "population<mypop+1>.bin" */
	sprintf(popfile_name, "population%d.bin", mypop + 1);
	popfile = fopen(popfile_name, "rb");
	start(0, popfile);
	if (popfile != NULL)
		fclose(popfile);
	
	/* Print status of GA parameters starting with population no. 1 */
	sprintf(pareto_name, "pareto%d.out", mypop + 1);
	if (mypop > 0)
	{
		MPI_Recv(&flag, 0, MPI_INT, mypop - 1, 0, MPI_COMM_WORLD, &status);
	}
	cout << endl << "Population " << mypop + 1 << endl << flush;
	report(stdout,"");
	if (generi)
	{
		cout << "Found file '" 
			<< popfile_name << "'." << endl
			<< "Used for initialization and continued at generation no. " 
			<< generi + 1 << "." << endl << flush;
		pareto = fopen(pareto_name, "a");
	}
	else
	{
		pareto = fopen(pareto_name, "w");
	}
	if (pareto != NULL)
	{
		if(generi) fprintf(pareto, "#\n# Evolution continued ...\n#\n");
		report(pareto, "# ");
		fclose(pareto);
	}
	if (mypop < npops - 1)
	{
		MPI_Send(&flag, 0, MPI_INT, mypop + 1, 0, MPI_COMM_WORLD);
	}
	else if (npops > 1)
	{
		MPI_Send(&flag, 0, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
	if (mypop == 0 && npops > 1)
		MPI_Recv(&flag, 0, MPI_INT, npops - 1, 0, MPI_COMM_WORLD, &status);
	
   /* Check that all populations continue at the same generation */
   MPI_Allreduce(&generi, &k, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
   if (k != generi)
   {
      cerr << "Population no. " << mypop + 1
         << ": generation counter adjusted to " << k
         << " by " << generi - k << " generations."
         << endl << flush;
      pareto = fopen(pareto_name, "a");
      if (pareto != NULL)
      {
         fprintf(stderr, "# Generation counter adjusted to %d by %d generations.\n",
            k, generi-k);
         fclose(pareto);
      }
      generi = k;
   }

	/* Print initial generation to file */
	pareto = fopen(pareto_name, "a");
	if (pareto != NULL)
	{
		print(pareto);
		fclose(pareto);
	}
	
	/* Dump initial population on binary file */
	popfile = fopen(popfile_name, "wb");
	flag = dump(popfile);
	if (!flag)
	{
		fprintf(stderr, "Error occured writing binary dumpfile %s.\n",
			popfile_name);
	}
	fclose(popfile);

	/* Adjust counter for populational exchange if continuing another evolutional run */
	if (generi)
	{
		igap -= generi;
		while (igap < 0) igap += gap;
	}

	/* Generation loop */
	while (generi < gener) 
	{
		if (generi + igap > gener) igap = gener - generi;

		/* Run evolution for igap generations */
		step(igap);

		/* Some counters */
		generi += igap;
		igap = gap;
		
		/* Dump population on binary file */
		popfile = fopen(popfile_name, "wb");
		flag = dump(popfile);
		if (!flag)
		{
			fprintf(stderr, "Error occured writing binary dumpfile %s.\n",
				popfile_name);
		}
		fclose(popfile);
		
		/* Print population to file */
		pareto = fopen(pareto_name, "a");
		if (pareto != NULL)
		{
			print(pareto);
			fclose(pareto);
		}

		if (mypop == 0) cout << endl << generi << " of " << gener << " generations completed."
			<< " - " << double(generi*100) / double(gener) << " % done." << endl << flush;

		/* Send current Pareto-front to next process and receive other
		front from previous process */
      if (npops > 1)
      {
         bufpos = pack(outbuf, bufsize, 1);
         MPI_Sendrecv(outbuf, bufpos, MPI_PACKED, neighbours[nextpop], 1,
            inbuf, bufsize, MPI_PACKED, neighbours[prevpop], 1, 
            MPI_COMM_WORLD,&status);
         if (nhypercube)
         {
            nextpop = (nextpop == nhypercube - 1) ? 0 : nextpop + 1;
            prevpop = (prevpop == nhypercube - 1) ? 0 : prevpop + 1;
         }
         k = unpack(inbuf, bufsize, popsize);
         keepalive(popsize + k);
      }
	}
	
	/* Collect all populations on process 0 */
	if (mypop > 0)
	{
		bufpos = pack(outbuf, bufsize, pop_ptr->maxrank);
		MPI_Send(outbuf, bufpos, MPI_PACKED, 0, 2,
			MPI_COMM_WORLD);
	}
	else
	{
      if (split)
      {
         transf(); // set transform matrix to unit -> use raw objective values
         // recalculate local fitness values (reduces to just copying since T is the unit matrix)
         for (i = 0; i < popsize; i++)
         {
            ptr1 = pop_ptr->ind[i].fitness;
            ptr2 = pop_ptr->ind[i].object;
            for (j = 0; j < nfunc; j++)
               *ptr1++ = *ptr2++;
         }
      }

		k = popsize;        // save old popsize since enlarge will change popsize
		enlarge(k * npops); // changes popsize to (npops * popsize)

		for (i = 1; i < npops; i++)
		{
			MPI_Recv(inbuf, bufsize, MPI_PACKED, MPI_ANY_SOURCE, 2, 
				MPI_COMM_WORLD, &status);
			unpack(inbuf, bufsize, i * k);
		}

		/* Find the global ranks */
		if (ncons == 0) rank(popsize);
		else            rankc(popsize);
		/* Sharing the fitness to get the dummy fitness */
		if (nsharepar)  share2(popsize);
		else            share(popsize);
	
	    /* Print results on output file  */
		FILE *stream = fopen("results.out","w");
		if (stream != NULL)
		{
			print(stream);
			fclose(stream);
		}
		else
		{
			cerr << "Error opening output file 'results.out'." << endl;
		}
	}
	
	/* Free send/receive buffers */
	free(inbuf);
	free(outbuf);
	delete [] pareto_name;
	delete [] popfile_name;
	delete [] neighbours;	
}

int nsga2::pack(int *buf, int bufsize, int rank)
{
  int i, n;
  int bufpos = 0;
  individual *ind_ptr;

  n = 0;
  ind_ptr = pop_ptr->ind;
  for (i = 0; i < popsize; i++)
  {
	  if (ind_ptr->rank <= rank)
	  {
		  flag[i] = 1;
		  n++;
	  }
	  else
	  {
		  flag[i] = 0;
	  }
	  ind_ptr++;
  }

  MPI_Pack(&n, 1, MPI_INT, 
	  buf, bufsize, &bufpos, MPI_COMM_WORLD);

  ind_ptr = pop_ptr->ind;
  for (i = 0; i < popsize; i++)
  {
	  if (flag[i] && chrom > 0)
		  MPI_Pack(ind_ptr->genes, chrom, MPI_INT,
		  buf, bufsize, &bufpos, MPI_COMM_WORLD);
	  ind_ptr++;
  }
  ind_ptr = pop_ptr->ind;
  for (i = 0; i < popsize; i++)
  {
	  if (flag[i])
	  {
		  if (nchrom > 0)
			  MPI_Pack(ind_ptr->xbin, nchrom, MPI_DOUBLE, 
			  buf, bufsize, &bufpos, MPI_COMM_WORLD);
		  if (nvar > 0)
			  MPI_Pack(ind_ptr->xreal, nvar, MPI_DOUBLE, 
			  buf, bufsize, &bufpos, MPI_COMM_WORLD);
		  if (nfunc > 0)
			  MPI_Pack(ind_ptr->object, nfunc, MPI_DOUBLE, 
			  buf, bufsize, &bufpos, MPI_COMM_WORLD);
		  if (ncons > 0)
		  {
			  MPI_Pack(ind_ptr->constr, ncons, MPI_DOUBLE, 
				  buf, bufsize, &bufpos, MPI_COMM_WORLD);
			  MPI_Pack(&(ind_ptr->error), 1, MPI_DOUBLE, 
				  buf, bufsize, &bufpos, MPI_COMM_WORLD);
		  }
	  }
	  ind_ptr++;
  }
  return(bufpos);
}

int nsga2::unpack(int *buf, int bufsize, int pos)
{
  int i, j, k, n;
  int bufpos = 0;
  double *ptr1, *ptr2;
  individual *ind_ptr;
  
  MPI_Unpack(buf, bufsize, &bufpos,
	  &n, 1, MPI_INT, MPI_COMM_WORLD);
  ind_ptr = pop_ptr->ind + pos;
  for (i = 0; i < n; i++)
  {
    if (chrom > 0)
      MPI_Unpack(buf, bufsize, &bufpos, 
        ind_ptr->genes, chrom, MPI_INT, MPI_COMM_WORLD);
    ind_ptr++;
  }
  ind_ptr = pop_ptr->ind + pos;
  for (i = 0; i < n; i++)
  {
    if (nchrom > 0)
      MPI_Unpack(buf, bufsize, &bufpos, 
        ind_ptr->xbin, nchrom, MPI_DOUBLE, MPI_COMM_WORLD);
    if (nvar > 0)
      MPI_Unpack(buf, bufsize, &bufpos, 
        ind_ptr->xreal, nvar, MPI_DOUBLE, MPI_COMM_WORLD);
    if (nfunc > 0)
    {
      ptr1 = ind_ptr->object;
      ptr2 = ind_ptr->fitness;
      MPI_Unpack(buf, bufsize, &bufpos,
        ptr1, nfunc, MPI_DOUBLE, MPI_COMM_WORLD);
      for(j = 0; j < nfunc; j++)
      {
         ptr2[j] = 0.0;
         for (k = 0; k < nfunc; k++)
            ptr2[j] += T[j][k] * (ptr1[k] - lim_f[k][0]) / (lim_f[k][1] - lim_f[k][0]);
      }
    }
    if (ncons > 0)
	{
      MPI_Unpack(buf, bufsize, &bufpos, 
        ind_ptr->constr, ncons, MPI_DOUBLE, MPI_COMM_WORLD);
      MPI_Unpack(buf, bufsize, &bufpos, 
        &(ind_ptr->error), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	}
    ind_ptr->change = 0;
    ind_ptr++;
  }
  return(n);
}

