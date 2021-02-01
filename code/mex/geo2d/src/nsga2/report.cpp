/*This program subroutine is used to print the report*/
#include "nsga2.h"

void nsga2::report(int t)
{
	int i,j,n;
    int *iptr;
	char name[256];
	double *ptr;
	individual *ind_ptr;
	static FILE *outrep;
	FILE *outfit, *outvar, *pareto;
	
	/*********************************************************************
	*                       Initial generation                           *
	*********************************************************************/
	if (t==0)
	{
		/*Opening the report file*/
		outrep = fopen("report.out", "w");
		
		/*Print the GA parameters and problem parameters in the file 
		report.out */
		fprintf(outrep,"# NSGA-2 PARAMETERS\n");
		fprintf(outrep,"# -----------------\n");
		fprintf(outrep,"# Population Size ->%d\n",popsize);
		fprintf(outrep,"# Chromosome Length ->%d\n",chrom);
		fprintf(outrep,"# No. of generations ->%d\n",gener);
		fprintf(outrep,"# No. of Functions ->%d\n",nfunc);
		fprintf(outrep,"# No. of Constraints ->%d\n",ncons);
		if (nchrom > 0)
		{
			fprintf(outrep,"# No. of binary-coded variables -> %d\n",nchrom);
		}
		if (nvar > 0)
		{
			fprintf(outrep,"# No. of real-coded variables -> %d\n",nvar);
		}
		fprintf(outrep,"# Selection Strategy is Tournament Selection\n");
		for(i = 0;i < nchrom;i++)
		{
			fprintf(outrep,"# Binary-coded variable no. %d\n", i+1);
			fprintf(outrep,"# No. of bits assigned to it -> %d\n", vlen[i]);
			fprintf(outrep,"# Lower limits on %d-th variable -> %f\n",
				i+1, lim_b[i][0]);
			fprintf(outrep,"# Upper limits on %d-th variable -> %f\n",
				i+1, lim_b[i][1]);
		}
		for(i = 0; i < nvar; i++)
		{
			fprintf(outrep, "# Real-coded variable no. %d\n", i+1);
			fprintf(outrep, "# Lower limits on %dth variable -> %f\n",
				i+1, lim_r[i][0]);
			fprintf(outrep, "# Upper limits on %dth variable -> %f\n",
				i+1, lim_r[i][1]);
			if (ans == 1)
			{
				fprintf(outrep, "# Variable bounds are rigid\n");
			}
			else
			{
				fprintf(outrep, "# Variable bounds are not rigid\n");
			}
		}
		if (nchrom > 0)
		{  
			if (optype == 1)
			{
				fprintf(outrep,
					"# X-over on binary string is SINGLE POINT X-OVER\n");
			}
			if (optype == 2)
			{
				fprintf(outrep,
					"# X-over on binary strings is UNIFORM X-OVER \n");
			}
			if (optype == 3)
			{
				fprintf(outrep,
					"# X-over on binary strings is MULTIPLE X-OVER \n");
			}
		}
		if (nvar > 0)
		{
			fprintf(outrep,
				"# Crossover parameter in the SBX operator is %f\n", di);
		}
		fprintf(outrep,"# Crossover probability -> %f\n", pcross);
		if (nchrom > 0)
		{
			fprintf(outrep,
				"# Mutation probability for binary strings -> %f\n", pmut_b);
		}
		if (nvar > 0)
		{
			fprintf(outrep,
				"# Mutation probability for real-coded vectors -> %f\n", pmut_r);
		}
		fprintf(outrep,"# Control parameter for elitism -> %f\n", relite);
        fprintf(outrep,"# Sharing with respect to distances in fitness and\n"
            "#  paremeter space performed in sequences of -> %d %d\n",
            nsharefit, nsharepar);
		fprintf(outrep,"# Random seed -> %ld\n",seed);
		
		fprintf(outrep,"\n"
			"# POPULATION REPORTS\n"
			"# ------------------");
		
        fflush(NULL);
	}
	
	/*********************************************************************
	*                    Intermediate generations                        *
	*********************************************************************/
	
	if (verb)
	{
		/* Headers for tables */
		fprintf(outrep,"\n"
			"# Generation no. %d\n"
			"# ------------------\n", t);
		if(ncons == 0)
		{
			fprintf(outrep,"# variables (real %d binary %d)  fitness (%d) "
				" rank cublen\n",nvar,nchrom,nfunc);
		}
		else
		{
			fprintf(outrep,"# variables (real %d binary %d)  fitness (%d)"
				" constraint (%d) penalty rank cublen\n",nvar,nchrom,nfunc,ncons);
		}
		
		/* Report of all individuals of current populations */
		ind_ptr = pop_ptr->ind;
		for(i = 0; i < popsize; i++)
		{
			ptr = ind_ptr->xreal;
			for(j = 0;j < nvar;j++)
			{
				fprintf(outrep,"%f ",*ptr++);
			}
			ptr = ind_ptr->xbin;
			for(j = 0;j < nchrom; j++)
			{
				fprintf(outrep,"%f ",*ptr++);
			}  
			ptr = ind_ptr->object;
			for(j = 0;j < nfunc;j++)
			{
				fprintf(outrep,"%f ",*ptr++);
			}  
			if(ncons != 0)
			{
				ptr = ind_ptr->constr;
				for(j = 0;j < ncons;j++)
				{
					fprintf(outrep,"%.2e ",*ptr++);
				}
				fprintf(outrep,"%.2e ",ind_ptr->error);
			}
			fprintf(outrep,"%d ", ind_ptr->rank);
			fprintf(outrep,"%f\n", ind_ptr->cub_len);
			ind_ptr++;
		}
		
        fflush(NULL);
	}
	
	/*********************************************************************
	*           Plot Pareto front at predefined intervalls               *
	*********************************************************************/
	if ( (t*popsize) % (popsize*gener/10) == 0 )
	{
        n = pop_ptr->rankno[0];
        iptr = *pop_ptr->rankar;
		sprintf(name,"pareto%d.out",t);
		pareto = fopen(name, "w");
		fprintf(pareto, "# Pareto front at generation %d (%d of %d)\n",
			t, n, popsize);
		for(j = 0; j < n; j++)
		{
			ind_ptr = pop_ptr->ind + iptr[j];
			if (ind_ptr->error <= 0.0)
			{
				ptr = ind_ptr->object;
				for(i = 0; i < nfunc; i++)
				{
					fprintf(pareto,"%f\t", *ptr++);
				}
				fprintf(pareto,"\n");
			}
		}
		fclose(pareto);
	}
	
	/*********************************************************************
	*                        Final generation                            *
	*********************************************************************/
	if(t == gener)
	{
		/* Print last news and close report file */
		fprintf(outrep,"\n"
			"# No. of crossovers:           %d\n"
			"# No. of mutations:            %d\n"
			"# No. of function evaluations: %d (%.1f%%)\n",
			ncross,nmut,neval,
			100.0 * ((double) neval) / ((double) ((gener+1)*popsize)) );
		fclose(outrep);
		
		/* Open additional output files */
		outfit = fopen("fitness.out","w");
		outvar = fopen("parameters.out","w");
		outrep = fopen("results.out", "w");
		
        n = pop_ptr->rankno[0];
        iptr = pop_ptr->rankar[0];
        for (i = 0; i < n; i++)
        {
            fit[i] = pop_ptr->ind[iptr[i]].fitness[0];
        }
        heapsort2(n, fit, iptr); 
		
		/* Print headers */
		if (ncons)
		{
			fprintf(outfit,"# last generation (%d)\n"
				"# fitness and constraints for feasible and non-dominated\n"
				"# solutions (%d of %d)\n"
				"# fitness vector (%d)  constraints (%d)\n",
				t,n,popsize,nfunc,ncons);
		}
		else
		{
			fprintf(outfit,"# last generation (%d)\n"
				"# fitness for feasible and non-dominated solutions (%d of %d)\n"
				"# fitness vector\n", t,n,popsize);
		}
		fprintf(outvar,"# last generation (%d)\n"
			"# feasible variable vectors for non-dominated solutions "
			"(%d of %d)\n"
			"# real-coded (%d) and binary-coded (%d) variables\n",
			t,n,popsize,nvar,nchrom);
		fprintf(outrep,"# last generation (%d)\n"
			"# all non-dominated solutions (%d of %d)\n", t, n, popsize);
		
		for(j = 0; j < n; j++)
		{
			ind_ptr = pop_ptr->ind + iptr[j];
			if ((ind_ptr->error <= 0.0) && (ind_ptr->rank == 1))
			{
				ptr = ind_ptr->object;
				for(i = 0; i < nfunc; i++)
				{
					f[i] = *ptr;
					fprintf(outfit,"%f\t", *ptr++);
				}
				ptr = ind_ptr->constr;
				for(i = 0; i < ncons; i++)
				{
					cstr[i] = *ptr;
					fprintf(outfit,"%f\t", *ptr++);
				}
				fprintf(outfit,"\n");
				
				ptr = ind_ptr->xreal;
				for(i = 0; i < nvar; i++)
				{
					x[i] = *ptr;
					fprintf(outvar,"%f\t", *ptr++);
				}
				ptr = ind_ptr->xbin;
				for(i = 0; i < nchrom; i++)
				{
					x[nvar+i] = *ptr;
					fprintf(outvar,"%f\t", *ptr++);
				}
				fprintf(outvar,"\n");
				
				app->report(x, nvar+nchrom, f, nfunc, cstr, ncons, (char) !j, outrep);
				
			}  /* feasibility check */
		} /* end of loop j */
		
		/*Closing the files*/
		fclose(outfit);
		fclose(outvar);
		fclose(outrep);
		
	} /* endif (last generation) */
}


/*
results(stream)
driver routine for user defined function app_report;
provides values of fitness, constraints, and variables of
all individuals in current pareto front
*/
void nsga2::print(FILE *stream)
{
	char dbl;
	int i, j, n;
	int *iptr;
	double *ptr;
	individual *ind_ptr, *prev_ptr;
	
	n = pop_ptr->rankno[0];
	iptr = pop_ptr->rankar[0];
	for (i = 0; i < n; i++)
	{
		fit[i] = pop_ptr->ind[iptr[i]].object[0];
	}
	heapsort2(n, fit, iptr); 
	
	for(j = 0; j < n; j++)
	{
		ind_ptr = pop_ptr->ind + iptr[j];
		if (j > 0)
		{
			dbl = 1;
			prev_ptr = pop_ptr->ind + iptr[j-1];
			for (i = 0; i < nvar; i++)
			{
				if (ind_ptr->xreal[i] != prev_ptr->xreal[i])
				{
					dbl = 0;
					break;
				}
			}
			for (i = 0; i < nchrom; i++)
			{
				if (ind_ptr->xbin[i] != prev_ptr->xbin[i])
				{
					dbl = 0;
					break;
				}
			}
		}
		else
		{
			dbl = 0;
		}
		if (dbl == 0)
		{
			ptr = ind_ptr->object;
			for(i = 0; i < nfunc; i++)
			{
				f[i] = *ptr++;
			}
			ptr = ind_ptr->constr;
			for(i = 0; i < ncons; i++)
			{
				cstr[i] = *ptr++;
			}
			ptr = ind_ptr->xreal;
			for(i = 0; i < nvar; i++)
			{
				x[i] = *ptr++;
			}
			ptr = ind_ptr->xbin;
			for(i = 0; i < nchrom; i++)
			{
				x[nvar+i] = *ptr++;
			}
			if (j) // print header only for first individual of current generation
				app->report(x, nvar+nchrom, f, nfunc, cstr, ncons, -1, stream);
			else 
				app->report(x, nvar+nchrom, f, nfunc, cstr, ncons, generi, stream);
		}
	}
}

int nsga2::dump(FILE *stream)
{
	int i, j, m, n, items = 0;
	int *iptr;
	individual *ind_ptr;
	
	items += fwrite(&generi, sizeof(generi), 1, stream);
	items += fwrite(&popsize, sizeof(popsize), 1, stream);
	items += fwrite(&chrom, sizeof(chrom), 1, stream);
	items += fwrite(&nvar, sizeof(nvar), 1, stream);

	m = pop_ptr->maxrank;
	for (i = 0; i < m; i++)
	{
		n = pop_ptr->rankno[i];
		iptr = pop_ptr->rankar[i];
		for (j = 0; j < n; j++)
		{
			fit[j] = pop_ptr->ind[iptr[j]].cub_len;
		}
		heapsort2(n, fit, iptr); 
		
		for(j = 0; j < n; j++)
		{
			ind_ptr = pop_ptr->ind + iptr[j];
			if (chrom)
				items += fwrite(ind_ptr->genes, sizeof(*(ind_ptr->genes)), chrom, stream);
			if (nvar)
				items += fwrite(ind_ptr->xreal, sizeof(*(ind_ptr->xreal)), nvar, stream);
		}
	}
	if (items != 4 + popsize * (chrom + nvar)) 
		items = 0;
	return (items);
}

void nsga2::report(FILE *out, const char *preline)
{
	int i, j;
	/* Print the GA parameters and problem parameters on file pointed at by out */
	fprintf(out,"%sN S G A - 2   P A R A M E T E R S\n",preline);
	fprintf(out,"%sPopulation size: %d\n",preline,popsize);
	if (nchrom > 0)
	{
		fprintf(out,"%sNo. of binary-coded variables: %d\n",preline,nchrom);
		fprintf(out,"%s Total chromosome length: %d bits\n",preline,chrom);
		for(i = 0, j = 1; j < nchrom; j++)
		{
			if (j == nchrom - 1 ||
				vlen[i] != vlen[j] || 
				lim_b[i][0] != lim_b[j][0] ||
				lim_b[i][1] != lim_b[j][1])
			{
				fprintf(out,"%s Binary-coded variables no. %d - %d:\n",preline, i+1, j+1);
				fprintf(out,"%s  Chromosome length: %d bits each variable\n",preline, vlen[i]);
				fprintf(out,"%s  Variable range: [%g,%g]\n",preline,lim_b[i][0],lim_b[i][1]);
				i = j;
			}
		}
		switch (optype)
		{
		case 1:
			fprintf(out,
				"%sCrossover operator for binary string: single point crossover\n",preline);
			break;
		case 2:
			fprintf(out,
				"%sCrossover operator for binary string: uniform crossover"
				" (using %g swap probability)\n",preline,di);
			break;
		case 3:
			fprintf(out,
				"%sCrossover operator for binary string: multiple crossover\n",preline);
		}
		fprintf(out, "%sCrossover probability for binary string: %g\n", preline, pcross);
		fprintf(out, "%sMutation probility for binary string: %g (%g for a single bit)\n",
			preline, pmut_b * double(chrom), pmut_b);
	}
	if (nvar > 0)
	{
		fprintf(out,"%sNo. of real-coded variables: %d\n",preline,nvar);
		for(i = 0, j = 1; j < nvar; j++)
		{
			if (j == nvar - 1 ||
				lim_r[i][0] != lim_r[j][0] ||
				lim_r[i][1] != lim_r[j][1])
			{
				fprintf(out,"%s Real-coded variables no. %d - %d:\n",preline, i+1, j+1);
				fprintf(out,"%s  Variable range: [%g,%g]",preline,lim_r[i][0],lim_r[i][1]);
				if (ans) fprintf(out, ", fixed\n");
				else     fprintf(out, ", at initialization\n");
				i = j;
			}
		}
		fprintf(out, 
			"%sCrossover operator for real-coded variables: simulated binary crossover (SBX)\n",
			preline);
		fprintf(out, "%s Distribution index for SBX: %g\n", preline, di);
		fprintf(out, "%s Crossover probability real-coded vector: %g\n", preline, pcross);
		fprintf(out, "%sMutation operator for real-coded variables: polynomial mutation\n", preline);
		fprintf(out, "%s Distribution index for polynomial mutation operator: %g\n", preline, dim);
		fprintf(out, "%s Mutation probility for real-coded vector: %g (%g for a single variable)\n",
			preline, pmut_r * double(nvar), pmut_r);
	}
	
	fprintf(out,"%sParameter of geometrical probability distribution "
		"used to control elitism: %g\n",preline, relite);
	
	if (nsharefit && !nsharepar)
		fprintf(out,"%sSharing with respect to distances in fitness space.\n", preline);
	else if(!nsharefit && nsharepar)
		fprintf(out,"%sSharing with respect to distances in parameter space.\n", preline);
	else
		fprintf(out,"%sSharing with respect to distances in fitness and\n"
		"%s paremeter space performed in sequences of: %d, %d\n",
		preline,preline,nsharefit, nsharepar);
	
	fprintf(out,"%sNo. of functions: %d\n",preline,nfunc);
   if (nfunc > 1)
   {
      fprintf(out,"%sLinear scaling of objective function values from estimated range of:\n",preline);
      for (i = 0; i < nfunc; i++)
      {
         fprintf(out, "%s f[%d]: %g to %g\n", preline, i, lim_f[i][0], lim_f[i][1]);
      }
      fprintf(out,"%sLinear combination of objective function values using matrix:\n",preline);
      for (i = 0; i < nfunc; i++)
      {
         fprintf(out,"%s %g",preline,T[i][0]);
         for (j = 1; j < nfunc; j++)
            fprintf(out,", %g",T[i][j]);
         fprintf(out,"\n");
      }
   }
	fprintf(out,"%sNo. of constraints: %d\n",preline,ncons);
				
	fflush(out);
}

