#include <stdio.h>
#include <stdlib.h>
#include <math.h>

unsigned long round(double x)
{
   if (x >= 0.0)
      return (x-double((unsigned long)(x)) >= 0.5) ? (unsigned long)(x)+1L : (unsigned long)(x);
   else
      return (double((unsigned long)(x))-x >= 0.5) ? (unsigned long)(x)-1L : (unsigned long)(x);
}

int main(int argc, char **argv)
{
   if (argc != 5)
   {
      printf("Usage: %s <# of parameters> <# of populations> <population size>"
         " <# of bits>\n",
         argv[0]);
      return 0;
   }
   
   int i, j, bit;
   int nvar = atoi(argv[1]);
   int npop = atoi(argv[2]);
   int popsize = atoi(argv[3]);
   int nbit = atoi(argv[4]);
   int nind = 1 << nbit;
   unsigned long maxbin = (1 << nbit) - 1;
   unsigned long bin, nib;
   int *code = new int [nbit];
   char outfilename[128];
   long pos = 0;
   int length, start;
   double dump = 0.0;
   FILE *outstream;
   
   int ivar[4] = {0,nind,nbit*nvar,0};
   int *dvar = new int [nind*nbit*nvar];
   int *ddvar = dvar;
   
   for (bin = 0; bin <= maxbin; bin++)
   {
      nib = bin;
      for (j = nbit-1; j >= 0; j--) // decimal to binary
      {
         if (bin >= (1UL << j))
         {
            code[nbit-1-j] = 1;
            bin -= (1L << j);
         }
         else
         {
            code[nbit-1-j] = 0;
         }
      }
      for (j = nbit-1; j > 0; j--) // binary to Gray
      {
         code[j] ^= code[j-1];
      }
      
      // Check reverse direction
      bit = 0;
      for (j = nbit-1; j >= 0; j--)
      {
         bit ^= code[nbit-1-j];
         if (bit) bin += (1 << j);
      }
      if (bin != nib) fprintf(stderr, "Error encoding %d -> %d\n", nib, bin);
      
      for (j = 0; j < nvar; j++)
      {
         for (i = 0; i < nbit; i++, dvar++)
         {
            *dvar = code[i];
            printf("%d", *dvar);
         }
         printf("\n");
      }
   }
   
   dvar = ddvar;
   nvar *= nbit;
   for (i = 0; i < npop; i++)
   {
      if (popsize >= nind)
      {
         start = 0;
         length = nind * nvar;
         ivar[1] = nind;
      }
      else
      {
         if (i > 0)
            start = (i * (nind - popsize) * nvar) / (npop - 1);
         else
            start = 0;
         length = popsize * nvar;
         ivar[1] = popsize;
      }
      sprintf(outfilename, "population%d.bin", i+1);
      outstream = fopen(outfilename, "wb");
      fwrite(ivar, sizeof(int), 4, outstream);
      fwrite(dvar+start, sizeof(int), length, outstream);
      fclose(outstream);
      printf("Entries %d - %d written on file %s.\n",
         start/nvar+1, (start+length)/nvar, outfilename);
   }
   
   delete [] dvar;
   delete [] code;
   return 0;
}