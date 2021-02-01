/*
 * main.cpp
 *
 * Example main file to demonstrate usage of class geo2d.
 * Reads 2D-model and electrode configurations from files.
 * Calculates apparent resistivity using class methods provided in
 * geo2d.dll (add geo2d.lib to library-modules for linking stage). 
 * Apparent resistivities are written to data file overwriting
 * existing input file.
 *
 * Usage: geo2d <modelfile> <datafile>
 *
 *
 * C. Schwarzbach, 04-07-2003
 * schwarzb@geophysik.tu-freiberg.de
 *
 */

#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>
#include <time.h>

#include "geo2d.h"

void read2dmodel(char *fname, int *NXGrid, double **XGrid, 
                 int *NZGrid, double **ZGrid, double **Conductivity);
void read2ddata(char *fname, int *NoSources, double **XSources, double **ZSources,
                int *NoData, int ***ArSources, double **DataVal, double **DataStd);
void signalhandler();


int main(int argc, char *argv[])
{
   clock_t t_start = clock(), t_end, t_delta;
   int i;
   char *ModelFile, *DataFile;
   int NXGrid, NZGrid;
   double *XGrid, *ZGrid, *Conductivity;
   int NoSources, NoData;
   int **ArSources;
   double *XSources, *ZSources, *DataVal, *DataStd;
   double err;
   
   /* Establish signal handler */
   signalhandler();

   /* Initialize variables */
   NXGrid = 0;
   NZGrid = 0;
   XGrid = NULL;
   ZGrid = NULL;
   Conductivity = NULL;
   NoSources = 0;
   XSources = NULL;
   ZSources = NULL;
   ArSources = NULL;
   DataVal = NULL;
   DataStd = NULL;
   
   /* Check command line arguments and get input file names */
   if (argc == 3)
   {
      ModelFile = argv[1];
      DataFile = argv[2];
   }
   else
   {
      cerr << "usage: " << argv[0] << " <modelfile> <datafile>" << endl;
      exit(0);
   }
   
   /* Read model parameters */
   read2dmodel(ModelFile, &NXGrid, &XGrid, &NZGrid, &ZGrid, &Conductivity);
   
   /* Read data */
   read2ddata(DataFile, &NoSources, &XSources, &ZSources,
      &NoData, &ArSources, &DataVal, &DataStd);
   
   /* Create object geo2d, standard constructor. */
   geo2d obj;

   /* Set up model geometry and conductivities */
   obj.setmodel(NXGrid, XGrid, NZGrid, ZGrid, Conductivity);
   
   /* Set up source positions */
   obj.setsources(NoSources, XSources, ZSources, NoData, ArSources);
   
   /* Calculate apparent resistivity (FD-Code) */
   err = obj.computerhoa(DataVal);

   /* Print apparent resistivities to output file: overwrite datafile */
   ofstream outfile(DataFile,ios::out);
   if (!outfile.good())
   {
      cerr << "Couldn't open file "
         << DataFile
         << " for output." << endl;
      exit(0);
   }
   outfile << NoSources << endl;
   for (i = 0; i < NoSources; i++)
   {
      outfile << XSources[i] << " " << ZSources[i] << endl;
   }
   outfile << NoData << endl;
   outfile.precision(12);
   for (i = 0; i < NoData; i++)
   {
      outfile << ArSources[0][i] << " "
         << ArSources[1][i] << " "
         << ArSources[2][i] << " "
         << ArSources[3][i] << " "
         << DataVal[i] << " "
         << 0 << endl;
   }
   outfile.close();
   
   cout << "Synthetic data written to file " << DataFile << endl
      << "Reciprocity error <= " << err*100.0 << "%" << endl;

   /* Free allocated memory */
   delete [] XSources;
   delete [] ZSources;
   delete [] XGrid;
   delete [] ZGrid;
   delete [] Conductivity;
   delete [] ArSources[0];
   delete [] ArSources[1];
   delete [] ArSources[2];
   delete [] ArSources[3];
   delete [] ArSources;
   delete [] DataVal;
   delete [] DataStd;
   
   t_end = clock();
   t_delta = (clock_t(1000)*(t_end-t_start)) / CLOCKS_PER_SEC;
   cout << "Time elapsed: ";
   if (t_delta > 60*60*1000)
   {
      cout << t_delta/(60*60*1000) << " h, ";
      t_delta %= 60;
   }
   if (t_delta > 60*1000)
   {
      cout << t_delta/(60*1000) << " min, ";
      t_delta %= 60;
   }
   cout << t_delta/1000 << "." << t_delta%1000 << " s." << endl;

   return(0);
}

