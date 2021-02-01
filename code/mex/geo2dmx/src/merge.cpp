#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int readline(ifstream *infile, double **x);


int main(int argc, char *argv[])
{
   double x1, x2, z1, z2;
   int i, j, n, nlayers, nparams;
   char *ModelFile, *ResultFile, *OutputFile;
   double *x, *params;
      
   /* Check command line arguments and get input file names */
   if (argc == 3)
   {
      ModelFile = argv[1];
      ResultFile = argv[2];
   }
   else
   {
      cerr << "usage: " << argv[0] << " <modelfile> <resultfile>" << endl;
      exit(0);
   }

   /* Read results and print merged modelfiles */

   /* Open result file */
   ifstream result(ResultFile,ios::in|ios::nocreate);
   if (!result.good())
   {
      cerr << "Couldn't open resultfile " << ResultFile << endl;
      exit(0);
   }

   /* Open model file */
   ifstream model(ModelFile,ios::in|ios::nocreate);
   if (!model.good())
   {
      cerr << "Couldn't open modelfile " << ModelFile << endl;
      exit(0);
   }

   /* Setup for output file */
   n = strlen(ResultFile);
   for (i = n-1; i >= 0; i--)
   {
      if (ResultFile[i] == '.')
      {
         ResultFile[i] = '\0';
         break;
      }
   }
   OutputFile = new char [n+9];

   for (i = 1; !result.eof(); i++)
   {
      /* Skip comments */
      while (result.peek() == '#')
         while(result.get() != '\n');
      /* Read one record */
      nparams = readline(&result, &params);
      nparams -= 3; // 2 objective function values and 1 constraint
      
      /* Open output file */
      sprintf(OutputFile, "%s_%.3d.out", ResultFile, i);
      ofstream output(OutputFile, ios::out);
      if (!output.good())
      {
         cerr << "Couldn't open outputfile " << OutputFile << endl;
         exit(0);
      }

      /* Read and print */
      model >> nlayers;
      output << nlayers << endl;

      /* x-grid */
      model.ignore(); // skip eol
      n = readline(&model, &x);
      output << x[0];
      for (j = 1; j < n; j++)
         output << " " << x[j];
      output << endl;
      delete [] x;

      /* z-grid */
      n = readline(&model, &x);
      output << x[0];
      for (j = 1; j < n; j++)
         output << " " << x[j];
      output << endl;
      delete [] x;

      /* layers */
      if (nlayers)
      {
         for (j = 0; j < nlayers; j++)
         {
            model >> n >> z2;
            if (n > nparams)
            {
               cerr << "index into parameter vector exceeds its length." << endl;
               break;
            }
            if (j) output << " ";
            output << params[2+n] << " " << z2;
         }
         output << endl;
      }

      /* blocks */
      for (j = 1; j && !model.eof(); ) // Skip white space
      {
         switch (model.peek())
         {
         case EOF: // need attempt to read beyond EOF to set ios::eofbit true
         case ' ':
         case '\t':
         case '\n':
            model.get();
            break;
         default:
            j = 0;
            break;
         }
      }
      while (!model.eof())
      {
         model >> x1 >> x2 >> z1 >> z2 >> n;
         if (n > nparams)
         {
            cerr << "index into parameter vector exceeds its length." << endl;
            break;
         }
         output << x1 << " " << x2 << " " << z1 << " " << z2 << " " << params[2+n] << endl;
         for (j = 1; j && !model.eof(); ) // Skip white space
         {
            switch (model.peek())
            {
            case EOF: // need attempt to read beyond EOF to set ios::eofbit true
            case ' ':
            case '\t':
            case '\n':
               model.get();
               break;
            default:
               j = 0;
               break;
            }
         }
      }

      delete [] params;
      output.close();
      model.clear();
      model.seekg(0);
      for (j = 1; j && !result.eof(); ) // Skip white space
      {
         switch (result.peek())
         {
         case EOF: // need attempt to read beyond EOF to set ios::eofbit true
         case ' ':
         case '\t':
         case '\n':
            result.get();
            break;
         default:
            j = 0;
            break;
         }
      }
   }
   result.close();
   model.close();

   delete [] OutputFile;

   return 0;
}