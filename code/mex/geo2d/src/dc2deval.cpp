#include "dc2d.h"
#include <stdio.h>
#include <iostream.h>


int main(int argc, char **argv)
{
   if (argc != 6)
   {
      cerr << "Usage: " << argv[0]
         << " <datafile> <modelfile> <resultfile> <obj1> <obj2>" << endl;
      return 0;
   }
   
   int i;
   int nf = 2;
   int nc = 1;
   int nx;
   char *dataname = argv[1];
   char *modelname = argv[2];
   char *resultname = argv[3];
   double obj1 = atof(argv[4]);
   int obj2 = atoi(argv[5]);
   double *f, *c, *x, *rho;
   ifstream resultfile;
   
   /* Get application ready */
   application *app = new application;
   nx = app->setup(dataname, modelname, obj1, obj2, 0.0, 0.0);
   f = new double [nf];
   c = new double [nc];
   x = new double [nx];
   rho = new double [nx];

   resultfile.open(resultname, ios::in|ios::nocreate);
   if (!resultfile.good())
   {
      cerr << "Error opening file " << resultname << "." << endl;
      delete app;
      return 0;
   }
   while (resultfile.peek() == '#')
   {
      while(resultfile.get() != '\n');
   }
   for (i = 0; i < nf; i++)
      resultfile >> f[i];
   for (i = 0; i < nc; i++)
      resultfile >> c[i];
   for (i = 0; i < nx; i++)
      resultfile >> rho[i];
   resultfile.close();

   cout  << "parameter no. | rho (Ohm*m)" << endl;
   for (i = 0; i < nx; i++)
      cout << i+1 << " | " << rho[i] << endl;
   cout << "f(rho) = (" << f[0];
   for (i = 1; i < nf; i++)
      cout << "," << f[i];
   cout << ")" << endl
      << "reciprocity error <= " << c[0] << endl;

   while (true)
   {
      cout << "Change parameter no. (0 will quit): " << flush;
      cin >> i;
      if (i <= 0 || i > nx) break;
      cout << "New resistivity value (Ohm*m): " << flush;
      cin >> rho[i-1];
   }

   for (i = 0; i < nx; i++)
      x[i] = -log10(rho[i]);

   app->eval(x, nx, f, nf, c, nc);

   cout << "f(rho) = (" << f[0];
   for (i = 1; i < nf; i++)
      cout << "," << f[i];
   cout << ")" << endl
      << "reciprocity error <= " << c[0] << endl;

   delete [] f;
   delete [] c;
   delete [] x;
   delete [] rho;
   delete app;

   return 0;
}