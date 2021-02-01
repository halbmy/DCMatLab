#include "dc2d.h"
#include <iostream.h>


int main(int argc, char **argv)
{
   if (argc != 6)
   {
      cerr << "Usage: " << argv[0]
         << " <datafile> <modelfile> <resultfile> <obj1> <obj2>" << endl;
      return 0;
   }
   

   char *dataname = argv[1];
   char *modelname = argv[2];
   char *resultname = argv[3];
   double obj1 = atof(argv[4]);
   int obj2 = atoi(argv[5]);
   double f[2], c[1];
   
   /* Get application ready */
   application *app = new application;
   app->setup(dataname, modelname, obj1, obj2, 0.0, 0.0);
   app->eval(resultname, f, 2, c, 1);

   cout << "data misfit = " << f[0] << endl
        << "model restriction = " << f[1] << endl
        << "reciprocity error <= " << c[0] << endl;

   delete app;

   return 0;
}