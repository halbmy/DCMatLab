#include <iostream.h>
#include <fstream.h>
#include "geo2d.h"

void geo2d::readconfiguration(char *configname)
{
   double fnum, f, g;
   int inum, jnum, i, nrx, nrz;
   char cnum;
   char ch;
   //ifstream configfile(configname, ios::in|ios::nocreate);
   ifstream configfile(configname, ios::in);
   if (configfile == NULL)
   {
      //cerr << "Error opening configuration file. Using default values." << endl;
   }
   else
   {
      /* Read configuration from file */
      for (ch = configfile.peek(); ch == '#' || ch == ' ' || ch == '\t' || ch == '\n'; ch = configfile.peek())
         while (configfile.get() != '\n');
      configfile >> inum;
      if (inum >= 0)
      {
         config.NoEquiX = inum;
      }
      else
      {
         cerr << inum << " is an invalid value for NoEquiX. Using default value " << config.NoEquiX << "." << endl;
      }
      
      for (ch = configfile.peek(); ch == '#' || ch == ' ' || ch == '\t' || ch == '\n'; ch = configfile.peek())
         while (configfile.get() != '\n');
      configfile >> inum;
      if (inum >= 0)
      {
         config.NoEquiZ = inum;
      }
      else
      {
         cerr << inum << " is an invalid value for NoEquiZ. Using default value " << config.NoEquiZ << "." << endl;
      }
            
      for (ch = configfile.peek(); ch == '#' || ch == ' ' || ch == '\t' || ch == '\n'; ch = configfile.peek())
         while (configfile.get() != '\n');
      configfile >> inum;
      if (inum >= 0)
      {
         config.NoProlX = inum;
      }
      else
      {
         cerr << inum << " is an invalid value for NoProlX. Using default value " << config.NoProlX << "." << endl;
      }
               
      for (ch = configfile.peek(); ch == '#' || ch == ' ' || ch == '\t' || ch == '\n'; ch = configfile.peek())
         while (configfile.get() != '\n');
      configfile >> inum;
      if (inum >= 0)
      {
         config.NoProlZ = inum;
      }
      else
      {
         cerr << inum << " is an invalid value for NoProlZ. Using default value " << config.NoProlZ << "." << endl;
      }
      
      for (ch = configfile.peek(); ch == '#' || ch == ' ' || ch == '\t' || ch == '\n'; ch = configfile.peek())
         while (configfile.get() != '\n');
      configfile >> fnum;
      if (fnum >= 1.0)
      {
         config.FacProlX = fnum;
      }
      else
      {
         cerr << fnum << " is an invalid value for FacProlX. Using default value " << config.FacProlX << "." << endl;
      }
      
      for (ch = configfile.peek(); ch == '#' || ch == ' ' || ch == '\t' || ch == '\n'; ch = configfile.peek())
         while (configfile.get() != '\n');
      configfile >> fnum;
      if (fnum >= 1.0)
      {
         config.FacProlZ = fnum;
      }
      else
      {
         cerr << fnum << " is an invalid value for FacProlZ. Using default value " << config.FacProlZ << "." << endl;
      }
      
      for (ch = configfile.peek(); ch == '#' || ch == ' ' || ch == '\t' || ch == '\n'; ch = configfile.peek())
         while (configfile.get() != '\n');
      configfile >> cnum;
      switch (cnum)
      {
      case 'a':
      case 'A':
         config.ModProlong = 'a';
         break;
      case 'b':
      case 'B':
         config.ModProlong = 'b';
         break;
      default:
         cerr << cnum << " is an invalid value for ModProlong. Using default value " << config.ModProlong << "." << endl;
         break;
      }
      
      for (ch = configfile.peek(); ch == '#' || ch == ' ' || ch == '\t' || ch == '\n'; ch = configfile.peek())
         while (configfile.get() != '\n');
      configfile >> inum;
      if (inum >= 0)
      {
         config.NoRefineX = inum;
      }
      else
      {
         cerr << inum << " is an invalid value for NoRefineX. Using default value " << config.NoRefineX << "." << endl;
      }
      
      for (ch = configfile.peek(); ch == '#' || ch == ' ' || ch == '\t' || ch == '\n'; ch = configfile.peek())
         while (configfile.get() != '\n');
      configfile >> inum;
      if (inum >= 0)
      {
         config.NoRefineZ = inum;
      }
      else
      {
         cerr << inum << " is an invalid value for NoRefineZ. Using default value " << config.NoRefineZ << "." << endl;
      }
      
      for (ch = configfile.peek(); ch == '#' || ch == ' ' || ch == '\t' || ch == '\n'; ch = configfile.peek())
         while (configfile.get() != '\n');
      configfile >> inum;
      for (ch = configfile.peek(); ch == '#' || ch == ' ' || ch == '\t' || ch == '\n'; ch = configfile.peek())
         while (configfile.get() != '\n');
      configfile >> jnum;
      if (inum >= 0 && jnum >= 0 && inum + jnum > 0)
      {
         config.NoLegendre = inum;
         config.NoLaguerre = jnum;
      }
      else
      {
         cerr << inum << "/" << jnum << " are invalid values for NoLegendre/NoLaguerre. Using default values " 
            << config.NoLegendre << "/" << config.NoLaguerre << "." << endl;
      }
      
      for (ch = configfile.peek(); ch == '#' || ch == ' ' || ch == '\t' || ch == '\n'; ch = configfile.peek())
         while (configfile.get() != '\n');
      configfile >> cnum;
      switch (cnum)
      {
      case 'd':
      case 'D':
         config.Boundary = 'd';
         break;
      case 'm':
      case 'M':
         config.Boundary = 'm';
         break;
      default:
         cerr << cnum << " is an invalid value for Boundary. Using default value " << config.Boundary << "." << endl;
         break;
      }
      
      for (ch = configfile.peek(); ch == '#' || ch == ' ' || ch == '\t' || ch == '\n'; ch = configfile.peek())
         while (configfile.get() != '\n');
      configfile >> cnum;
      switch (cnum)
      {
      case 'a':
      case 'A':
         config.Background = 'a';
      case 'm':
      case 'M':
         config.Background = 'm';
      case 's':
      case 'S':
         config.Background = 's';
         break;
      default:
         cerr << cnum << " is an invalid value for Background. Using default value " << config.Background << "." << endl;
         break;
      }
      
      /* print list of chosen program parameters to screen */
      cout << endl << "* * * * * * * * * * * * * * * * GEO2D "
         << "* * * * * * * * * * * * * * * *" << endl
         << "                          program parameters" << endl;
      cout 
         << "no. of prolongations " << endl
         << " equidistantly in x-direction ... " << config.NoEquiX << endl
         << " equidistantly in z-direction ... " << config.NoEquiZ << endl
         << " increasingly in x-direction .... " << config.NoProlX << endl
         << " increasingly in z-direction .... " << config.NoProlZ << endl
         << "factor for prolongation" << endl 
         << " in x-direction ................. " << config.FacProlX << endl
         << " in z-direction ................. " << config.FacProlZ << endl
         
         << "modus of prolongation ........... ";
      switch (config.ModProlong) 
      {
      case 'a':
         cout << "average value" << endl;
         break;
      case 'b':
         cout << "boundary values" << endl;
         break;
      }
      cout << "no. of refinements " << endl
         << " in x-direction ................. " << config.NoRefineX << endl
         << " in z-direction ................. " << config.NoRefineZ << endl
         << "no. of Gauss-Legendre abscissas . " << config.NoLegendre << endl
         << "no. of Gauss-Laguerre abscissas . " << config.NoLaguerre << endl
         << "boundary condition .............. ";
      switch (config.Boundary) 
      {
      case 'd':
         cout << "Dirichlet" << endl;
         break;
      case 'm':
         cout << "mixed" << endl;
         break;
      }
      cout << "background conductivity ......... ";
      switch (config.Background) 
      {
      case 'a':
         cout << "average conductivity" << endl;
         break;
      case 'm':
         cout << "average conductivity at source location" << endl;
         break;
      case 's':
         cout << "conductivity of source location" << endl;
         break;
      }
   }

	// New number for prolongation depending on number of refinements.
   nrx = config.NoRefineX + 1;
   nrz = config.NoRefineZ + 1;
   if (i = config.NoProlX)
   {
      f = config.FacProlX;
      for (; i > 1; i--)
      {
         f = config.FacProlX * (1.0 + f);
      }
      f *= (double) nrx;
      g = config.FacProlX;
      for (; g < f; i++)
      {
         g = config.FacProlX * (1.0 + g);
      }
   }
   config.NoProlX = i;
   config.NoEquiX *= nrx;

   if (i = config.NoProlZ)
   {
      f = config.FacProlZ;
      for (; i > 1; i--)
      {
         f = config.FacProlZ * (1.0 + f);
      }
      f *= (double) nrz;
      g = config.FacProlZ;
      for (; g < f; i++)
      {
         g = config.FacProlZ * (1.0 + g);
      }
   }
   config.NoProlZ = i;
   config.NoEquiZ *= nrz;
}
