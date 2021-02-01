/*
 * read2dmodel.cpp
 *
 * void read2dmodel(char *fname, int *NXGrid, double **XGrid,
 *                  int *NZGrid, double **ZGrid, double **Conductivity)
 *
 * reads geometry and conductivities for a 2D-resistivity-model from a file 
 * given by string 'fname'
 *
 * Coordinate system: 0-------> x (profile direction; perpendicular to strike 
 *                    |            of conductivity structures)
 *                    |
 *                    |
 *                    |
 *                    V
 *                    z (positively downward; earth surface at z = 0)
 *
 *
 * C. Schwarzbach, 04-07-2003
 * schwarzb@geophysik.tu-freiberg.de
 *
 */

#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>
#include <math.h>

int readline(ifstream *infile, double **x);

void read2dmodel(char *fname, int *NXGrid, double **XGrid,
                 int *NZGrid, double **ZGrid, int **Conductivity)
{
   int i, j, k, m, n, nx, nz, idum;
   double xstart, xend, zstart, zend, dum;
   int *sigma;
   double *x, *z;
   
   /* Open input file */
   ifstream file2d(fname,ios::in|ios::nocreate);   // Open input file
   if (!file2d.good())
   {
      cerr << "Couldn't open modelfile " << fname << endl;
      exit(0);
   }
   
   /* Read number of horizontal layers */ 
   file2d >> n;
   file2d.ignore(); // skip eol

   /* Read start, increment and end of x-grid-coordinates and create x-grid */
   nx = readline(&file2d, &x);
   if (nx == 3)
   {
      xstart = x[0];
      dum = x[1];
      xend = x[2];
      delete x;
      nx = int(ceil((xend-xstart)/dum)) + 1;
      x = new double[nx];
      for (i = 0; i < nx; i++)
      {
         x[i] = xstart;
         xstart += dum;
      }
   }

   /* Read start, increment and end of z-grid-coordinates and create z-grid */
   nz = readline(&file2d, &z);
   if (nz == 3)
   {
      zstart = z[0];
      dum = z[1];
      zend = z[2];
      delete z;
      nz = int(ceil((zend-zstart)/dum)) + 1;
      z = new double[nz];
      for (i = 0; i < nz; i++)
      {
         z[i] = zstart;
         zstart += dum;
      }
   }

   /* Read resistivity and depth of layers */
   sigma = new int [(nx-1)*(nz-1)];
   j = 0;
   for (i = 1;i <= n;i++)
   {
      file2d >> idum >> zend;
      if (i < n)
      {
         while (z[j] < zend && j < nz-1)
         {
            for (k = 0; k < nx-1; k++)
            {
               sigma[k*(nz-1)+j] = idum; // Array Conductivity contains conductivities
            }                           // in columnwise order
            j++;
         }
      }
      else
      {
         while (z[j] <= zend && j < nz-1)
         {
            for (k = 0; k < nx-1; k++)
            {
               sigma[k*(nz-1)+j] = idum;
            }
            j++;
         }
      }
   }
   
   /* Read resistivity and coordinates of blocks */
   for (i = 1; i && !file2d.eof(); ) // Skip white space
   {
      switch (file2d.peek())
      {
      case EOF: // need attempt to read beyond EOF to set ios::eofbit true
      case ' ':
      case '\t':
      case '\n':
         file2d.get();
         break;
      default:
         i = 0;
         break;
      }
   }
   m = 0;
   while (!file2d.eof())
   {
      file2d >> xstart >> xend >> zstart >> zend >> idum;
      for (i = 0; i < nx; i++)
      {
         if (x[i] >= xstart)
         {
            m = i;
            break;
         }
      }
      for (i = 0; i < nz; i++)
      {
         
         if (z[i]>=zstart)
         {
            n = i;
            break;
         }
      }
      for (j = m; x[j] < xend && j < nx-1; j++)
      {
         for (k = n; z[k] < zend && k < nz-1; k++)
         {
            sigma[j*(nz-1)+k] = idum;
         }
      }
      for (i = 1; i && !file2d.eof(); ) // Skip white space
      {
         switch (file2d.peek())
         {
         case EOF: // need attempt to read beyond EOF to set ios::eofbit true
         case ' ':
         case '\t':
         case '\n':
            file2d.get();
            break;
         default:
            i = 0;
            break;
         }
      }
   }

   /* Close input file */
   file2d.close();
   
   /* Assign temporary pointers to input/output variables */
   *NXGrid = nx;
   *XGrid = x;
   *NZGrid = nz;
   *ZGrid = z;
   *Conductivity = sigma;
}

void read2dmodel(char *fname, int *NXGrid, double **XGrid,
                 int *NZGrid, double **ZGrid, double **Conductivity)
{
   int i, j, k, m, n, nx, nz;
   double xstart, xend, zstart, zend, dum;
   double *x, *z, *sigma;
   
   /* Open input file */
   ifstream file2d(fname,ios::in|ios::nocreate);   // Open input file
   if (!file2d.good())
   {
      cerr << "Couldn't open modelfile " << fname << endl;
      exit(0);
   }
   
   /* Read number of horizontal layers */ 
   file2d >> n;
   file2d.ignore(); // skip eol

   /* Read start, increment and end of x-grid-coordinates and create x-grid */
   nx = readline(&file2d, &x);
   if (nx == 3)
   {
      xstart = x[0];
      dum = x[1];
      xend = x[2];
      delete x;
      nx = int(ceil((xend-xstart)/dum)) + 1;
      x = new double[nx];
      for (i = 0; i < nx; i++)
      {
         x[i] = xstart;
         xstart += dum;
      }
   }

   /* Read start, increment and end of z-grid-coordinates and create z-grid */
   nz = readline(&file2d, &z);
   if (nz == 3)
   {
      zstart = z[0];
      dum = z[1];
      zend = z[2];
      delete z;
      nz = int(ceil((zend-zstart)/dum)) + 1;
      z = new double[nz];
      for (i = 0; i < nz; i++)
      {
         z[i] = zstart;
         zstart += dum;
      }
   }

   /* Read resistivity and depth of layers */
   sigma = new double [(nx-1)*(nz-1)];
   j = 0;
   for (i = 1;i <= n;i++)
   {
      file2d >> dum >> zend;
      dum = 1.0 / dum;
      if (i < n)
      {
         while (z[j] < zend && j < nz-1)
         {
            for (k = 0; k < nx-1; k++)
            {
               sigma[k*(nz-1)+j] = dum; // Array Conductivity contains conductivities
            }                           // in columnwise order
            j++;
         }
      }
      else
      {
         while (z[j] <= zend && j < nz-1)
         {
            for (k = 0; k < nx-1; k++)
            {
               sigma[k*(nz-1)+j] = dum;
            }
            j++;
         }
      }
   }
   
   /* Read resistivity and coordinates of blocks */
   for (i = 1; i && !file2d.eof(); ) // Skip white space
   {
      switch (file2d.peek())
      {
      case EOF: // need attempt to read beyond EOF to set ios::eofbit true
      case ' ':
      case '\t':
      case '\n':
         file2d.get();
         break;
      default:
         i = 0;
         break;
      }
   }
   m = 0;
   n = 0;
   while (!file2d.eof())
   {
      file2d >> xstart >> xend >> zstart >> zend >> dum;
      dum = 1.0 / dum;
      for (i = 0; i < nx; i++)
      {
         if (x[i] >= xstart)
         {
            m = i;
            break;
         }
      }
      for (i = 0; i < nz; i++)
      {
         
         if (z[i]>=zstart)
         {
            n = i;
            break;
         }
      }
      for (j = m; x[j] < xend && j < nx-1; j++)
      {
         for (k = n; z[k] < zend && k < nz-1; k++)
         {
            sigma[j*(nz-1)+k] = dum;
         }
      }
      for (i = 1; i && !file2d.eof(); ) // Skip white space
      {
         switch (file2d.peek())
         {
         case EOF: // need attempt to read beyond EOF to set ios::eofbit true
         case ' ':
         case '\t':
         case '\n':
            file2d.get();
            break;
         default:
            i = 0;
            break;
         }
      }
   }

   /* Close input file */
   file2d.close();
   
   /* Assign temporary pointers to input/output variables */
   *NXGrid = nx;
   *XGrid = x;
   *NZGrid = nz;
   *ZGrid = z;
   *Conductivity = sigma;

   /* Check if conductivity array has only valid values */
   n = (nx-1)*(nz-1);
   for (i = 0; i < n; i++)
   {
      if (sigma[i] <= 0.0 || sigma[i] == HUGE_VAL)
      {
         cerr << "Invalid conductivity values. Aborting." << endl;
         exit(0);
      }
   }
}

