/*
 * read2ddata.cpp
 *
 * void read2ddata(char *fname, int *NoSources, double **XSources, double **ZSources,
 *                 int *NoData, int ***ArSources, double **DataVal, double **DataStd);
 *
 * reads data from a 2D-resistivity survey collected in a file  given by string 'fname'
 *
 * Coordinate system: 0-------> x (profile direction; perpendicular to strike 
 *                    |            of conductivity structures)
 *                    |
 *                    |
 *                    |
 *                    V
 *                    z (positively downward; earth surface at z = 0)
 *
 * Example file:
 * 7                    no. of electrodes used (excluding electrodes at infinity)
 * -30 0                (x,z)-coordinates of 1st electrode [m]
 * -20 0                (x,z)-coordinates of 2nd electrode [m]
 * -10 0                                      :
 * 0 0                                        :
 * 10 0                                       :
 * 20 0                                       :
 * 30 0                 (x,z)-coordinates of 7th electrode [m]
 * 5                    no. of data points
 * 1 2 3 0 100.0 1.0    1st datum: indices of current/potential electrodes (4 entries;
 * 2 3 4 0 99.9 1.1     2nd datum:    0 indicates electrode at infinity); value and standard
 * 3 4 5 0 100.2 0.9     :            deviation of measured apparent resistivity [Ohm*m]
 * 4 5 6 0 100.1 1.1     :
 * 5 6 7 0 99.8 1.2     5th datum:
 *
 *
 * C. Schwarzbach, 04-07-2003
 * schwarzb@geophysik.tu-freiberg.de
 *
 */

#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>

/* void read2ddata(char *fname) */
void read2ddata(char *fname, int *NoSources, double **XSources, double **ZSources,
                int *NoData, int ***ArSources, double **DataVal, double **DataStd)
{
  int i, ns, nd;
  int *is1, *is2, *ip1, *ip2, **ar;
  double *x, *z, *val, *std;

  /* open input file */
  ifstream stream(fname, ios::in|ios::nocreate);
  if (!stream.good())
  {
    cerr << "Couldn't open datafile " << fname << endl;
    exit(0);
  }

  /* get no. of source positions and allocate arrays for source coordinates */
  stream >> ns;
  x = new double [ns];
  if (x == NULL)
  {
    cerr << "Error allocating memory for XSources." << endl;
    exit(0);
  }
  z = new double [ns];
  if (z == NULL)
  {
    cerr << "Error allocating memory for ZSources." << endl;
    exit(0);
  }
  *NoSources = ns;
  *XSources = x;
  *ZSources = z;

  /* read coordinates for sources */
  for (i = 0; i < ns; i++)
  {
    stream >> *x++ >> *z++;
  }

  /* read no. of measurements and allocate memory for electrode configurations and data */
  stream >> nd;
  ar = new int *[4];
  if (ar == NULL)
  {
    cerr << "Error allocating memory for ArSources." << endl;
    exit(0);
  }
  is1 = new int [nd];
  if (is1 == NULL)
  {
    cerr << "Error allocating memory for ArSources[0]." << endl;
    exit(0);
  }
  is2 = new int [nd];
  if (is2 == NULL)
  {
    cerr << "Error allocating memory for ArSources[1]." << endl;
    exit(0);
  }
  ip1 = new int [nd];
  if (ip1 == NULL)
  {
    cerr << "Error allocating memory for ArSources[2]." << endl;
    exit(0);
  }
  ip2 = new int [nd];
  if (ip2 == NULL)
  {
    cerr << "Error allocating memory for ArSources[3]." << endl;
    exit(0);
  }
  val = new double [nd];
  if (val == NULL)
  {
    cerr << "Error allocating memory for DataVal." << endl;
    exit(0);
  }
  std = new double [nd];
  if (std == NULL)
  {
    cerr << "Error allocating memory for DataStd." << endl;
    exit(0);
  }
  *NoData = nd;
  *ArSources = ar;
  ar[0] = is1;     /* Index of 1st source electrode */
  ar[1] = is2;     /* Index of 2nd source electrode */
  ar[2] = ip1;     /* Index of 1st potential electrode */
  ar[3] = ip2;     /* Index of 2nd potential electrode */
  *DataVal = val;     /* Measured apparent resistivity */
  *DataStd = std;     /* Standard deviation of apparent resistivity */

  /* read measured data */
  for (i = 0; i < nd; i++)
  {
    stream >> *is1 >> *is2 >> *ip1 >> *ip2 >> *val >> *std;
    if (*is1 < 0 || *is1 > ns) *is1 = 0; // electrode at infinity
    if (*is2 < 0 || *is2 > ns) *is2 = 0;
    if (*ip1 < 0 || *ip1 > ns) *ip1 = 0;
    if (*ip2 < 0 || *ip2 > ns) *ip2 = 0;
    is1++; is2++; ip1++; ip2++; val++; std++;
  }

  /* close input file */
  stream.close();

  /* print no. of sources and data points */
  cout << "Data set consists of " << nd << " data points using "
    << ns << " electrode locations." << endl;
}
