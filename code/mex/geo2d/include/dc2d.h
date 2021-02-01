/* evaluation for 2d inversion */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "geo2d.h"

void read2ddata(char *fname, int *NoSources, double **XSources, double **ZSources,
                int *NoData, int ***ArSources, double **DataVal, double **DataStd);
void read2dmodel(char *fname, int *NXGrid, double **XGrid,
                 int *NZGrid, double **ZGrid, int **IConductivity);
void read2dmodel(char *fname, int *NXGrid, double **XGrid,
                 int *NZGrid, double **ZGrid, double **Conductivity);

class application
{
private:
   geo2d *geo2dobj;
   struct array
   {
      int n;
      int *cnd;
      int *nnb;
      int **inb;
      int *flag;
      int *list;
   } block;
   double obj1;
   double C_min, C_max;
   int NoData;
   int NXGrid, NZGrid;
   int obj2;
   double *DataVal, *DataStd, *DataSyn;
   double *XGrid, *ZGrid, *Conductivity;
   double *CXM, *CXP, *CZM, *CZP;
   int *IConductivity;

   void findblocks(int nx, int nz, int *C);
   void datamisfit(double *f, double *e);
   void modelrestriction(double *f);
public:
   application();
   /* Everything that needs to be done before ga is run */
   int setup(char *dataname, char *modelname, double key1, int key2, double xmin, double xmax);
   /* Set estimated range of objective functions */
   void setrange(double **range, int nf);
   /* Given vector of decoded values x, calculate function evaluations and
   constraints and put them to vectors f and c */
   void eval(double *x, int nx, double *f, int nf, double *c, int nc);
   void eval(char *filename, double *f, int nf, double *c, int nc);
   /* Make your own report with given vectors x, f, c as above and print it to stream file */
   void report(double *x, int nx, double *f, int nf, double *c, int nc, int gen, FILE *file);
   /* Everything that needs to be done after ga is finished */
   ~application();
};

