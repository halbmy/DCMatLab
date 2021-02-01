#include <math.h>
#define IA 16807 
#define IM 2147483647 
#define AM (1.0/IM) 
#define IQ 127773 
#define IR 2836 
#define NTAB 32 
#define NDIV (1+(IM-1)/NTAB) 
#define EPS 1.0e-16 
#define RNMX (1.0-EPS) 

double ran1(long *idum);
double gasdev(long *idum);

