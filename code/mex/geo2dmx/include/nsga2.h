/* nsga2.h */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "random.h"
#include "dc2d.h"

#ifndef min
#define min(a,b) (a < b) ? a : b
#endif
#ifndef max
#define max(a,b) (a > b) ? a : b
#endif

/*** STRUCTURE DEFINITIONS ***/

/* Individual structure */
typedef struct
{
  int *genes,      /*binary chromosome*/
    rank,          /*Rank of the individual*/
    change;        /*Flag for evaluation if changed during x-over or mutation*/
  double *xbin,    /*list of decoded value of the chromosome */
    *xreal,        /*list of real variables*/
    *object,       /*Objective function values*/
    *fitness,      /*Fitness values*/
    *constr,       /*Constraints values*/
    cub_len,       /*crowding distance of the individual*/
    error;         /*overall constraint violation for the individual*/
} individual;            


/* Population Structure */
typedef struct
{
  int maxrank;     /*Maximum rank present in the population*/
  int *rankno;     /*No. of individuals at different ranks*/
  int **rankar;    /*Indizes of individuals at different ranks*/
  individual *ind; /*Array of individuals*/
} population;            

/*** CLASS DEFINITION FOR NSGA2 ***/

class nsga2 
{
private:
  /* variables */
  int popsize,      /*Population Size*/
    popsize2,       /*Population Size * 2*/
	generi,         /*Actual no of generation*/
    gener,          /*Max. no of generations*/
    nvar,           /*No of real variables*/
    nchrom,         /*No of binary variables*/
    chrom,          /*Chromosome size*/
    *vlen,          /*Array to store no of bits for each variable*/
    ncons,          /*No of Constraints*/
    nmut,           /*No of Mutations */
    ncross,         /*No of crossovers*/
    ans,            /*Rigid/flexible variable bounds*/
    optype,         /*Cross-over type*/
    nfunc,          /*No of objective functions*/
    *flag,          /*Array of flags used in rank, rankc, and keepalive*/
    *index,         /*Array of indices used in keepaliven and share*/
    verb,           /*Switch for more or less detailed report*/
    nsharefit,      /*Number of generations to share in fitness space*/
    nsharepar,      /*Number of generations to share in parameter space*/
    nshare,         /*Counter for sharing*/
    neval,          /*Number of function evaluations*/
    noldfit,
    conv;
  long seed;        /*Random Seed*/
  double pcross,    /*Cross-over Probability*/
    pmut_b,         /*Binary mutation probability*/
    pmut_r,         /*Real mutation probability*/
	relite,         /*Controlls elitism*/
    **lim_b,        /*Limits of binary variables in array*/
    **lim_r,        /*Limits of real variables in array*/
    di,             /*Distribution Index for the Cross-over*/
    dim,            /*Distribution Index for the Mutation*/
    *x,             /*Vector of decoded parameters for func*/
    *f,             /*Vector objective functions for func*/
    **lim_f,        /*Range of objective function values of current population*/
    **T,            /*Matrices to transform fitness values*/
    *cstr,          /*Vector of constraints for func*/
    *fit,           /*Array of fitness/dummy fitness values used
                      in keepaliven and share*/
    **oldfit;       /*Array for fitness of old pareto front used 
                      in converge*/
  population *pop_ptr;

  /* functions */
  void memory();
  void transf(int k = 1, int n = 1);
  void report(int t);
  void realmutate(individual *ind_ptr);
  void sbcross(individual *ind_ptr);
  int indcmp(double *fit1, double *fit2);
  void rank(int size);
  void rankc(int size);
  void mutate(individual *ind_ptr);
  void keepalive(int size);
  int iround(double x);
  void heapsort(unsigned long n, double *ra);
  void heapsort2(unsigned long n, double *ra, int *ia);
  void share(int size);
  void share2(int size);
  void input();
  void init(individual *ind_ptr);
  void init(individual *ind_ptr, FILE *stream);
  int dump(FILE *stream);
  void decode(individual *ind_ptr);
  void simplecross(individual *ind_ptr);
  void multicross(individual *ind_ptr);
  void unicross(individual *ind_ptr);
  void selection();
  void printtime(double secs);
  void converge();
  void converge(int gen);
  void func(individual *ind_ptr);
  int pack(int *buf, int bufsize, int rank);
  int unpack(int *buf, int bufsize, int pos);
  void report(FILE *out, const char *preline);
  void enlarge(int newpopsize);
  application *app;

public:
  nsga2();
  nsga2(int verbose);
  nsga2(int nind, int ngen, int ox, double px, double dx,
             double pm, double dm, double r, int nreal, int nbin,
             int bits, double vmin, double vmax, int vfix, int nf,
             int nc, int nsf, int nsp, long is, application *appl);
  ~nsga2();
  void run(int verbose);
  void runs(int igap, int gap, const bool split);
  void start(int verbose, FILE *stream = NULL);
  void step(int ngen);
  void statistics(double *mean);
  void print(FILE *stream);
};

