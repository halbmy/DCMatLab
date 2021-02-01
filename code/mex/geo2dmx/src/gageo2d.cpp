#include <mpi.h>
#include <string.h>
#ifdef WIN32
#include <Windows.h>
#else
#include <sys/resource.h>
#endif
#include "nsga2.h"

void signalhandler();

int main(int argc, char **argv)
{
/* Initialize MPI - all command line arguments will be sent to 
    each process */
    MPI_Init(&argc, &argv);
    
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    /* Run programm as a background process with low priority */
#ifdef WIN32
    SetPriorityClass(GetCurrentProcess(), IDLE_PRIORITY_CLASS);
#else
    setpriority(PRIO_PROCESS,0,PRIO_MAX);
#endif
    /* Establish signal handler */
    signalhandler();

    if (argc != 16 && argc != 15)
    {
        if (myid == 0)
        {
            cerr << "Usage: " << argv[0] << " <nind> <igen> <gap> <ngen> "
                << "<vmin> <vmax> <nfuncshare> <nparshare> <datafile> <modelfile> "
                << "<obj1> <obj2> <constrain> <split>"
                << endl;
            cerr << "       " << argv[0] << " <nind> <igen> <gap> <ngen> "
                << "<vmin> <vmax> <nfuncshare> <nparshare> <datafile> <modelfile> "
                << "<obj1> <obj2> <constrain> <split> <nbits>"
                << endl;
        }
    }
    else
    {
		bool binary = ((argc == 16) ? 1 : 0); /* 0 ... real encoding; 1 ... binary encoding */
		int nind = atoi(argv[1]);
		int igap = atoi(argv[2]);
		int gap = atoi(argv[3]);
		int ngen = atoi(argv[4]);
		long is = -myid * nind - 1;  /* different seed for each subpopulation */
		ran1(&is);
		int ox = 2;
		double px = 0.5 + 0.3 * ran1(&is); /* (0.5 ... 0.8) */
		double dx = (binary ? 0.5 : 20.0);
		double pm = (1.0 - px) * ran1(&is);
		double dm = 50.0;
		double r = 0.5;
		int nreal = 0;
		int nbin = 0;
		int bits = (binary ? atoi(argv[15]) : 0);
		double vmin = atof(argv[5]);
		double vmax = atof(argv[6]);
		int vfix = 1;
		int nf = 2;
		int nc = (atoi(argv[13])) ? (1) : (0);
		int nsf = atoi(argv[7]);
		int nsp = atoi(argv[8]);
		char *dataname = argv[9];
		char *modelname = argv[10];
		double obj1 = atof(argv[11]);
		int obj2 = atoi(argv[12]);
      const bool split = (atoi(argv[14]) != 0);
		
		
		/* Get application ready */
		application *app = new application;
		if (binary)
			nbin = app->setup(dataname, modelname, obj1, obj2, vmin, vmax);
		else
			nreal = app->setup(dataname, modelname, obj1, obj2, vmin, vmax);
		
		/* Create subpopulation no. myid */
		nsga2 *ga;
		ga = new nsga2(nind, ngen, ox, px, dx, pm, dm, r,
			nreal, nbin, bits, vmin, vmax, vfix, nf, nc, nsf, nsp, is, app);
		
		/* Evolutionary run, synchronous mode */
		ga->runs(igap,gap,split);
		
		
		/* Destroy GA-object */
		delete ga;
		delete app;
    }
    
    /* Quit MPI */
    MPI_Finalize();
    return(0);
}

