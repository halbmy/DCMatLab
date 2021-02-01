/*
 * geo2dmx.cpp
 *
 * mex-gateway function for geo2d, the finite difference modelling code
 * for the potential of dc-sources on 2D conductivity distributions in 
 * the earth by C. Schwarzbach and R.-U. Börner.
 *
 */

#include "geo2d.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int NXGrid, NZGrid, NoSources, NoReceivers;
	double *XGrid, *ZGrid, *Conductivity,
		*XSources, *ZSources, *XReceivers, *ZReceivers;
  double *Potential;

	/* BEGIN OF MATLAB SPECIFIC PRE-STUFF */

  /* Check number of input and output arguments */
  if(nlhs != 1)
    mexErrMsgTxt("Improper number of output arguments.");
  if(nrhs != 7)
    mexErrMsgTxt("Improper number of input arguments.");

  /* Check type and size of input arguments */
  if(!mxIsDouble(prhs[0]))
	{
    mexErrMsgTxt("First argument 'XGrid' must be a double vector.");
	}
	NXGrid = mxGetM(prhs[0]) * mxGetN(prhs[0]);
  
	if(!mxIsDouble(prhs[1]))
  {
		mexErrMsgTxt("Second argument 'ZGrid' must be a double vector.");
	}
	NZGrid = mxGetM(prhs[1]) * mxGetN(prhs[1]);

  if(!mxIsDouble(prhs[2]) || 
		 mxGetM(prhs[2]) != NZGrid-1 || mxGetN(prhs[2]) != NXGrid-1)
  {
		mexErrMsgTxt("Third argument 'XGrid' must be a double vector "
			"of size (NZGrid-1)x(NXGrid-1).");
	}
  if(!mxIsDouble(prhs[3]))
	{
    mexErrMsgTxt("Fourth argument 'XSources' must be a double vector.");
	}
  NoSources = mxGetM(prhs[3]) * mxGetN(prhs[3]);
  
	if(!mxIsDouble(prhs[4]) || 
		 mxGetM(prhs[4]) * mxGetN(prhs[4]) != NoSources)
	{
    mexErrMsgTxt("Fifth argument 'ZSources' must be a double vector "
			"of the same length as 'XSources'.");
	}
  if(!mxIsDouble(prhs[5]))
	{
    mexErrMsgTxt("Sixth argument 'XReceivers' must be a double vector.");
	}
  NoReceivers = mxGetM(prhs[5]) * mxGetN(prhs[5]);
  
	if(!mxIsDouble(prhs[6]) || 
		 mxGetM(prhs[6]) * mxGetN(prhs[6]) != NoReceivers)
	{
    mexErrMsgTxt("Seventh argument 'ZReceivers' must be a double vector "
			"of the same length as 'XReceivers'.");
	}

  /* Put input arguments to local variables and pointers */
  XGrid        = (double *) mxGetPr(prhs[0]);
  ZGrid        = (double *) mxGetPr(prhs[1]);
  Conductivity = (double *) mxGetPr(prhs[2]);
  XSources     = (double *) mxGetPr(prhs[3]);
  ZSources     = (double *) mxGetPr(prhs[4]);
  XReceivers   = (double *) mxGetPr(prhs[5]);
  ZReceivers   = (double *) mxGetPr(prhs[6]);

  /* Create array for output of phi to MATLAB */
  plhs[0] = mxCreateDoubleMatrix(NoReceivers, NoSources, mxREAL);
	Potential = (double *) mxGetPr(plhs[0]);

	/* END OF MATLAB SPECIFIC PRE-STUFF */

	/* Create object geo2d, standard constructor. */
  geo2d obj;

  /* Set up model geometry and conductivities */
  obj.setmodel(NXGrid, XGrid, NZGrid, ZGrid, Conductivity);
  
	/* Set up source positions */
  obj.setsources(NoSources, XSources, ZSources,
		NoReceivers, XReceivers, ZReceivers);
	
	/* Calculate potential (FD-Code) */
	obj.computephi(Potential);
}