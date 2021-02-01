% DCMATLAB::base functions
% b2r           - Blue-white-red colormap for difference plots
% cbar          - draw colorbar
% chi2          - calculate weighted L2 residual
% cleandirrec   - clean directory recursively
% convadd       - convert hardware address
% curvature     - curvature of a parametric function
% destrip	- strips string from comments (#...)
% epsprint      - export figure into eps (and pdf) file
% equify        - restrict an array to definite values
% exportfig     - general export function for figures
% exportpng     - export figure as png image
% filfak        - calculate and plot filter factors
% findmindist   - find minimum distance from point to point list
% gaulag        - Gauss-Laguerre quadreature points/weights
% gauleg        - Gauss-Legendre quadreature points/weights
% getcindex     - get color indices for graphics patching
% getgrfile     - UI function for exporting graphics file (eps/png)
% gps2xyz       - converts GPS coordinates (HHMMSS) to xyz values
% harmfit       - least squares fit of time series by harmonic functions
% iconify       - use icon file for figure  
% lkurv         - L-curve criterion for optimized regularization
% ilogtrans	- inverse log transformation (+upper/lower bound)
% int2hum	- converts number in human readable (12.3M) string
% interperc	- computes interpercentiles from vector histogram
% inter2dpoint  - interpolate 2d to points by delaunay interpolation
% irls          - iteratively reweighted least squares (L1-norm) weights
% irlsgen       - generalized IRLS (Lp-norm) weights
% itantrans	- inverse tangens transformation with upper/lower bound
% loadprimpot/loadsensmat - load matrix from binary file
% loadsurfergrid- load surfer GRD file into matrix
% logdroptrans  - two-sided logarithmic transformation (with droptolerance)
% logtrans	- logarithmic transformation (+lower/upper bound)
% logtransdiff	- logarithmic difference transformation (+lower/upper bound)
% message       - common message function for output window
% minmax	- yields minimum/maximum value [min(x) max(x)]
% mytextscan	- textscan function for matlab version < R14
% num2strcell   - converts numerical array into string cell (for labels)
% putonxaxis	- adds string onto x axis or replaces an existing one
% readconfigfile- reads confic file of (token=value) into structure
% rms           - (relative) root mean squared misfit
% rndig		- round value(s) to number of counting digits
% savesensmat	- save matrix to binary file
% showvolume    - show 3d volume plot
% showxyz	- plot xyz point list 
% snapline	- snap 2d points onto straight line
% spy2		- improved version of spy
% tantrans	- tangens transformation with upper/lower bound
% tantransdiff	- tangens difference transformation with upper/lower bound
% turnmatrix	- compute turn matrix out of three angles
% veltocrim	- converts velocity into watercontent by CRIM formula
% viereck       - draw rectangle into existing plot
% vol           - volume plotting
% writeconfigfile - writes data structure into config file (token=value)
% writepng      - write figure into png image
% writesurferfile - write matrix into surfer (GRD) file
% zeil2vec	- converts line into vector of numerical values

%deprecated functions for backward compatibility
% antiindex - overall index out of I,J,K position in a 3d grid
% cbar2.m   - alternative version of cbar
% finex1,finex2,finex3 - specialized export functions
% finexbw.m - export function in grayscale
% kruemm    - alternative version of curvature
% tomkurv   - alternative version of curvature
% message3  - alternative version of message 
% rms2.m    - alternative version of rms
