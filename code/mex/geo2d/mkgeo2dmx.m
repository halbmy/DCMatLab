outdir=['..' filesep '..' filesep 'dc2dinvres' filesep 'geo2dmx'];
% mex([' -o ' outdir ' -I./include src/math/bessel.cpp src/geo2d/computephi.cpp src/geo2d/coupling.cpp src/math/gammln.cpp src/math/gaulag.cpp src/math/gauleg.cpp src/geo2d/geo2d.cpp src/geo2dmx.cpp src/geo2d/grid.cpp src/math/hygf.cpp src/geo2d/integration.cpp src/geo2d/locate.cpp src/geo2d/model2d.cpp src/geo2d/potential.cpp src/math/pythag.cpp src/geo2d/readconfig.cpp src/geo2d/righthandside.cpp src/geo2d/solvelineq.cpp src/geo2d/sourcepos.cpp src/geo2d/superpos.cpp src/geo2d/wavenumber.cpp -llapack']);
mex -o ../../dc2dinvres/geo2dmx -I./include src/math/bessel.cpp src/geo2d/computephi.cpp src/geo2d/coupling.cpp src/math/gammln.cpp src/math/gaulag.cpp src/math/gauleg.cpp src/geo2d/geo2d.cpp src/geo2dmx.cpp src/geo2d/grid.cpp src/math/hygf.cpp src/geo2d/integration.cpp src/geo2d/locate.cpp src/geo2d/model2d.cpp src/geo2d/potential.cpp src/math/pythag.cpp src/geo2d/readconfig.cpp src/geo2d/righthandside.cpp src/geo2d/solvelineq.cpp src/geo2d/sourcepos.cpp src/geo2d/superpos.cpp src/geo2d/wavenumber.cpp -llapack
%src/geo1d/potential1d.cpp 
