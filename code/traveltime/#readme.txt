dijkstra.m - Dijkstra's algorithm (Shortest path problem)
getstartmodel.m - Function to create gradient starting model
grid2mesh.m - Converts FD grid into FE mesh
gridconst2d.m - Constraint matrix for 2d FD mesh
meshsmoothness - Constraint matrix for unstructured mesh
plotshot.m - Plot RA data (and responses)
readgli.m - read GLI3D data file
readgrm.m - read GREMIX data file
readsensvectors.m - read sensitivity matrix by vectors
showrays.m - show rays starting from one shot position
tripatchmod.m - show mesh attributes (m/s)
waymatrix.m - compute way matrix for traveltime calculation
writepoly.m - write poly file as input for dctriangle
dctriangle(.exe) - create (para) mesh from poly file

ratomo.fig/.m - GUI for refraction tomography
Options to be set:
  Mesh: depth, quality of Mesh, smoothing
  Start: vmin/vmax, Gradient or 2layer-model
  Inversion: Parameter (log,v,s), lower/upper bound, smoothness
             Flatness factor, robust, blocky, ...
