%DCMATLAB data section
%
%File reading
% check2dfile   - check file type of 2d data
% check3dfile   - check file type of 3d data
% read2dfile    - read 2d data file (checks type)
% read2drawfile - reads 2d raw (column based) data
% read3dfile    - read 3d data file (checks type)
% read3drawfile - read 3d raw (column based) data
% read4ch       - read 4-channel ABEM configuration file
% readampfile   - read multi-purpose ABEM data file
% readcollect   - read collect file (potential matrix)
% readgemfile   - read sensinv3d data file
% readinv2dfile - read standard 2d data file
% readinv3dfile - read standard 3d data file
% readpro       - read 3d data by use of profiles
% readres2dinvfile - read res2dinv data file
% readres3dinvfile - read res3dinv data file
% readsnd2d     - read 2d data by use of sounding files
% readresecsfile - read resecs ascii output file
%
%File writing
% exportres2dinv - export res2dinv data file
% saveinv2dfile - save standard 2d data file
% saveinv3dfile - save standard 3d data file
% savepro       - save profile data
% saveres2dinvfile - save res2dinv data file
% saveres3dinvfile - save res3dinv data file
%
%Conversion and calculation
% abmn2n        - convert abmn list into data struct
% collectrhoa   - collect apparent resistivities out of potential matrix
% combdata      - combine two data structs (2d or 3d)
% combdata2d    - combine two 2d data structs
% combdata3d    - combine two 3d data structs
% createsxhdata - create surface-hole and cross-hole data
% deldeadelecs  - delete dead (unused) electrodes in data struct
% electrodeerr  - estimate electrode positioning error (2d)
% electrodeerr3d - estimate electrode positioning error (3d) 
% estimateerror - error estimation
% getkonf       - compute configuration factor (2d or 3d)
% getkonf2d     - compute 2d configuration factor
% getkonf3d     - compute 2d configuration factor 
% getpseudos    - extract pseudosections out of 3d data struct
% linesearch    - line search procedure
% mbm2xz        - converts tape measure coodinates into x,z
%
%Plotting
% plotprofiles  - show 3d data set as profiles (subplots)
% midkonf2d     - determine 2d data plot parameter (needed by showdata2d)
% showdata2d    - show 2d data pseudosection
% showdata3d    - show 3d data (pseudosections if available)
% showel        - show electrodes

%deprecated versions
% geomerr       - geometrical error
% linesearch1   - line search
% massbandkorr  - correct Tape measure
% read3draw.m   - read 3d raw (column based) data
% readfeinv3dfile - old version of inv3d file
% topowickel    - old method of tape measure correction