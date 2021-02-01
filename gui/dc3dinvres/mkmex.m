mex -I../src/taucs/src -I../src/taucs/build/linux -L../lib -lg2c -lblas -lcblas -lf77blas -llapack -lmetis -ltaucs fd3dmea.cpp %obj/taucs_ccs*.o obj/taucs_vec*.o obj/readhb.o
