rm -rf ./*.o
rm -rf ./*.mod
module load intel/2017u5
ifort -c Vartype.f90
ifort -c Funcs.f90
ifort -c gidpost.F90
ifort -c StressTrans.f90
ifort -o StressTrans *.o -L./LibL libgidpost.a -lhdf5 -lz -lm

