compile using gcc 4.4.7 (standard on HPC, don't load any extra modules)

load asegpaw/3.9.1.4567-0.11.0.13004.py26

(only needed if you changed ins/outs in python-fortran interface)
f2py -m scme2015 -h scme2015.pyf main.f --overwrite-signature

make clean if old .o's are there

make -f Makefile_f2py

f2py --fcompiler=gfortran -c scme2015.pyf *.o
