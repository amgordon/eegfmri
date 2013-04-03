% This make.m is used under Windows



mex -largeArrayDims -c ../blas/daxpy.c 

mex -largeArrayDims -c ../blas/ddot.c 

mex -largeArrayDims -c ../blas/dnrm2.c 

mex -largeArrayDims -c ../blas/dscal.c
mex -largeArrayDims -c ../linear.cpp

mex -largeArrayDims -c ../tron.cpp

mex -largeArrayDims -c ./linear_model_matlab.c -I../

 

mex -largeArrayDims train.c -I../ tron.o linear.o linear_model_matlab.o daxpy.o ddot.o dnrm2.o dscal.o

mex -largeArrayDims predict.c -I../ tron.o linear.o linear_model_matlab.o daxpy.o ddot.o dnrm2.o dscal.o


mex -largeArrayDims libsvmread.c

mex -largeArrayDims libsvmwrite.c