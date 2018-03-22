# Add source files here
EXECUTABLE	:= POSTERR

# Cuda source files (compiled with cudacc)
CUFILES		:= harness.cu

# C/C++ source files (compiled with gcc / c++)
#CCFILES		:= matrixmul_gold.cpp

CUDEPS		:= ODE_kernel.cu;


hellomake: harness.cu ODE_kernel.cu; nvcc -Xptxas -O1 -arch=sm_20 -o $(EXECUTABLE) $(CUFILES) $(CUDEPS) 

clean: 
	rm POSTERR #*.o
