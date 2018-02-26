# Add source files here
EXECUTABLE	:= hello

# Cuda source files (compiled with cudacc)
CUFILES		:= hello.cu

# C/C++ source files (compiled with gcc / c++)
#CCFILES		:= matrixmul_gold.cpp

CUDEPS		:= hello_kernel.cu;


hellomake: hello.cu hello_kernel.cu; nvcc -arch=sm_20 -o $(EXECUTABLE) $(CUFILES) $(CUDEPS) 

