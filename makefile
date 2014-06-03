# This makefile borrows a lot from Eric Aldrich's VFI repo on GitHub. Kudos to him.
# We try to do a separation compilation here because nvcc doesn't fully support c++11 <random>

# Paths for Linux CUDA
ICUDA    = /usr/local/cuda-5.5/include
LCUDA    = /usr/local/cuda-5.5/lib64
ICUDA_MAC = /Developer/NVIDIA/CUDA-5.5/include
LCUDA_MAC = /Developer/NVIDIA/CUDA-5.5/lib
ICPP_MAC = /usr/local/include
LCPP_MAC = /usr/local/lib
ILAPACK = /usr/include/lapacke

SDIR     = .
IDIR     = .
LDIR     = .

# Compiler for CUDA
NVCC      = nvcc

# CUDA compiling options
NVCCFLAGS =  -arch sm_35

# Compiler for C code
CXX       = g++

# Compiler for Markdown, requires Discount from yum
MD		  = markdown

# Standard optimization flags to C++ compiler
CXXFLAGS  = -O3 -I$(ICUDA) -I$(ICUDA_MAC) -I$(ICPP_MAC) -I$(ILAPACK)

# Add CUDA libraries to C++ compiler linking process
LDFLAGS  += -lcublas -lcurand -lcudart -L$(LCUDA) -L$(LCUDA_MAC) -L$(LCPP_MAC) -lnlopt -larmadillo -llapacke -lopenblas

# List Executables and Objects
EXEC = adrian 
OBJECTS = cppcode.o
OBJECTS1 = cppcodemulti.o

all : $(EXEC)

# Link objects from CUDA and C++ codes
$(EXEC) : adrian.o $(OBJECTS) 
	$(CXX) -o $@ $? $(LDFLAGS)

# Compile CUDA code
adrian.o : adrian.cu 
	$(NVCC) $(NVCCFLAGS) $(CXXFLAGS) -c $<  

# Compile C++ code
cppcode.o : cppcode.cpp
	$(CXX) $(CXXFLAGS) -c -std=c++11  $<

clean :
	rm -f *.o
	rm -f core core.*

veryclean :
	rm -f *.o
	rm -f core core.*
	rm -f $(EXEC)

run : all
	./$(EXEC)

doc : README.md
	$(MD) -o README.html README.md
