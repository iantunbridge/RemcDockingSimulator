INCLUDE=-I/$(HOME)/NVIDIA_GPU_Computing_SDK/C/common/inc -I/usr/include/nvidia-current/cuda -I/usr/local/cuda/include 
LIBS=-L/$(HOME)/NVIDIA_GPU_Computing_SDK/C/lib -L/usr/lib64/nvidia-current -L/usr/local/cuda/lib64  -L/usr/lib64
LINKS=-lpthread -lgsl -lgslcblas -lcudart -lcutil -lcuda
#TAU_MAKEFILE=$(HOME)/tau2/x86_64/lib/Makefile.tau-pthread-pdt
#COMPILER=tau_cxx.sh -tau_makefile=${TAU_MAKEFILE} -optPDTInst
#COMPILER=g++-4.3
COMPILER=g++
#config for zero.ch.cam.ac.uk
#INCLUDE=-I/usr/local/shared/rhel-4.5/x86_64/cuda/2.3/include -I/usr/local/shared/rhel-4.5/x86_64/cuda/2.3/C/common/inc -I/$(HOME)/include
#LIBS=-L/usr/local/shared/rhel-4.5/x86_64/cuda/2.3/cuda/lib64 -L/usr/local/shared/rhel-4.5/x86_64/cuda/2.3/C/lib -L/$(HOME)/lib

#set CUDA_PROFILE=1

#DEFINE=-DEnableOPENGL -DEnableCUDA -DEnableStreams 
### + debug symbols
#CFLAGS=-g 
CFLAGS=-fno-omit-frame-pointer -O2 #-ffast-math #

### the following line enables debug output and emulation
#NVCC_COMPILER_FLAGS=-g -deviceemu -D_EMU
NVCC_COMPILER_FLAGS=-g 
#NVCC_ARCH=-arch=sm_13

#if using ubuntu: link to -lgslcblas for gsl rng types
help:
	@echo Make using the following arguments
	@echo "  'make DEFINE="-DEnableCUDA" GPU'"
	@echo "  'make DEFINE="-DEnableCUDA -DEnableStreams" GPU' for async functionality"
	@echo "  'make CPU' for CPU"
	@echo "  'make DEFINE="-DEnableCUDA -DEnableOPENGL" GPUGL' for OPENGL, dont use it for simulations"
	
# compile this with  "make DEFINE=-DEnableCUDA GPU"
# compile this with  "make DEFINE=-DEnableCUDA -DEnableStreams GPU" for async functionality
GPU: obj clean.o AminoAcid.o AminoAcids.o Link.o main.o CudaFunctions.co Molecule.o Quaternion.o Replica.o Residue.o TorsionalLookupMatrix.o vector3f.o
	@echo Linking GPU Executable: REMCDockingGPU
	${COMPILER} ${INCLUDE} ${DEFINE} ${CFLAGS} ${LIBS} -o REMCDockingGPU obj/*.o ${LINKS}

# compile this with  DEFINE=-DEnableCUDA -DEnableOPENGL GPU
GPUGL: clean.o AminoAcid.o AminoAcids.o Camera.o Link.o main.o CudaFunctions.co Molecule.o Quaternion.o Replica.o Residue.o TorsionalLookupMatrix.o vector3f.o
	@echo Linking GPU Executable: REMCDockingGPU_OpenGL
	${COMPILER} ${INCLUDE} ${DEFINE} ${CFLAGS} ${LIBS} -o REMCDockingGPU_OpenGL obj/*.o ${LINKS} -lglut -lGL -lGLU

	
CPU: clean.o AminoAcid.o AminoAcids.o Camera.o Link.o main.o Molecule.o Quaternion.o Replica.o Residue.o TorsionalLookupMatrix.o vector3f.o
	@echo Linking CPU Executable: REMCDockingCPU
	${COMPILER} ${INCLUDE} ${DEFINE} ${CFLAGS} ${LIBS} -o REMCDockingCPU obj/*.o ${LINKS}

#DRMS: AminoAcid.o AminoAcids.o Camera.o Link.o Molecule.o Quaternion.o Replica.co Residue.o rngfunctions.o TorsionalLookupMatrix.o vector3f.o drms.o
#	@echo Making DRMS Executable: DRMS
#	g++ ${INCLUDE} ${DEFINE} ${CFLAGS} -o DRMS obj/*.o -lpthread -lgsl -lglut
#	rm -f drms.o
 
AminoAcid.o: src/AminoAcid.cpp src/AminoAcid.h 
	${COMPILER} ${CFLAGS} -c src/AminoAcid.cpp -o obj/AminoAcid.o

AminoAcids.o: src/AminoAcids.cpp src/AminoAcids.h src/AminoAcid.h 
	${COMPILER} ${CFLAGS} -c src/AminoAcids.cpp -o obj/AminoAcids.o

Camera.o: src/Camera.cpp src/Camera.h src/vector3f.h
	${COMPILER} ${CFLAGS} -c src/Camera.cpp -o obj/Camera.o

Link.o: src/Link.cpp src/Link.h
	${COMPILER} ${CFLAGS} -c src/Link.cpp -o obj/Link.o

main.o: src/main.cpp src/definitions.h  
	${COMPILER} ${CFLAGS} ${INCLUDE} ${DEFINE} -c src/main.cpp -o obj/main.o
	
Molecule.o: src/Molecule.cpp src/Molecule.h 
	${COMPILER} ${CFLAGS} ${INCLUDE} ${DEFINE} -c src/Molecule.cpp -o obj/Molecule.o

Quaternion.o: src/Quaternion.cpp src/Quaternion.h 
	${COMPILER} ${CFLAGS} -c src/Quaternion.cpp -o obj/Quaternion.o

Replica.o: src/Replica.cpp src/Replica.h src/vector3f.h
	${COMPILER} ${CFLAGS} ${INCLUDE} ${DEFINE} -c src/Replica.cpp -o obj/Replica.o

CudaFunctions.co: src/CudaFunctions.cu src/CudaFunctions.h
	@echo Making CUDA files.
	nvcc ${NVCC_ARCH} ${NVCC_COMPILER_FLAGS} ${INCLUDE} ${DEFINE} -c src/CudaFunctions.cu -o obj/CudaFunctions.o
	nvcc ${NVCC_COMPILER_FLAGS} ${INCLUDE} ${DEFINE} -cubin src/CudaFunctions.cu

Residue.o: src/Residue.cpp src/Residue.h src/vector3f.h
	${COMPILER} ${CFLAGS} -c src/Residue.cpp -o obj/Residue.o

TorsionalLookupMatrix.o: src/TorsionalLookupMatrix.cpp 
	${COMPILER} ${CFLAGS} -c src/TorsionalLookupMatrix.cpp -o obj/TorsionalLookupMatrix.o

vector3f.o: src/vector3f.cpp src/vector3f.h 
	${COMPILER} ${CFLAGS} -c src/vector3f.cpp -o obj/vector3f.o

configure:
	mkdir obj

clean:
	@echo Cleaning.
	@rm -f obj/*
	@rm -f REMCDockingGPU
	@rm -f REMCDockingCPU
	@rm -f *.cubin	
	@rm -f *.linkinfo
	@echo Done.
	
clean.o:
	@echo Cleaning object files for safety sake.
	@rm -f obj/*
