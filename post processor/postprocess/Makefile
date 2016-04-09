INCLUDE=-I/$(HOME)/NVIDIA_GPU_Computing_SDK/C/common/inc -I/usr/include/nvidia-current/cuda -I/usr/local/cuda/include 
LIBS=-L/$(HOME)/NVIDIA_GPU_Computing_SDK/C/lib -L/usr/lib64/nvidia-current -L/usr/lib/nvidia-current -L/usr/local/cuda/lib64 -L/usr/lib
LINKS=-lpthread -lgsl -lgslcblas -lcudart -lcutil -lcuda
COMPILER=g++-4.3

#set CUDA_PROFILE=1

#DEFINE=-DEnableOPENGL -DEnableCUDA -DEnableStreams 
### + debug symbols
CFLAGS=-g 
#CFLAGS=-fno-omit-frame-pointer -O3 #-ffast-math #

### the following line enables debug output and emulation
#NVCC_COMPILER_FLAGS=-g -deviceemu -D_EMU
NVCC_COMPILER_FLAGS=-g 
#NVCC_ARCH=-arch=sm_13

default: *.cpp
	${COMPILER} ${INCLUDE} ${CFLAGS} -c *.cpp
	${COMPILER} ${INCLUDE} ${DEFINE} ${CFLAGS} ${LIBS} -o postprocess *.o ${LINKS}
	mv postprocess ../processTraj

CudaFunctions.co: CudaFunctions.cu CudaFunctions.h
	@echo Making CUDA files.
	nvcc ${NVCC_ARCH} ${NVCC_COMPILER_FLAGS} ${INCLUDE} ${DEFINE} -c CudaFunctions.cu -o CudaFunctions.o
	#nvcc ${NVCC_COMPILER_FLAGS} ${INCLUDE} ${DEFINE} -cubin src/CudaFunctions.cu

configure:
	mkdir obj

clean:
	@echo Cleaning.
	@rm -f *.o
	@rm -f *.cubin	
	@rm -f *.linkinfo
	@echo Done.
	
clean.o:
	@echo Cleaning object files for safety sake.
	@rm -f *.o
