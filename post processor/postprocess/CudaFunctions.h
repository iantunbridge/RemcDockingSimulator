
#ifndef CUDAFUNCTIONS_H_
#define CUDAFUNCTIONS_H_
// cudaFunctions.h

using namespace std;

#include "definitions.h"

#include <cuda.h>
#include <cutil.h>  // CUDA c util package
#include <cuda_runtime_api.h>
#include <cutil_inline.h>
#include <cutil_inline_runtime.h>

#include "Replica.h"
#include "Molecule.h"

//#include "cudaExterns.h"

float * LJPotentialDataToDevice (AminoAcids *a);
void copyLJPotentialDataToDevice (float * dev_LJPotentialData, AminoAcids *a);

void MCSearchOnDevice();

// rotations and translations on gpu memory
void CUDA_rotateMolecule (float4 *d_residuePositions, int *d_startPosition, int *d_moleculeLength, int moleculeLength, float4* d_rotationVector, float4* d_center, cudaStream_t stream);
void CUDA_translateMolecule (float4 *d_residuePositions, int *d_startPosition, int *d_moleculeLength, int moleculeLength, float4* d_translation, float4* d_center, cudaStream_t stream);


void CUDA_EonDevice(float4 *residuePositions, float4 *residueMeta, int * residueCount, int *moleculePositions, int *moleculeCount, float* LJPotentials, double* result, int blockSize, int datasetSize, int sharedMemSize);
void CUDA_EonDeviceTest(float *d_x, float *d_y, float *d_z, int *d_id, float4 *residueMeta, float* LJPotentials, double* result, int blockSize, int datasetSize);


void cudaInfo();
//copy data to the device
int CUDA_memcpy_to_device(void * destination, void * source, int mem_size);
int CUDA_memcpy_to_host(void * destination, void * source, int mem_size);

// Asynchronous calls for use with streams
int CUDA_memcpy_to_device_async(void * destination, void * source, int mem_size, cudaStream_t stream);
int CUDA_memcpy_to_host_async(void * destination, void * source, int mem_size, cudaStream_t stream);
void CUDA_EonDevice_async(float4 *residuePositions, float4 *residueMeta, int *residueCount, int *moleculePositions, int *moleculeCount, float* LJPotentials, float* device_result, int resultSize, int blockSize, int datasetSize, int sharedMemSize, cudaStream_t stream);
void CUDA_Esum_async(float* result, float *d_resultMatrix, int resultSize, int datasetSize, cudaStream_t stream);

int CUDA_cudaMalloc(void ** destination, int mem_size);
int CUDA_cudaMemset(void * destination, int value, int mem_size);

void printCudaError(int code);
int CUDA_cudaFree(void ** destination);
int CUDA_createCudaStream(cudaStream_t *stream);

#if LJ_LOOKUP_METHOD == TEXTURE_MEM

void bindLJTexture(float * ljp);
void unbindLJTexture();
void bindLJTexture2D(float *ljp);
void unbindLJTexture2D ();

#endif

#if LJ_LOOKUP_METHOD == TEXTURE_MEM
	texture<float, 1, cudaReadModeElementType> LJTexture;
	//texture<float, 2, cudaReadModeElementType> LJTexture2D;
#elif LJ_LOOKUP_METHOD == CONST_MEM
	__constant__ float *LJPotentialData;  // store the LJ postentials in constant memory
#elif LJ_LOOKUP_METHOD == GLOBAL_MEM || LJ_LOOKUP_METHOD == SHARED_MEM
	float *LJPotentialData;  			// store the LJ postentials in global memory
#endif

// configure the threads to load into shared memory, ensures, minimum work done.
#if LJ_LOOKUP_METHOD == SHARED_MEM
	#if TILE_DIM == 32
		#define SM_LOAD_SIZE 50	// sizeof(float)*20*20 / 32 = 50
	#elif TILE_DIM == 64
		#define SM_LOAD_SIZE 25	// sizeof(float)*20*20 / 64 = 25
	#elif TILE_DIM == 128
		#define SM_LOAD_SIZE 13	// sizeof(float)*20*20 / 128 = 12.5, round up and pad, 13*128= 1664 shared mem used
	#elif TILE_DIM == 256
		#define SM_LOAD_SIZE 7	// sizeof(float)*20*20 / 256 = 6.25, round up and pad, 7*256= 1792 shared mem used
	#elif TILE_DIM == 512
		#define SM_LOAD_SIZE 4	// sizeof(float)*20*20 / 64 = 3.125, -> 4*512 => 2048
	#endif
#endif

#if METADATA_MEMORY == TEXTURE_MEM
	texture<float4, 1, cudaReadModeElementType> residueMetaTex;
#endif
#if POSITIONDATA_MEMORY == TEXTURE_MEM
	texture<float4, 1, cudaReadModeElementType> residuePositionTex;
#endif




__global__ void parallelSum_kernel(float * values);
__global__ void E_SimpleKernel(float4 * residuePositions, float4 * residueMeta, int * residueCount, int * moleculePositions, int * moleculeCount, float* LJPotentialData, float* result);
__global__ void E_TiledKernel(float4 * residuePositions, float4 * residueMeta, int * residueCount, int * moleculePositions, int * moleculeCount, float* LJPotentialData, float* result);
__global__ void E_TestTiledKernel(float *x, float *y, float *z, int *id , float4 * residueMeta, float* LJPotentialData, float* result);
__global__ void rotateMolecule_kernel (float4 *residuePositions, int *startPosition, int *length, float4* rotationVector, float4* center);
__global__ void translateMolecule_kernel (float4 *residuePositions, int *startPosition, int *length, float4* translationVector, float4* center);


#endif

