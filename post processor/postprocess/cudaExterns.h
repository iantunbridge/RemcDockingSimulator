/*
 * cudaExterns.h
 *
 * Prototype all cuda global functions such that they can be used in code compiled with only gcc
 *
 */

#ifndef CUDAEXTERNS_H_
#define CUDAEXTERNS_H_

extern "C"
{
	float * LJPotentialDataToDevice (AminoAcids *a);
	void copyLJPotentialDataToDevice (float * dev_LJPotentialData, AminoAcids *a);

	void MCSearchOnDevice();

	// rotations and translations on gpu memory
	void CUDA_rotateMolecule (float4 *d_residuePositions, int *d_startPosition, int *d_moleculeLength, int moleculeLength, float4* d_rotationVector, float4* d_center, cudaStream_t stream);
	void CUDA_translateMolecule (float4 *d_residuePositions, int *d_startPosition, int *d_moleculeLength, int moleculeLength, float4* d_translation, float4* d_center, cudaStream_t stream);


	void CUDA_EonDevice(float4 *residuePositions, float4 *residueMeta, int * residueCount, int *moleculePositions, int *moleculeCount, float* LJPotentials, double* result, int blockSize, int datasetSize, int sharedMemSize);
	void CUDA_EonDeviceTest(float *d_x, float *d_y, float *d_z, int *d_id, float4 *residueMeta, float* LJPotentials, double* result, int blockSize, int datasetSize);

	// initialise device
	void cudaInfo();
	//copy data to the device
	int CUDA_memcpy_to_device(void * destination, void * source, int mem_size);
	int CUDA_memcpy_to_host(void * destination, void * source, int mem_size);

	// Asynchronous calls for use with streams
	int CUDA_memcpy_to_device_async(void * destination, void * source, int mem_size, cudaStream_t stream);
	int CUDA_memcpy_to_host_async(void * destination, void * source, int mem_size, cudaStream_t stream);
	void CUDA_EonDevice_async(float4 *residuePositions, float4 *residueMeta, int *residueCount, int *moleculePositions, int *moleculeCount, float* LJPotentials, float* device_result, int resultSize, int blockSize, int datasetSize, int sharedMemSize, cudaStream_t stream);
	void CUDA_Esum_async(float* result, float *d_resultMatrix, int resultSize, int datasetSize, cudaStream_t stream);

	void printCudaError(int code);


	#if LJ_LOOKUP_METHOD == TEXTURE_MEM
	void bindLJTexture(float * ljp);
	void unbindLJTexture();
	void bindLJTexture2D(float *ljp);
	void unbindLJTexture2D ();
	#endif

#if METADATA_MEMORY == TEXTURE_MEM
	int bindMetaDataToTexture(void* deviceMemory, size_t size);
	int freeMetaDataTexture();
#endif

#if POSITIONDATA_MEMORY == TEXTURE_MEM
	int bindPositionDataToTexture(void* deviceMemory, size_t size);
	int freePositionDataTexture();
#endif

	cudaError_t cudaStreamCreate(cudaStream_t*);
	cudaError_t cudaStreamSynchronize(cudaStream_t);
	cudaError_t cudaStreamDestroy(cudaStream_t);
	cudaError_t cudaStreamQuery(cudaStream_t);
	cudaError_t cudaThreadSynchronize();
	const char* cudaGetErrorString(cudaError_t);
	cudaError_t cudaGetLastError();
	cudaError_t cudaGetDeviceCount(int*);
	cudaError_t cudaMalloc( void** devPtr, size_t size );
	cudaError_t cudaMemset( void* devPtr, int value ,size_t size );
	cudaError_t cudaMemcpy(void *, const void *, size_t, cudaMemcpyKind);
	cudaError_t cudaMemcpyAsync(void *, const void *, size_t, cudaMemcpyKind,cudaStream_t);
	cudaError_t cudaMallocHost( void** hostPtr, size_t size );
	cudaError_t cudaFreeHost( void* hostPtr );
	cudaError_t cudaFree( void* devPtr );
	cudaError_t cudaGetDevice( int* dev );
	cudaError_t cudaSetDevice( int dev );

};

#endif /* CUDAEXTERNS_H_ */
