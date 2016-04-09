// Replica Kernel functions to call cuda kernels for energy calculations in a replica
#ifndef CUDA_FUNCTIONS_CU
#define CUDA_FUNCTIONS_CU

using namespace std;

#include "CudaFunctions.h"
#include <cutil.h>
#include <cutil_inline.h>
// Block index
#define bx blockIdx.x
#define by blockIdx.y

// Thread index
#define tx threadIdx.x
#define ty threadIdx.y

#if USING_CUDA

#if LJ_LOOKUP_METHOD == TEXTURE_MEM

void bindLJTexture(float * ljp)
{
	size_t offset;
	cudaBindTexture(&offset, LJTexture, ljp, sizeof(float)*AA_COUNT*AA_COUNT);
	if (int(offset)>0)
		printf("!!! Texture offset %d\n",int(offset));
}

void unbindLJTexture()
{
	cudaUnbindTexture(LJTexture);
}

/*void bindLJTexture2D(float *ljp)
{
	const textureReference *texRefPtr;
	cudaGetTextureReference(&texRefPtr, "LJTexture2D");
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();

	LJTexture2D.addressMode[0] = cudaAddressModeWrap;
	LJTexture2D.addressMode[1] = cudaAddressModeWrap;
	LJTexture2D.filterMode = cudaFilterModePoint;
	LJTexture2D.normalized = false;

	size_t offset;
	cudaBindTexture2D(&offset, texRefPtr, ljp, &channelDesc, AA_COUNT, AA_COUNT, AA_COUNT*sizeof(float)+1);

	if (int(offset)>0)
		printf("!!! Texture offset %d\n",int(offset));
}*/


/*void unbindLJTexture2D ()
{
	const textureReference *texRefPtr;
	cudaGetTextureReference(&texRefPtr, "LJTexture2D");
	cudaUnbindTexture(texRefPtr);
}*/

#endif

#if METADATA_MEMORY == TEXTURE_MEM
int bindMetaDataToTexture(void* deviceMemory, size_t size)
{
	return cudaBindTexture(0, residueMetaTex, deviceMemory, size);
}

int freeMetaDataTexture()
{
	return cudaUnbindTexture(residueMetaTex);
}
#endif

#if POSITIONDATA_MEMORY == TEXTURE_MEM
int bindPositionDataToTexture(void* deviceMemory, size_t size)
{
	return cudaBindTexture(0, residuePositionTex, deviceMemory, size);
}

int freePositionDataTexture()
{
	return cudaUnbindTexture(residuePositionTex);
}
#endif

float * LJPotentialDataToDevice (AminoAcids *a)
{
	float * _LJPotentialData;
	long tableSize(sizeof(float[AA_COUNT*AA_COUNT]));
	long containerSize(tableSize);  // pad the container such that shared memory tables align nicely when preloading

	#if LJ_LOOKUP_METHOD == SHARED_MEM
		containerSize = TILE_DIM * ceil(float(tableSize)/float(TILE_DIM));
	#endif

	float * safeData = new float[AA_COUNT*AA_COUNT];

	for (int i=0;i<AA_COUNT;i++)
	{
		memcpy (&safeData[i*AA_COUNT],&a->LJpotentials[i][0],sizeof(float)*AA_COUNT);
	}

	CUDA_SAFE_CALL(cudaMalloc((void**)&_LJPotentialData,containerSize));
	CUDA_SAFE_CALL(cudaMemset(_LJPotentialData,0,containerSize));
	// copy host memory to device
	CUDA_SAFE_CALL(cudaMemcpy(_LJPotentialData, safeData, tableSize, cudaMemcpyHostToDevice) );
	delete [] safeData;
   	return _LJPotentialData;
}

void copyLJPotentialDataToDevice (float * dev_LJPotentialData, AminoAcids *a)
{
	#define mem_size int(sizeof(float)*AA_COUNT*AA_COUNT)
	#define table_size (AA_COUNT*AA_COUNT)
	float * safeData = new float[table_size];

	for (int i=0;i<table_size;i++)
		safeData[i] = 0.0f;

	for (int i=0;i<AA_COUNT;i++)
	{
		memcpy (&safeData[i*AA_COUNT],a->LJpotentials[i],sizeof(float)*AA_COUNT);
	}

	//CUDA_SAFE_CALL(cudaMalloc((void**)&dev_LJPotentialData,mem_size));
	// copy host memory to device
	cudaMemcpy(dev_LJPotentialData, safeData, mem_size, cudaMemcpyHostToDevice);
	cutilCheckMsg("Failed to copy contact potentials");
	delete [] safeData;
}

void cudaInfo()
{
	int deviceCount;
	CUDA_SAFE_CALL_NO_SYNC(cudaGetDeviceCount(&deviceCount));
	if (deviceCount == 0)
	{
		printf("ERROR: no devices supporting CUDA.\n");
	    exit(EXIT_FAILURE);
	}
	#if OUTPUT_LEVEL > 0
	struct cudaDeviceProp devInfo;
	printf("\n------------------------------------\n");
	for (int i=0;i<deviceCount;i++)
	{
		cudaGetDeviceProperties( &devInfo, 0 );
		printf ("CUDA Device %d Info:\n",i);
		printf("Name: %s\n", devInfo.name );
		printf("totalGlobalMem: %i\n",(int)devInfo.totalGlobalMem );
		printf("sharedMemPerBlock: %i\n",(int)devInfo.sharedMemPerBlock );
		printf("regsPerBlock: %i\n",devInfo.regsPerBlock );
		printf("warpSize: %i\n",devInfo.warpSize );
		printf("memPitch %i\n",(int)devInfo.memPitch );
		printf("maxThreadsPerBlock: %i\n",devInfo.maxThreadsPerBlock);
		if (i < deviceCount-1)
			printf("\n");
	}
	printf("------------------------------------\n");
	#endif
}

void MCSearchOnDevice()
{
	printf("MCSearchOnDevice() does nothing\n");
};

// stream kernel for overlapping calls
void CUDA_EonDevice_async(float4 *residuePositions, float4 *residueMeta, int * residueCount, int *moleculePositions, int *moleculeCount, float* LJPotentials, float* kernelResult, int resultSize, int blockSize, int datasetSize, int sm_size, cudaStream_t stream)
{
	dim3 dimBlock(blockSize,1,1);   // threads
	dim3 dimGrid(datasetSize/blockSize,datasetSize/blockSize,1);    	// blocks

	// Launch the device computation
	// result must become an array of dimentions the same as grid so that the multiple thread blocks can sum in it properly
	//cudaStreamSynchronize(stream);

	E_TiledKernel<<< dimGrid,dimBlock,sm_size,stream >>>(residuePositions, residueMeta, residueCount, moleculePositions, moleculeCount, LJPotentials, kernelResult);

	cutilCheckMsg("Kernel execution failed");
	//cudaMemcpyAsync(hostResult,kernelResult, sizeof(float)*resultSize*resultSize, cudaMemcpyDeviceToHost,stream);
	return;
};

void CUDA_Esum_async(float* result, float *d_resultMatrix, int resultSize, int datasetSize, cudaStream_t stream)
{
	// parallel sum on gpu if required
}
/*
void CUDA_EonDevice_B(float4 *residuePositions, float4 *residueMeta, int * residueCount, int *moleculePositions, int *moleculeCount, float* LJPotentials, float* result, int blockSize, int datasetSize)
{
	// result stored on the device
	// gridSize can be arbitrary
	int gridSize = datasetSize/blockSize;
	// the parallel sum of the results must have an array with size that is a power of 2
	int resultSize = int(pow(2,ceil(log(float(gridSize))/log(2.0))));
	//resultSize *= resultSize;

	//for a parallel sum each grid must have one cell in the array of results from all the threads
	float *d_result;
	cudaMalloc((void **)&d_result,sizeof(float)*resultSize*resultSize);
	cudaMemset(d_result,0,sizeof(float)*resultSize*resultSize);


	// Launch the device computation
	// result must become an array of dimentions the same as grid so that the multiple thread blocks can sum in it properly

	dim3 dimBlock(blockSize,blockSize,1);   // threads
	dim3 dimGrid(gridSize,gridSize,1);    	// blocks

	E_SimpleKernel<<< dimGrid,dimBlock >>>(residuePositions, residueMeta, residueCount, moleculePositions, moleculeCount, LJPotentials, d_result);


	//sum the results from the thread blocks on the GPU
	 fix this for multiple reduction sums!!!

	//    float * sumResults;
	//
	//    if (resultSize*resultSize<=512) //options are 1 ... 256 == 16x16 thread blocks, ie fits in one kernel
	//    {
	//    	parallelSum_kernel<<<1,resultSize*resultSize>>>(d_result,d_stride);  // THIS WORKS!
	//	}
	//    else   // resultsize > 512 == 32x32,64x64,128x128,256x256 .... etc
	//    {
	//    	int threadBlockSize = 512;				// maximum coverage by threads
	//    	gridSize = resultSize*resultSize/512;   // how many grids we need
	//	   	cudaMalloc((void **)&sumResults,sizeof(float)*gridSize);
	//		cudaMemset(sumResults,0,sizeof(float)*gridSize);
	//	   	parallelSum_kernel<<<gridSize,threadBlockSize>>>(d_result,sumResults);
	//
	//	}
	//	cudaFree(d_stride);
	//	cudaMemcpy(result,d_result, sizeof(float), cudaMemcpyDeviceToHost);

	#ifdef _EMU

	printf("After || sum\n");
	for (int i=0;i<resultSize;i++)
	{
		for (int j=0;j<resultSize;j++)
			printf (" %10.7f",d_result[i*resultSize+j]);
		printf("\n");
	}
	#endif
	// copy the result back to host memory
	// sum it on the CPU, seems faster than a kernel for 64x64 at least

		float *tmpSums = new float [resultSize*resultSize];
		cudaMemcpy(tmpSums,d_result, sizeof(float)*resultSize*resultSize, cudaMemcpyDeviceToHost);
		result[0] = 0;
		for (int i=0;i<resultSize*resultSize;i++)
			result[0] += tmpSums[i];
		result[0] *= KBTConversionFactor;

	// free the memory assigned for this iteration.
	cudaFree(d_result);

	return;
};


void CUDA_rotateMolecule (float4 *d_residuePositions, int *d_startPosition, int *d_moleculeLength, int moleculeLength, float4* d_rotationVector, float4* d_center, cudaStream_t stream)
{
	int grids = int(ceil(float(moleculeLength)/512.0f));
	int threads = int(ceil(float(moleculeLength)/float(grids)));
	dim3 dimBlock(moleculeLength,1,1);   // threads
	dim3 dimGrid(grids,1,1);    	// blocks

	rotateMolecule_kernel<<<dimGrid,dimBlock,0,stream>>>(d_residuePositions, d_startPosition, d_moleculeLength, d_rotationVector, d_center);
}
*/
/*
void CUDA_translateMolecule (float4 *d_residuePositions, int *d_startPosition, int *d_moleculeLength, int moleculeLength, float4* d_translation, float4* d_center, cudaStream_t stream)
{
	int grids = int(ceil(float(moleculeLength)/512.0f));
	int threads = int(ceil(float(moleculeLength)/float(grids)));
	dim3 dimBlock(moleculeLength,1,1);   // threads
	dim3 dimGrid(grids,1,1);    	// blocks

	translateMolecule_kernel <<< dimGrid,dimBlock,0,stream >>> (d_residuePositions, d_startPosition, d_moleculeLength, d_translation, d_center);
}

*/
// alternative summation algorithm will work with up to 33554432 residues
void CUDA_EonDevice(float4 *residuePositions, float4 *residueMeta, int * residueCount, int *moleculePositions, int *moleculeCount, float* LJPotentials, double* result, int blockSize, int datasetSize, int sm_size)
{

	// result stored on the device
	// gridSize can be arbitrary
	int gridSize = datasetSize/blockSize;
	// the parallel sum of the results must have an array with size that is a power of 2
	int resultSize = gridSize; //int(pow(2,ceil(log(float(gridSize))/log(2.0))));
	float *tmpSums = new float[resultSize*resultSize];

	//for a parallel sum each grid must have one cell in the array of results from all the threads
	float *d_result;
	cudaMalloc((void **)&d_result,sizeof(float)*resultSize*resultSize);
	cudaMemset(d_result,0,sizeof(float)*resultSize*resultSize);

	//float4 testmem;
	//cudaMemcpy(&testmem,residuePositions, sizeof(float4), cudaMemcpyDeviceToHost);


	// Launch the device computation
	// result must become an array of dimensions the same as grid so that the multiple thread blocks can sum in it properly

	dim3 dimBlock(blockSize,1,1);   		// threads
	dim3 dimGrid(gridSize,gridSize,1);    	// blocks


	E_TiledKernel<<< dimGrid,dimBlock,sm_size,0 >>>(residuePositions, residueMeta, residueCount, moleculePositions, moleculeCount, LJPotentials, d_result);
    	// check if kernel invocation generated an error

#ifdef _EMU

	printf("After kernel:\n");
	for (int i=0;i<resultSize;i++)
	{
		for (int j=0;j<resultSize;j++)
			printf (" %10.6f",d_result[i*resultSize+j]);
		printf("\n");
	}
#endif

#if PARALLEL_SUM
	parallelSum_kernel<<<1,gridSize>>>(d_result);
	cudaMemcpy(result,d_result, sizeof(float), cudaMemcpyDeviceToHost);
	result[0] *= KBTConversionFactor;
#else
	cudaMemcpy(tmpSums,d_result, sizeof(float)*resultSize*resultSize, cudaMemcpyDeviceToHost);
	result[0] = 0.0f;
	float c(0.0f);
	/*
			function kahanSum(input)
			 var potential = 0
			 var c = 0
			 for i = 0 to blockdim-1
			  y = lj+dh - c
			  t = potential + y
			  c = (t - potential) - y
			  potential = t
			 next i
			return sum
			*/

	for (int i=0;i<resultSize*resultSize;i++)
	{
		#if COMPENSATE_KERNEL_SUM
		float y(tmpSums[i] -c);
		float t((*result) + y);
		c = (t-(*result)) -y;
		*result = t;
		#else
		result[0] += tmpSums[i];
		#endif
	}
	result[0] *= KBTConversionFactor;

	delete [] tmpSums;
#endif
	// free the memory assigned for this iteration.
	cudaFree(d_result);
	cutilCheckMsg("Kernel execution failed");
	return;
};

void CUDA_EonDeviceTest(float *d_x, float *d_y,float *d_z, int *d_id, float4 *residueMeta, float* LJPotentials, float* result, int blockSize, int datasetSize)
{

	// result stored on the device
	// gridSize can be arbitrary
	int gridSize(datasetSize/blockSize);
	float *tmpSums = new float[gridSize*gridSize];

	//for a parallel sum each grid must have one cell in the array of results from all the threads
	float *d_result;
	cudaMalloc((void **)&d_result,sizeof(float)*gridSize*gridSize);
	cudaMemset(d_result,0,sizeof(float)*gridSize*gridSize);

	//float4 testmem;
	//cudaMemcpy(&testmem,residuePositions, sizeof(float4), cudaMemcpyDeviceToHost);


	// Launch the device computation
	// result must become an array of dimensions the same as grid so that the multiple thread blocks can sum in it properly

	dim3 dimBlock(blockSize,1,1);   		// threads
	dim3 dimGrid(gridSize,gridSize,1);    	// blocks

	printf ("+++++++++++++++++++ NO KERNEL EXECUTED +++++++++++++++++++++++++++++");
	//E_TestTiledKernel<<< dimGrid,dimBlock >>>(d_x,d_y,d_z,d_id, residueMeta, LJPotentials, d_result);

	cudaMemcpy(tmpSums,d_result, sizeof(float)*gridSize*gridSize, cudaMemcpyDeviceToHost);
	result[0] = 0.0f;
	for (int i=0;i<gridSize*gridSize;i++)
	{
		result[0] += tmpSums[i];
	}
	result[0] *= KBTConversionFactor;
	delete [] tmpSums;

	// free the memory assigned for this iteration.
	cudaFree(d_result);
    // check if kernel invocation generated an error
    cutilCheckMsg("Kernel execution failed");

	return;
};

/*void CUDA_EonDeviceLoop(float4 *residuePositions, float4 *residueMeta, int * residueCount, int *moleculePositions, int *moleculeCount, float* LJPotentials, float* result, int blockSize, int datasetSize)
{
	// result stored on the device
	// gridSize can be arbitrary
	int gridSize = datasetSize/blockSize;
	// the parallel sum of the results must have an array with size that is a power of 2
	int resultSize = gridSize; //int(pow(2,ceil(log(float(gridSize))/log(2.0))));

	//for a parallel sum each grid must have one cell in the array of results from all the threads
	float *d_result;
	cudaMalloc((void **)&d_result,sizeof(float)*resultSize);
	cudaMemset(d_result,0,sizeof(float)*resultSize);

	// Launch the device computation
	// result must become an array of dimensions the same as grid so that the multiple thread blocks can sum in it properly
	dim3 dimBlock(blockSize,1,1);   		// threads
	dim3 dimGrid(gridSize,1,1);    	// blocks



#ifdef _EMU

	printf("After kernel:\n");
	for (int i=0;i<resultSize;i++)
	{
		for (int j=0;j<resultSize;j++)
			printf (" %10.6f",d_result[i*resultSize+j]);
		printf("\n");
	}
#endif

#if PARALLEL_SUM
	parallelSum_kernel<<<1,gridSize>>>(d_result);
	cudaMemcpy(result,d_result, sizeof(float), cudaMemcpyDeviceToHost);
	result[0] *= KBTConversionFactor;
#else
	float *tmpSums = new float[resultSize];
	cudaMemcpy(tmpSums,d_result, sizeof(float)*resultSize, cudaMemcpyDeviceToHost);
	result[0] = 0.0f;
	for (int i=0;i<resultSize;i++)
	{
		result[0] += tmpSums[i];
	}
	result[0] *= KBTConversionFactor;
	delete [] tmpSums;
#endif
	// free the memory assigned for this iteration.
	cudaFree(d_result);
	return;
};


__global__ void MCSearchOnDevice_kernel()
{
	return;
};
*/

__device__ float length(const float4 a, const float4 b)
{
	return sqrtf((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z));

};

__device__ float crowderPotential(float r)
{
	return __powf(6.0f/r,6.0f);
};

__device__ float dhPotential(float q1, float q2, float r)
{
	return  q1 * q2 *__expf(-r/Xi) / r ;
};
/*
__global__ void E_SimpleKernel(float4 * residuePositions, float4 * residueMeta, int * residueCount, int * moleculePositions, int * moleculeCount, float* LJPotentialData, float* result)
{
	#define	 residueIp residuePositions[bx*blockDim.x+tx]
	#define residueIm residueMeta[bx*blockDim.x+tx]

	#define	 residueJp residuePositions[by*blockDim.x+ty]
	#define residueJm residueMeta[by*blockDim.x+ty]

	__shared__ float sharedmem_results[SHARED_MEM_RESULT_SIZE];
			// array of temp results to sum, can never be > 512
			// because of the number of threads allowed on a gpu
			// and is the only fixed memoy requirement I have in
			// this code so its set at 512


	// set each tmp_result to 0 if we are going to be using it.
	// problems here are that if the blockdim is not a power of 2 or too small, not all the
	// required cells are set to 0, resulting in an incorrect summation. the second line fixes
	// this provided SHARED_MEM_RESULT_SIZE is eq the next power of 2 up from the blocksize
	// squared. eg: BLOCK_SIZE == 12 then 2^ ceil((log(12)/log(2)) == 16 => SHARED_MEM_RESULT_SIZE == 256
	sharedmem_results[tx*blockDim.x+ty] = 0.0f;
	sharedmem_results[SHARED_MEM_RESULT_SIZE-(tx*blockDim.x-ty)] = 0.0f;

	#if LJ_LOOKUP_METHOD == SHARED_MEM_LJ
	// copy the LJ lookup table into shared memory
	__shared__ float sharedLJ[400];
	for (int i=0;i< ceil(400.0f/(blockDim.x*blockDim.y)) ;i++)
		sharedLJ[tx*blockDim.x+ty +blockDim.x*blockDim.y*i] = LJPotentialData[tx*blockDim.x+ty +blockDim.x*blockDim.y*i]; // basically there are nxn threads so use this to quickly fill the shared memory array
	#endif
	#if CULL_LOWER_TRIANGULAR_SUM  // only works if we throw away the entire block, otherwise results cannot be guaranteed correct due to the parallel sum
	if (bx < by)
	{
		return;
	}
	#endif
	__syncthreads();  // very important so that the entire shared mem array is assigned before any threads use it

	//p.w < 0 == padding element
	//m.x = aminoacid index
	//m.y = charge
	//m.z = vdw radius
	//m.w = crowder if == CROWDER_IDENTIFIER

	float4 rIp = residuePositions[bx*blockDim.x+tx];
	float4 rJp = residuePositions[by*blockDim.x+ty];

	if (rIp.w == rJp.w || rIp.w < CROWDER_IDENTIFIER || rJp.w < CROWDER_IDENTIFIER )		// same molecule || padding residues || padding residues
	{
		// residueIp.w == residueJp.w means that they are the same molecule
		// if either is -1 then its a padding residue and must be ignored
	}
	else
	{

		float4 rIm = residueMeta[bx*blockDim.x+tx];
		float4 rJm = residueMeta[by*blockDim.x+ty];

		float r = length(rIp,rJp) + EPS;  // add eps so that r is never 0, can happen in a collision


		//if there are crowders the LJ is replaced by the repulsive formula.
		float LJ;   // needs to be here for scoping issues
		float DH;

		#if REPULSIVE_CROWDING
		if (rIm.w == CROWDER_IDENTIFIER || rJm.w == CROWDER_IDENTIFIER) // repulsive crowder interaction
		{
			LJ = crowderPotential(r);
		}
		else  // normal LJ interaction
		{
		#endif
			// do the debye huckel long range intereactions
			DH = rIm.y * rJm.y * 1.602176487f * 1.602176487f * expf(-r/Xi) / r * DH_CONVERSION_FACTOR;
			//int ijX =  rint(rIm.x + float(AA_COUNT)*rJm.x);
			float sigmaij = (rIm.z + rJm.z ) * 0.5f;

		#if	LJ_LOOKUP_METHOD == SHARED_MEM
			float Eij = lambda*(sharedLJ[ijX] - e0);
		#elif LJ_LOOKUP_METHOD == TEXTURE_MEM
			float Eij = lambda*(tex2D(LJTexture,rIm.x,rJm.x) - e0);
		#elif LJ_LOOKUP_METHOD == CONST_MEM || LJ_LOOKUP_METHOD == GLOBAL_MEM
			float Eij = lambda*(LJPotentialData[ijX] - e0);
		#endif

			// sigmaij is the average atomic radius determined by the van der waals radius in kim2008
			//float r0 = sigmaij*1.122462048309372981433533049679f; //pow(2.0,(1.0/6.0));
			// the following 2 lines are needed for preserving numerical accuracy on a cpu. not needed on GPU
			float LJtmp = __powf(sigmaij/r, 6.0f);

			LJ = -4.0f*Eij*LJtmp*(LJtmp-1.0f);
			if (Eij>0.0f && r<sigmaij*1.122462048309372981433533049679f)
			{
				LJ = -LJ + 2.0f*Eij;
			}

		#if REPULSIVE_CROWDING
		}  // end conditional branch for LJ or repulsive short-range energy
		#endif
		//LJ = LJ * LJ_CONVERSION_FACTOR;  // convert from units used to K_bT

		// store the sum of all potentials for this thread so that we dont store to global mem until we have to.
#if CULL_LOWER_TRIANGULAR_SUM
		if (bx==by)
		{
			sharedmem_results[tx*blockDim.x+ty] = 0.5f*(LJ * LJ_CONVERSION_FACTOR + DH); // convert to kcal/mol
		}
		else // if (bx > by)
		{
			sharedmem_results[tx*blockDim.x+ty] = (LJ * LJ_CONVERSION_FACTOR + DH); // convert to kcal/mol
		}
#else
		sharedmem_results[tx*blockDim.x+ty] = 0.5f*(LJ * LJ_CONVERSION_FACTOR + DH) ; // convert to kcal/mol
#endif //  CULL_LOWER_TRIANGULAR_SUM
	}

	__syncthreads();  // all threads must have computed values so that they can be read concurrently from shared memory
	for(int bitMask=SHARED_MEM_RESULT_SIZE/2; bitMask>0; bitMask/=2)
	{
		float t=sharedmem_results[tx*blockDim.x+ty] + sharedmem_results[(tx*blockDim.x+ty)^bitMask];
		__syncthreads(); // make sure all values are read before modifying them
		sharedmem_results[tx*blockDim.x+ty]=t;
		__syncthreads(); // make sure all values are written back before the next iteration
	}

	// write results of the calculation back to global memory using a reduction summation

	result[bx*gridDim.x+by] = sharedmem_results[0];

	return;
};
*/
/*
__global__ void E_TiledKernelSMopt(float4 * residuePositions, float4 * residueMeta, int * residueCount, int * moleculePositions, int * moleculeCount, float* LJPotentialData, float* result)
{
	extern __shared__ float sharedmem_results[];	// array of temp results to sum, one per thread run.
	const int positionOfPosition = blockDim.x/4;
	extern __shared__ float4 X_tile_residuePositions[];
	const int positionOfMeta = __mul24(5,blockDim.x)/4;
	extern __shared__ float4 X_tile_residueMeta[];

	//__shared__ float sharedmem_results[TILE_DIM];	// array of temp results to sum, one per thread run.
	sharedmem_results[tx] = 0.0f;//EPS;
	#if  CULL_LOWER_TRIANGULAR_SUM  // only works if we throw away the entire block, otherwise results cannot be guaranteed correct due to the parallel sum
	if (bx < by)
		return;
	#endif

	#if POSITIONDATA_MEMORY == TEXTURE_MEM
		float4 yresiduep = tex1Dfetch(residuePositionTex, by*blockDim.x+tx);
	#else
		float4 yresiduep = residuePositions[by*blockDim.x+tx];
	#endif
	#if METADATA_MEMORY == TEXTURE_MEM
		float4 yresiduem = tex1Dfetch(residueMetaTex, by*blockDim.x+tx);
	#else
		float4 yresiduem = residueMeta[by*blockDim.x+tx];
	#endif
	// think of the position as a column of residues so 0 == by:0,tx:0, N = by:gridDimY*bDX+tx

	// copy residues from GRAM to shared memory tiles

	//__shared__ float4 X_tile_residuePositions[TILE_DIM];
	//__shared__ float4 X_tile_residueMeta[TILE_DIM];
	#if POSITIONDATA_MEMORY == SHARED_MEM
		X_tile_residuePositions[positionOfPosition + tx] = residuePositions[bx*blockDim.x+tx];
	#endif
	#if METADATA_MEMORY == SHARED_MEM
		X_tile_residueMeta[positionOfMeta + tx] = residueMeta[bx*blockDim.x+tx];
	#endif

	__syncthreads();  // very important so that the entire shared mem array is assigned before any threads use it

	if (yresiduep.w > PADDER_IDENTIFIER) // the residue in question is not a padding element
	{
		float c_lj(0.0f);   // Kahan summation variables
		float c_dh(0.0f);   // Kahan summation variables
		float lj_subtotal(0.0f);
		float dh_subtotal(0.0f);
		float y; // Kahan summation variables
		float t; // Kahan summation variables

		//#pragma unroll 4 - cant unroll dynamic variable
		for (int i = 0; i < blockDim.x; i++)
		{
			//#if POSITIONDATA_MEMORY == TEXTURE_MEM
			//	float4 POS = tex1Dfetch(residuePositionTex,(bx*blockDim.x+i));
			//#endif

			if (yresiduep.w == X_tile_residuePositions[positionOfPosition + i].w || X_tile_residuePositions[positionOfPosition + i].w < 0.0f )		// same molecule || padding residue
			{
				// Xpos.w == Ypos.w means that they are the same molecule
				// if either is -2 then its a padding residue and must be ignored
				// DO NOTHING
			}
			else
			{
				//#if METADATA_MEMORY == TEXTURE_MEM
				//float4 meta = tex1Dfetch(residueMetaTex,(bx*blockDim.x+i));
				//#endif

				//m.x = aminoacid index
				//m.y = charge
				//m.z = vdw radius
				//m.w = crowder if == CROWDER_IDENTIFIER
				float r(length(yresiduep,X_tile_residuePositions[positionOfPosition + i]) + EPS);  // add eps so that r is never 0, can happen in a collision

				float LJ(0.0f);
				float DH(0.0f);

				#if REPULSIVE_CROWDING
					//#if METADATA_MEMORY == TEXTURE_MEM
					//if (yresiduem.w == CROWDER_IDENTIFIER || meta.w == CROWDER_IDENTIFIER) // repulsive crowder interaction
					//#else
					if (yresiduem.w == CROWDER_IDENTIFIER || X_tile_residueMeta[positionOfMeta + i].w == CROWDER_IDENTIFIER) // repulsive crowder interaction
					//#endif
					{
						LJ = crowderPotential(r);//__powf(6.0f/r,6);
					}
					else  // normal LJ interaction
					{
				#endif

				//do the texture fetch first
				#if LJ_LOOKUP_METHOD == TEXTURE_MEM
					int ijX(rint( AA_COUNT*yresiduem.x + X_tile_residueMeta[positionOfMeta + i].x));
					float Eij(lambda*(tex1Dfetch(LJTexture,ijX) - e0));
					//float Eij(lambda*(e0 - e0));
					//float Eij(lambda*(tex2D(LJTexture2D,rint(X_tile_residueMeta[positionOfMeta + i].x),rint(yresiduem.x)) - e0));
				#else
					int ijX(rint( AA_COUNTf*X_tile_residueMeta[positionOfMeta + i].x + yresiduem.x));
					#if LJ_LOOKUP_METHOD == SHARED_MEM
					float Eij = lambda*(sharedLJ[ijX] - e0);
					#else  // __global__ or __constant__
					float Eij = lambda*(LJPotentialData[ijX] - e0);
					#endif
				#endif

				//#if METADATA_MEMORY == TEXTURE_MEM
				//	int ijX =  rint( yresiduem.x + AA_COUNTf*meta.x);

				//	DH = dhPotential(yresiduem.y,meta.y,r);
				//#else

					//old int ijX(rint( yresiduem.x + AA_COUNTf*X_tile_residueMeta[i].x));
					//old lookup for 1d tex
					//int ijX(rint( AA_COUNTf*X_tile_residueMeta[positionOfMeta + i].x + yresiduem.x));

					DH = dhPotential(yresiduem.y,X_tile_residueMeta[positionOfMeta + i].y,r);
				//#endif

				// compensate sum
				y = DH - c_dh;
				t = dh_subtotal + y;
				c_lj = (t-dh_subtotal)-y;
				dh_subtotal = t;




					// sigmaij is the average atomic radius determined by the van der waals radius in kim2008
//					float sigmaij = __fmul_rn(__fadd_rn(yresiduem.z,X_tile_residueMeta[i].z), 0.5f);
					//#if METADATA_MEMORY == TEXTURE_MEM
					//	float sigmaij((yresiduem.z + meta.z) * 0.5f);
					//#else
						float sigmaij((yresiduem.z + X_tile_residueMeta[positionOfMeta + i].z) * 0.5f);
					//#endif

					float LJtmp(__powf(sigmaij/r,6.0f));
					//float r0 = sigmaij*1.122462048309372981433533049679f; //sigma*pow(2.0,(1.0/6.0));

					LJ = -4.0f*Eij*LJtmp*(LJtmp-1.0f);
					//LJ = __fmul_rn(-4.0f,__fmul_rn(Eij,__fmul_rn(LJtmp,__fadd_rn(LJtmp,-1.0f))));
					if (Eij>0.0f && r < (sigmaij*1.122462048309372981433533049679f))
					{
						LJ = -LJ + 2.0f*Eij;
					}

				#if REPULSIVE_CROWDING
				}  // end conditional branch for LJ or repulsive short-range energy
				#endif
				//sharedmem_results[tx] += (LJ * LJ_CONVERSION_FACTOR) + (DH * DH_CONVERSION_FACTOR * 1.602176487f * 1.602176487f );
				y = LJ - c_lj;
				t = lj_subtotal + y;
				c_lj = (t-lj_subtotal)-y;
				lj_subtotal = t;


			} // if !(X_tile_residuePositions.w == Y_tile_residuePositions.w || X_tile_residuePositions.w < CROWDER_IDENTIFIER )
			sharedmem_results[tx] = (lj_subtotal * LJ_CONVERSION_FACTOR) + (dh_subtotal * DH_CONVERSION_FACTOR * 1.602176487f * 1.602176487f);
		} // for i = 0..Bdx

		#if CULL_LOWER_TRIANGULAR_SUM
		if (bx==by)
			sharedmem_results[tx] = sharedmem_results[tx]*0.5f;
		#else
		sharedmem_results[tx] *= 0.5f
		#endif //  CULL_LOWER_TRIANGULAR_SUM
	} // if !padding

	__syncthreads();  // all threads must have computed values so that they can be read concurrently from shared memory

	// do reduction, see SDK reduction example
	// FULLY unrolled reduction

	// do reduction in shared mem
    if (blockDim.x >= 512) { if (tx < 256) { sharedmem_results[tx] += sharedmem_results[tx + 256]; } __syncthreads(); }
    if (blockDim.x >= 256) { if (tx < 128) { sharedmem_results[tx] += sharedmem_results[tx + 128]; } __syncthreads(); }
    if (blockDim.x >= 128) { if (tx <  64) { sharedmem_results[tx] += sharedmem_results[tx + 64];  } __syncthreads(); }

	if (tx < 32)
	{
		if (blockDim.x >=  64) { sharedmem_results[tx] += sharedmem_results[tx + 32]; }
		if (blockDim.x >=  32) { sharedmem_results[tx] += sharedmem_results[tx + 16]; }
		if (blockDim.x >=  16) { sharedmem_results[tx] += sharedmem_results[tx +  8]; }
		if (blockDim.x >=   8) { sharedmem_results[tx] += sharedmem_results[tx +  4]; }
		if (blockDim.x >=   4) { sharedmem_results[tx] += sharedmem_results[tx +  2]; }
		if (blockDim.x >=   2) { sharedmem_results[tx] += sharedmem_results[tx +  1]; }
	}

	// write results of the calculation back to global memory, only thread 0 has the correct answer
	// only one warp processes < 32, therefore no syncthreads required here
	if (tx == 0)
		result[bx*gridDim.x+by] = sharedmem_results[tx];

	return;
};*/

///* full kernel with optional texture containers
__global__ void E_TiledKernel(float4 * residuePositions, float4 * residueMeta, int * residueCount, int * moleculePositions, int * moleculeCount, float* LJPotentialData, float* result)
{
	extern __shared__ float sharedmem_results[];	// array of temp results to sum, one per thread run.
	#if POSITIONDATA_MEMORY == SHARED_MEM
	const int positionOfPosition = blockDim.x/4;
	extern __shared__ float4 X_tile_residuePositions[];
	#endif
	#if METADATA_MEMORY == SHARED_MEM
	#if POSITIONDATA_MEMORY == SHARED_MEM
	const int positionOfMeta = __mul24(5,blockDim.x)/4;
	#else
	const int positionOfMeta = blockDim.x/4;
	#endif
	extern __shared__ float4 X_tile_residueMeta[];
	#endif

	//__shared__ float sharedmem_results[TILE_DIM];	// array of temp results to sum, one per thread run.
	sharedmem_results[tx] = 0.0f;//EPS;
	#if  CULL_LOWER_TRIANGULAR_SUM  // only works if we throw away the entire block, otherwise results cannot be guaranteed correct due to the parallel sum
	if (bx < by)
		return;
	#endif

	#if POSITIONDATA_MEMORY == TEXTURE_MEM
		float4 yresiduep = tex1Dfetch(residuePositionTex, by*blockDim.x+tx);
	#else
		float4 yresiduep = residuePositions[by*blockDim.x+tx];
	#endif
	#if METADATA_MEMORY == TEXTURE_MEM
		float4 yresiduem = tex1Dfetch(residueMetaTex, by*blockDim.x+tx);
	#else
		float4 yresiduem = residueMeta[by*blockDim.x+tx];
	#endif
	// think of the position as a column of residues so 0 == by:0,tx:0, N = by:gridDimY*bDX+tx

	// copy residues from GRAM to shared memory tiles

	//__shared__ float4 X_tile_residuePositions[TILE_DIM];
	//__shared__ float4 X_tile_residueMeta[TILE_DIM];
	#if POSITIONDATA_MEMORY == SHARED_MEM
		X_tile_residuePositions[positionOfPosition + tx] = residuePositions[bx*blockDim.x+tx];
	#endif
	#if METADATA_MEMORY == SHARED_MEM
		X_tile_residueMeta[positionOfMeta + tx] = residueMeta[bx*blockDim.x+tx];
	#endif

	__syncthreads();  // very important so that the entire shared mem array is assigned before any threads use it

	if (yresiduep.w > PADDER_IDENTIFIER) // the residue in question is not a padding element
	{
		float lj_subtotal(0.0f);
		float dh_subtotal(0.0f);
		#if COMPENSATE_KERNEL_SUM
		float c_lj(0.0f);   // Kahan summation variable
		float c_dh(0.0f);
		float y;
		float t;

		//		function kahanSum(input)
		//		 var potential = 0
		//		 var c = 0
		//		 for i = 0 to blockdim-1
		//		  y = lj+dh - c
		//		  t = potential + y
		//		  c = (t - potential) - y
		//		  potential = t
		//		 next i
		//		return sum

		#endif
		//#pragma unroll 4 - cant unroll dynamic variable
		for (int i = 0; i < blockDim.x; i++)
		{
			#if POSITIONDATA_MEMORY == TEXTURE_MEM
			float4 pos = tex1Dfetch(residuePositionTex,(bx*blockDim.x+i));
			if (yresiduep.w == pos.w || pos.w < 0.0f )		// same molecule || padding residue
			#else
			if (yresiduep.w == X_tile_residuePositions[positionOfPosition + i].w || X_tile_residuePositions[positionOfPosition + i].w < 0.0f )		// same molecule || padding residue
			#endif
			{
				// Xpos.w == Ypos.w means that they are the same molecule
				// if either is -2 then its a padding residue and must be ignored
				// DO NOTHING
			}
			else
			{
				#if METADATA_MEMORY == TEXTURE_MEM
				float4 meta = tex1Dfetch(residueMetaTex,(bx*blockDim.x+i));
				#endif
				//m.x = aminoacid index
				//m.y = charge
				//m.z = vdw radius
				//m.w = crowder if == CROWDER_IDENTIFIER
				#if POSITIONDATA_MEMORY == TEXTURE_MEM
				float r(length(yresiduep,pos) + EPS);  // add eps so that r is never 0, can happen in a collision
				#else
				float r(length(yresiduep,X_tile_residuePositions[positionOfPosition + i]) + EPS);  // add eps so that r is never 0, can happen in a collision
				#endif
				float LJ(0.0f);
				float DH(0.0f);

				#if REPULSIVE_CROWDING
					#if METADATA_MEMORY == TEXTURE_MEM
					if (yresiduem.w == CROWDER_IDENTIFIER || meta.w == CROWDER_IDENTIFIER) // repulsive crowder interaction
					#else
					if (yresiduem.w == CROWDER_IDENTIFIER || X_tile_residueMeta[positionOfMeta + i].w == CROWDER_IDENTIFIER) // repulsive crowder interaction
					#endif
					{
						LJ = crowderPotential(r);//__powf(6.0f/r,6);
					}
					else  // normal LJ interaction
					{
				#endif

				#if METADATA_MEMORY == TEXTURE_MEM
					int ijX(rint( AA_COUNT*yresiduem.x + meta.x));
				#else
					int ijX(rint( AA_COUNT*yresiduem.x + X_tile_residueMeta[positionOfMeta + i].x));
				#endif
				//do the texture fetch first
				#if LJ_LOOKUP_METHOD == TEXTURE_MEM
					float Eij(lambda*(tex1Dfetch(LJTexture,ijX) - e0));
					//float Eij(lambda*(e0 - e0));
				#elif LJ_LOOKUP_METHOD == SHARED_MEM
					float Eij(lambda*(sharedLJ[ijX] - e0));
				#else  // __global__ or __constant__
					float Eij(lambda*(LJPotentialData[ijX] - e0));
				#endif

				#if METADATA_MEMORY == TEXTURE_MEM
					DH = dhPotential(yresiduem.y,meta.y,r);
				#else
					DH = dhPotential(yresiduem.y,X_tile_residueMeta[positionOfMeta + i].y,r);
				#endif

				#if COMPENSATE_KERNEL_SUM
				// compensate sum

				y = DH - c_dh;
				t = dh_subtotal + y;
				c_lj = (t-dh_subtotal)-y;
				dh_subtotal = t;
				#else
				dh_subtotal += DH;
				#endif



					// sigmaij is the average atomic radius determined by the van der waals radius in kim2008
//					float sigmaij = __fmul_rn(__fadd_rn(yresiduem.z,X_tile_residueMeta[i].z), 0.5f);
					#if METADATA_MEMORY == TEXTURE_MEM
						float sigmaij((yresiduem.z + meta.z) * 0.5f);
					#else
						float sigmaij((yresiduem.z + X_tile_residueMeta[positionOfMeta + i].z) * 0.5f);
					#endif

					float LJtmp(__powf(sigmaij/r,6.0f));
					//float r0 = sigmaij*1.122462048309372981433533049679f; //sigma*pow(2.0,(1.0/6.0));

					LJ = -4.0f*Eij*LJtmp*(LJtmp-1.0f);
					//LJ = __fmul_rn(-4.0f,__fmul_rn(Eij,__fmul_rn(LJtmp,__fadd_rn(LJtmp,-1.0f))));
					if (Eij>0.0f && r < (sigmaij*1.122462048309372981433533049679f))
					{
						LJ = -LJ + 2.0f*Eij;
					}

				#if REPULSIVE_CROWDING
				}  // end conditional branch for LJ or repulsive short-range energy
				#endif
				//sharedmem_results[tx] += (LJ * LJ_CONVERSION_FACTOR) + (DH * DH_CONVERSION_FACTOR * 1.602176487f * 1.602176487f ); // moved outside loop

				#if COMPENSATE_KERNEL_SUM
				y = LJ - c_lj;
				t = lj_subtotal + y;
				c_lj = (t-lj_subtotal)-y;
				lj_subtotal = t;
				#else
				lj_subtotal += LJ;
				#endif

			} // if !(X_tile_residuePositions.w == Y_tile_residuePositions.w || X_tile_residuePositions.w < CROWDER_IDENTIFIER )
			sharedmem_results[tx] = (lj_subtotal * LJ_CONVERSION_FACTOR) + (dh_subtotal * DH_CONVERSION_FACTOR * 1.602176487f * 1.602176487f);
		} // for i = 0..Bdx

		#if CULL_LOWER_TRIANGULAR_SUM
		if (bx==by)
			sharedmem_results[tx] = sharedmem_results[tx]*0.5f;
		#else
		sharedmem_results[tx] *= 0.5f
		#endif //  CULL_LOWER_TRIANGULAR_SUM
	} // if !padding

	__syncthreads();  // all threads must have computed values so that they can be read concurrently from shared memory

	// do reduction, see SDK reduction example
	// FULLY unrolled reduction

	// do reduction in shared mem
    if (blockDim.x >= 512) { if (tx < 256) { sharedmem_results[tx] += sharedmem_results[tx + 256]; } __syncthreads(); }
    if (blockDim.x >= 256) { if (tx < 128) { sharedmem_results[tx] += sharedmem_results[tx + 128]; } __syncthreads(); }
    if (blockDim.x >= 128) { if (tx <  64) { sharedmem_results[tx] += sharedmem_results[tx + 64];  } __syncthreads(); }

	if (tx < 32)
	{
		if (blockDim.x >=  64) { sharedmem_results[tx] += sharedmem_results[tx + 32]; }
		if (blockDim.x >=  32) { sharedmem_results[tx] += sharedmem_results[tx + 16]; }
		if (blockDim.x >=  16) { sharedmem_results[tx] += sharedmem_results[tx +  8]; }
		if (blockDim.x >=   8) { sharedmem_results[tx] += sharedmem_results[tx +  4]; }
		if (blockDim.x >=   4) { sharedmem_results[tx] += sharedmem_results[tx +  2]; }
		if (blockDim.x >=   2) { sharedmem_results[tx] += sharedmem_results[tx +  1]; }
	}

	// write results of the calculation back to global memory, only thread 0 has the correct answer
	// only one warp processes < 32, therefore no syncthreads required here
	if (tx == 0)
		result[bx*gridDim.x+by] = sharedmem_results[tx];

	return;
};



/*__global__ void E_TiledLoopKernel(float4 * residuePositions, float4 * residueMeta, int * residueCount, int * moleculePositions, int * moleculeCount, float* LJPotentialData, float* result)
{
	// alternative method
	// we can also double buffer, will be better for small problem sizes
	__shared__ float sharedmem_results[TILE_DIM];	// array of temp results to sum, one per thread run.
	__shared__ float4 X_tile_residuePositions[TILE_DIM];
	__shared__ float4 X_tile_residueMeta[TILE_DIM];


	#if LJ_LOOKUP_METHOD == SHARED_MEM_LJ
		__shared__ float sharedLJ[AA_COUNT*AA_COUNT];
		for (int i=0;i<AA_COUNT*AA_COUNT;i+=AA_COUNT)
			memcpy(&(sharedLJ[i]),&(LJPotentialData[i]),sizeof(float)*AA_COUNT);
	#endif

	//__shared__ float sharedmem_results[TILE_DIM];	// array of temp results to sum, one per thread run.
	sharedmem_results[tx] = 0.0f;

	// save the local residue
	float4 yresiduep = residuePositions[bx*TILE_DIM+tx];
	float4 yresiduem = residueMeta[bx*TILE_DIM+tx];

	#define gx gridDim.x
	float mulfactor(0.5f);

	for (int address = bx*TILE_DIM; address <  TILE_DIM*gx; address+=TILE_DIM)
	{
		// copy residues from GRAM to shared memory tiles
		//__shared__ float4 X_tile_residuePositions[TILE_DIM];
		//__shared__ float4 X_tile_residueMeta[TILE_DIM];

		X_tile_residuePositions[tx] = residuePositions[address+tx];
		X_tile_residueMeta[tx]      = residueMeta[address+tx];

		// make sure no interactions are counted twice, when the preloaded tile address == block address * 1/2 to remove dupes
		// placed before the syncthreads because it does not rely on shared mem data

		__syncthreads();

		if (yresiduep.w > PADDER_IDENTIFIER) // the residue in question is not a padding element
		{
			for (unsigned int i = 0; i < TILE_DIM; i++)
			{
				if (yresiduep.w == X_tile_residuePositions[i].w || X_tile_residuePositions[i].w < 0.0f )		// same molecule || padding residue
				{
					// Xpos.w == Ypos.w means that they are the same molecule
					// if either is -2 then its a padding residue and must be ignored
					// DO NOTHING
				}
				else
				{
					//m.x = aminoacid index
					//m.y = charge
					//m.z = vdw radius
					//m.w = crowder if == CROWDER_IDENTIFIER

					float r = length(yresiduep,X_tile_residuePositions[i]) + EPS;  // add eps so that r is never 0, can happen in a collision
					float LJ = 0.0f;
					float DH = 0.0f;
					#if REPULSIVE_CROWDING
					if (yresiduem.w == CROWDER_IDENTIFIER || X_tile_residueMeta[i].w == CROWDER_IDENTIFIER) // repulsive crowder interaction
					{
						LJ = crowderPotential(r);//__powf(6.0f/r,6);
					}
					else  // normal LJ interaction
					{
					#endif


					//DH = yresiduem.y * X_tile_residueMeta[i].y  * __expf(-r/Xi) / r;
					DH = dhPotential(yresiduem.y,X_tile_residueMeta[i].y,r);


					int ijX =  rint( yresiduem.x + AA_COUNTf*X_tile_residueMeta[i].x);

					#if LJ_LOOKUP_METHOD == TEXTURE_MEM_LJ
						float Eij = lambda*(tex1Dfetch(LJTexture,ijX) - e0);
					#elif LJ_LOOKUP_METHOD == SHARED_MEM_LJ
						float Eij = lambda*(sharedLJ[ijX] - e0);
					#else  // __global__ or __constant__
						float Eij = lambda*(LJPotentialData[ijX] - e0);
					#endif


					// sigmaij is the average atomic radius determined by the van der waals radius in kim2008
					float sigmaij = (yresiduem.z + X_tile_residueMeta[i].z) * 0.5f;
					float LJtmp = __powf(sigmaij/r,6.0f);
					//float r0 = sigmaij*1.122462048309372981433533049679f; //sigma*pow(2.0,(1.0/6.0));
					//LJ = -4.0f*Eij*LJtmp*(LJtmp-1.0f);
					LJ = __fmul_rn(-4.0f,__fmul_rn(Eij,__fmul_rn(LJtmp,__fadd_rn(LJtmp,-1.0f))));
					if (Eij>0.0f && r < __fmul_rn(sigmaij,1.122462048309372981433533049679f))
					{
						LJ = -LJ + 2.0f*Eij;
					}

					#if REPULSIVE_CROWDING
					}  // end conditional branch for LJ or repulsive short-range energy
					#endif
					sharedmem_results[tx] += ((LJ * LJ_CONVERSION_FACTOR) + (DH * DH_CONVERSION_FACTOR * 1.602176487f * 1.602176487f ) )* mulfactor;
				} // if !(X_tile_residuePositions.w == Y_tile_residuePositions.w || X_tile_residuePositions.w < CROWDER_IDENTIFIER )
			} // for (unsigned int i = 0; i < TILE_DIM; i++)
		} // if (yresiduep.w > PADDER_IDENTIFIER) // the residue in question is not a padding element
		// sync to load next tile
		mulfactor = 1.0f;  // only the first set will be dupicated, so multiply by 1 each time now onwards
		__syncthreads();
	} // for (unsigned int address = bx*gx; address <  TILE_DIM*gx; address+=TILE_DIM)


	// do reduction, see SDK reduction example
    // FULLY unrolled reduction
    if (TILE_DIM >= 512) { if (tx < 256) { sharedmem_results[tx] += sharedmem_results[tx + 256]; } __syncthreads(); }
    if (TILE_DIM >= 256) { if (tx < 128) { sharedmem_results[tx] += sharedmem_results[tx + 128]; } __syncthreads(); }
    if (TILE_DIM >= 128) { if (tx <  64) { sharedmem_results[tx] += sharedmem_results[tx +  64]; } __syncthreads(); }
    if (tx < 32)
    {
        if (TILE_DIM >=  64) { sharedmem_results[tx] += sharedmem_results[tx + 32]; }
        if (TILE_DIM >=  32) { sharedmem_results[tx] += sharedmem_results[tx + 16]; }
        if (TILE_DIM >=  16) { sharedmem_results[tx] += sharedmem_results[tx +  8]; }
        if (TILE_DIM >=   8) { sharedmem_results[tx] += sharedmem_results[tx +  4]; }
        if (TILE_DIM >=   4) { sharedmem_results[tx] += sharedmem_results[tx +  2]; }
        if (TILE_DIM >=   2) { sharedmem_results[tx] += sharedmem_results[tx +  1]; }
    }

	// write results of the calculation back to global memory, only thread 0 has the correct answer
    if (tx == 0) result[bx] = sharedmem_results[tx];
	return;
};
*/

// performs a parallel sum on the data in values on the gpu with dimentions blockDim.x*blockDim.y
__global__ void parallelSum_kernel(float * values)
{
	#if REDUCTION_BLOCK_DIM > 256
		if (REDUCTION_BLOCK_DIM == 512) { if (tx < 256) { values[tx] += values[tx + 256]; } __syncthreads(); }
	#endif
	#if REDUCTION_BLOCK_DIM > 128
		if (REDUCTION_BLOCK_DIM >= 256) { if (tx < 128) { values[tx] += values[tx + 128]; } __syncthreads(); }
	#endif
	#if REDUCTION_BLOCK_DIM > 64
		if (REDUCTION_BLOCK_DIM >= 128) { if (tx <  64) { values[tx] += values[tx +  64]; } __syncthreads(); }
	#endif
		if (tx < 32)
		{
			if (REDUCTION_BLOCK_DIM >=  64) { values[tx] += values[tx + 32]; }
			if (REDUCTION_BLOCK_DIM >=  32) { values[tx] += values[tx + 16]; }
			if (REDUCTION_BLOCK_DIM >=  16) { values[tx] += values[tx +  8]; }
			if (REDUCTION_BLOCK_DIM >=   8) { values[tx] += values[tx +  4]; }
			if (REDUCTION_BLOCK_DIM >=   4) { values[tx] += values[tx +  2]; }
			if (REDUCTION_BLOCK_DIM >=   2) { values[tx] += values[tx +  1]; }
		}
};

__global__ void rotateMolecule_kernel (float4 *residuePositions, int *startPosition, int *length, float4* rotationVector, float4 *center)
{
	// cache frequently used variables
	__shared__ float4 centerV;
	__shared__ float4 q;

	centerV.x = center->x;
	centerV.y = center->y;
	centerV.z = center->z;

	// create rotation quaternion
	q.w = cos(rotationVector->w*0.5f);
	q.x = rotationVector->x*sin(rotationVector->w*0.5f);
	q.y = rotationVector->y*sin(rotationVector->w*0.5f);
	q.z = rotationVector->z*sin(rotationVector->w*0.5f);

	__syncthreads();

	// made it a unit quaternion
	float rq_len = rsqrt(q.x*q.x + q.y*q.y + q.z*q.z + q.w*q.w);

	__syncthreads();

	q.w = q.w*rq_len;
	q.x = q.x*rq_len;
	q.y = q.y*rq_len;
	q.z = q.z*rq_len;

	__syncthreads();

	//rotate the molecule using the quaternion made earlier

	#define vx (residuePositions[bx*blockDim.x+tx+startPosition[0]].x - centerV.x)
	#define vy (residuePositions[bx*blockDim.x+tx+startPosition[0]].y - centerV.y)
	#define vz (residuePositions[bx*blockDim.x+tx+startPosition[0]].z - centerV.z)

	if (bx*blockDim.x+tx < length[0])  // not addressing out of bounds
	{
		residuePositions[bx*blockDim.x+startPosition[0]+tx].x = ( vx*(q.w*q.w+q.x*q.x-q.y*q.y-q.z*q.z) + vy*(q.x*q.y-q.w*q.z)*2.0f  		+ vz*(q.w*q.y+q.x*q.z)*2.0f )        	 + centerV.x;
		residuePositions[bx*blockDim.x+startPosition[0]+tx].y = ( vx*(q.w*q.z+q.x*q.y)*2.0f   		 + vy*(q.w*q.w-q.x*q.x+q.y*q.y-q.z*q.z) + vz*(q.y*q.z-q.w*q.x)*2.0f ) 		     + centerV.y;
		residuePositions[bx*blockDim.x+startPosition[0]+tx].z = ( vx*(q.x*q.z-q.w*q.y)*2.0f   		 + vy*(q.y*q.z+q.w*q.x)*2.0f   		    + vz*(q.w*q.w-q.x*q.x-q.y*q.y+q.z*q.z) ) + centerV.z;
	}

};

__global__ void translateMolecule_kernel (float4 *residuePositions, int *startPosition, int *length, float4* translationVector, float4 *center)
{
	if (bx*blockDim.x+tx < length[0])  // not addressing out of bounds
	{
		residuePositions[bx*blockDim.x+startPosition[0]+tx].x = (residuePositions[bx*blockDim.x+startPosition[0]+tx].x + translationVector->x);
		residuePositions[bx*blockDim.x+startPosition[0]+tx].y = (residuePositions[bx*blockDim.x+startPosition[0]+tx].y + translationVector->y);
		residuePositions[bx*blockDim.x+startPosition[0]+tx].z = (residuePositions[bx*blockDim.x+startPosition[0]+tx].z + translationVector->z);
	}

	float3 tmpCenter;
	tmpCenter.x = center->x + translationVector->x;
	tmpCenter.y = center->y + translationVector->y;
	tmpCenter.z = center->z + translationVector->z;

	__syncthreads();

	center->x = tmpCenter.x;
	center->y = tmpCenter.y;
	center->z = tmpCenter.z;
};




#endif //USING CUDA
#endif
