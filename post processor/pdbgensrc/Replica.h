#ifndef REPLICA_H_
#define REPLICA_H_
#include "vector3f.h"
#include "definitions.h"
#include "Quaternion.h"
#include <fstream>
#include <vector>
#include "Molecule.h"
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
//#include "TorsionalLookupMatrix.h"
#include <pthread.h>
//#include "rngfunctions.h"


using namespace std;

#if USING_CUDA

#include <cutil.h>
#include <cuda.h>
#include <builtin_types.h>

#include "cudaExterns.h"

#endif  // USING_CUDA

class Replica
{
	public:
	Replica();
	~Replica();
	Replica(const Replica& r);
	void setAminoAcidData(AminoAcids a);
	void reserveConiguousMoleculeArray(int size);
	void initRNGs();
	void freeRNGs();
	void copy(const Replica& r); 						// assignment copy, duplicate it.
	//bool checkCollisions();							// check that molecules dont collide (for loading and explicit moves/rotates)
	//bool translateMolecule(int moleculeIndex);
	//bool rotateMolecule(int moleculeIndex);

	void exchangeReplicas(Replica &r);

	float f(double& surfaceAccessableSolventAreaRatio);
	float phi(Residue& i, Residue& j);
	double E();
	double Ef();

	float E_sse(); // using sse
	float E(Molecule *a,Molecule *b);   // use for fraction bound calcs
	float Eopt();  // cuda like optimisations to E
	float Eorig();

	// the search method that mutates the replicas molecules.
	void MCSearch(int steps);//, Replica * replaceReplica);

	bool savePDB(const char *filename);
	void saveAsSinglePDB(const char *filename);

	int loadMolecule(const char *pdbfilename);
	int loadMolecule(const char *pdbfilename, Vector3f position, Vector3double rotationAxis, double rotationAmount);

	AminoAcids aminoAcids;
	//TorsionalLookupMatrix torsions;

	Residue * contiguousResidues;
	int contiguousResiduesSize;
	Molecule * molecules;
	int moleculeCount;
	int moleculeArraySize;
	int residueCount;

	float boundingValue;
	void setBoundingValue(float value);
	float fullPotentialLoop(const uint y, const uint x, const uint cacheTile, double * LJacc, double * DHacc );
	float crowderLoop(const uint y, const uint x, const uint cacheTile, double * LJacc);



	float temperature;
	float temperatureBeta;
	short label;
	double potential;
	float E_DH;
	float E_LJ;
	float min_rotation;
	float min_translation;

	// rngs for this object
	gsl_rng * replacePRN;		// random number generator, replace if in/not in the boltzman dist
	gsl_rng * rng_translate;	// the translation vector rngs
	gsl_rng * rng_translateAmount;	// the translation vector rngs
	gsl_rng * rng_rotate;		// the rotation vector rng
	gsl_rng * rng_rotateAmount;	// the translation vector rngs
	gsl_rng * rng_moleculeSelection; // molecule selection rng
	gsl_rng  * RERng;  // used for replica exchange
	gsl_rng  * MCRng;  // used to determine what change to make
	gsl_rng  * MCKbRng;	// boltzmann acceptance rng in MC

	int threadCount;
	int nonCrowderCount;
	int maxMoleculeSize;

	int accept;
	int acceptA;
	int reject;
	int totalAcceptReject;
	int totalAccept;

	bool recordCrowders;


	int boundSamples;
	int samplesSinceLastExchange;

	int totalBoundSamples;
	int totalSamples;

	bool sample(SimulationData *data, int current_step, float boundEnergyThreshHold, pthread_mutex_t *writeFileMutex);
	float geometricDistance(Molecule *a,Molecule *b, float cutoff);
	void save(char* filename);

	float setTemperatureBeta() { temperatureBeta = temperature/300.0f; return temperatureBeta;};
	float getTemperatureBeta() { return temperatureBeta; };



	int countpairs();
	int paircount;
	float acc_err;

#if USING_CUDA

	float4 *device_float4_residuePositions; // x,y,z,w == position.x,position.y,position.z,id
	float4 *device_float4_residueMeta;		// x,y,z,w == index,charge,vdwr,temperature


	/*float *d_x;
	float *d_y;
	float *d_z;
	int *d_id;*/

	float4 *host_float4_residuePositions;
	float4 *host_float4_residueMeta;

	/*float *h_x;
	float *h_y;
	float *h_z;
	int *h_id;*/


	float *device_LJPotentials;

	int *device_residueCount;
	int *device_moleculeCount;
	int *device_moleculeStartPositions;		// position offset where each molecule starts in the above arrays
	int *host_moleculeStartPositions; 		// tracker for the above on the hose side so that individual molecules can be updated without additional calculation

	int blockSize;
	int dataSetSize;	// size of the data on gpu including padding residues
	int paddedSize;		// size taken in gpu memory because of alignment padding
	int gridSize;
	int resultSize;
	int sharedMemSize;

	bool replicaIsOnDevice;


	#if CUDA_MC
		int *device_moleculeLenghts;
		float4 *device_moleculeCenters;				// rotational centers of each molecule

		float4 *device_translationVector;
		float4 *device_reverseTranslationVector;
		float4 *device_rotationVector;  // vector(x,y,z)|amount(w)
		float4 *device_reverseRotationVector;
		float4 *host_translationVector;
		float4 *host_reverseTranslationVector;
		float4 *host_rotationVector;  // vector(x,y,z)|amount(w)
		float4 *host_reverseRotationVector;

		bool rotateOnDevice(int moleculeID, Vector3f vector, float amount);
		bool rotateOnDeviceM(int moleculeID, RotationalMatrix);
		bool translateOnDevice(int moleculeID, Vector3f translation);
		bool lastMutationWasATranslate;
		Vector3f reverseTranslate;
		RotationalMatrix reverseRotationalMatrix;
		int lastMutatedMolecule;
		void cudaRollbackMutation();
		// do all the mc functions on the device at once
		void MCSearchOnDevice(int steps);
	#endif

	void ReplicaDataToDevice();	// copy the replica to the device
	void initialiseDeviceLJPotentials();
	float * getDeviceLJPotentials() { return device_LJPotentials; };
	void setDeviceLJPotentials(float * ljp) { device_LJPotentials = ljp; }
	void ReplicaDataToHost();
	void UpdateDeviceData();
	void MoleculeDataToDevice(int MoleculeID); // update a molecule on the device
	void MoleculeDataToDeviceAsync(int MoleculeID); // update a molecule on the device
	double EonDevice();
	void EonDeviceAsync();
	void setLJpotentials(float *ljp);
	int getBlockSize();
	void setBlockSize(int blockSize);
	int getDataSetSize();
	void setDataSetSize(int dataSetSize);
	int getPaddedSize();
	void setPaddedSize(int paddedSize);
	int getGridSize();
	void setGridSize(int gridSize);
	int getResultSize();
	void setResultSize(int resultSize);


	void FreeDevice();
	void freeLJpotentials();


	#if CUDA_STREAMS
		cudaStream_t cudaStream;
		float *device_kernelResult;  	// holds the result of a kernel on the device
		float *kernelResult;			// as above but on host, sync cudaStream to make them the same
		void ReserveSumSpace();	// space on host and device to allow for the sum of kernel returns
		void FreeSumSpace();		// frees the above
		float SumGridResults();	//sums the grid returned by the potential kernel
		uint lastMutationIndex;
		Molecule savedMolecule;
		float oldPotential;
		float newPotential;
		//split the functions to allow for latency hiding and overlapping calls
		void MCSearchMutate();
		void MCSearchEvaluate();
		void MCSearchAcceptReject();

#endif



#endif

#if INCLUDE_TIMERS
	//timers for profiling the cuda functions
	uint replicaMCTimer;
	uint replicaMCCounter;
	uint replicaToGPUTimer;
	uint replicaToGPUCounter;
	uint replicaToHostTimer;
	uint replicaToHostCounter;
	uint replicaUpdateGPUTimer;
	uint replicaUpdateGPUCounter;
	uint replicaKernelTimer;
	uint replicaKernelCounter;
	uint replicaECUDATimer;
	uint replicaECUDACounter;
	uint replicaMoleculeUpdateTimer;
	uint replicaMoleculeUpdateCounter;
	uint replicaDeviceMCTimer;
	uint replicaDeviceMCCounter;
	uint initGPUMemoryTimer;
	uint initGPUMemoryCounter;
	uint replicaEHostTimer;
	uint replicaEHostCounter;
	bool timersInit;
#endif

	void initTimers();
	void printTimers();

	Vector3f createNormalisedRandomVector(gsl_rng * r);
	Vector3double createNormalisedRandomVectord(gsl_rng * r);

	Replica * GLreplica;


private:
	size_t molecules_max_size;
	void rotate(const int m, const double step);
	void translate(const int m, const float step);
};

//void bindLJtoTexture(int GPUID, float *ljp);
//void FreeLJfromTexture(int GPUID);


#endif /*REPLICA_H_*/
