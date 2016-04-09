#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

// Constant Data Input Files

#define AMINOACIDDATASOURCE 	"data/AminoAcids"
#define TORSIONALPAIRDATA 		"data/torsional_pair_potentials"

// global constants

#define K_b 		1.3806504240e-23f	// boltzmann constant  J/K
#define E_charge 	1.602176487e-19f	// elemental charge in coulombs
#define Angstrom	1e-10f				// one angstrom
#define D_water	  	80.0f				// dialectric constant for water
#define N_A			6.02214179e23f		// avogadros number
#define Rgas 		8.314472f			// universal gas constant
#define PI 			float(M_PIl)
#define KBTConversionFactor  (1.0f/(294.0f*Rgas/4184.0f))  // use 294K for verification
#define K_bTtoKcalmol  KBTConversionFactor
#define KcalmoltoK_bt  (1.0f/K_bTtoKcalmol)
#define BTU_to_J             0.948f
#define Xi 10.0f
#define LJ_CONVERSION_FACTOR 0.3507221006079f  // 1/(KBTConversionFactor)^2
#define DH_CONVERSION_FACTOR BTU_to_J //10e7f //2.0f
#define CF KBTConversionFactor

//#define LJPDSOURCE 	"data/pair_potentials_BT"
//#define e0			0.0372f		// K_b T : experimental fitted value, see [K & B 2008], BT model
//#define lambda		1.700f		// scalar  : experimental fitted value, see [K & B 2008]
//#define LJPDSOURCE 	"data/pair_potentials_MJ1999"
//#define e0			0.0328f		// K_b T : experimental fitted value, see [K & B 2008], MJ1999 model
//#define lambda		1.480f		// scalar  : experimental fitted value, see [K & B 2008]
#define LJPDSOURCE 		"data/pair_potentials_MJ1996"
#define e0				-2.27f			// K_b T : experimental fitted value, see [K & B 2008], MJ1996 model
#define lambda			0.159f			// scalar  : experimental fitted value, see [K & B 2008]
#define AA_COUNT 		20
#define AA_COUNTf 		float(AA_COUNT)
#define LJArraySize		(sizeof(float)*AA_COUNT*AA_COUNT)

// pseudobond constants

#define K_spring 	378 			// 378 kcal / (mol.A^2)
#define R0			3.81f			// reference bond length in atomic units

// pseudo angle constants

#define GammaAngle 		0.1f		// mol/kcal
#define EpsilonAlpha	4.3f		// kcal/mol
#define ThetaAlpha		1.6f 		// radians
#define ThetaBeta		2.27f 		// radians
#define KAlpha			106.4f		// kcal/(mol rad)^2
#define KBeta			26.3f		// kcal/(mol rad)^2
// membrane interaction constants

#define Zreference  20.0f
#define ephi0 		-30.0f 			//mV check units
#define invDebye  	1e9f		    // inverse debye length water = 1/(10 angtrom)
#define kappa 		invDebye

// Simulation Values
#define EPS 1e-38f
#define BOUNDING_RADIUS		100.0f   // angstrom. size of sphere for spherical boundary conditions
#define BOUNDING_VALUE		150.0f   // size of box edge for periodic conditions
#define BOUNDING_SPHERE  	0
#define PERIODIC_BOUNDRY 	1
#define BOUNDING_METHOD	  	PERIODIC_BOUNDRY
#define MEMBRANE_PRESENT 0 // use membrane calculations, molecule 0 is the embedded molecule in this case.

#define REPULSIVE_CROWDING 	0

#define REPLICA_COUNT 						20
#define MC_STEPS 							5000
#define HIGHTEMP 							500.0f
#define LOWTEMP 							250.0f
#define REMC_STEPS							10
#define STEPS_BEFORE_SAMPLE 				MC_STEPS
#define SAMPLE_FREQ							MC_STEPS
#define INITIAL_TRANSLATIONAL_STEP			0.5f // angstrom
//#define MIN_TRANSLATION						0.25f
#define INITIAL_ROTATIONAL_STEP				0.2f // rads
//#define MIN_ROTATION						0.05f
#define BOUND_ENERGY_VALUE					-1.1844f	// 2k_bT in kcal/mol at 294K
#define translate_rotate_bernoulli_bias 	0.5f // translations / mutations (rotations+translations) ratio
#define GEOMETRIC_CUTOFF					8.0f

#define PADDER_IDENTIFIER	-2.0f	//id values for identifying a padding residue
#define CROWDER_IDENTIFIER	-1.0f 	//id values for identifying a crowder residue
#define PADDER_AATYPE 		20.0f

// code things

#define PRINT_RE_STEP_COUNT	   1
#define PRINT_RE_STEP_RESULTS  2
#define PRINT_MC_STEP_COUNT    3
#define PRINT_MC_MUTATIONS     4
#define OUTPUT_LEVEL		   1

#ifdef EnableOPENGL			//set here to override the makefile
	#define GLVIS 1
#else
	#define GLVIS 0
#endif

#define GL_AXES 			1
#define THREADING 			1
#define THREAD_COUNT		1
#define STREAM_COUNT		1
#define INCLUDE_TIMERS		0	//include timers
//#define CHECKPOINTING		1
#define CHECKPOINTFREQUENCY 10000

#ifdef EnableCUDA			//set here to override the makefile
	#define USING_CUDA 1
#else
	#define USING_CUDA 0
#endif

#define COMPENSATE_KERNEL_SUM	  0

// override for the above definition
//#define USING_CUDA 1

#if USING_CUDA
	#define PERFORM_GPU_AND_CPU_E  0  // do both cpu and gpu sums for comparison

	#define GPU_COUNT          		1
	#define USING_TILE_KERNEL  		1
	#define CUDA_E      	  		1	// use cuda to do the energy summations (vs cpu if not included)
	#define PARALLEL_SUM			0

	#define SHARED_MEM				0
	#define CONST_MEM				1
	#define GLOBAL_MEM				2
	#define TEXTURE_MEM				3

	#define LJ_LOOKUP_METHOD        TEXTURE_MEM // type of memory to use for LJ lookups
	#define METADATA_MEMORY			SHARED_MEM  // type of memory to use for residue metadata
	#define POSITIONDATA_MEMORY		SHARED_MEM  // type of memory to use for residue position data

	#define USE_LOOP_UNROLL	          0
	#define CULL_LOWER_TRIANGULAR_SUM 1		// only computer the upper triangular sum for the summations, saves 0.5*n(n-1) work
	#ifdef EnableStreams
		#define CUDA_STREAMS 1
	#else
		#define CUDA_STREAMS 0
	#endif


	#define CUDA_MC		               0	// use cuda to do the MC operations

	#define TILE_DIM 			64 // block size for looping kernel, must be a power of 2 <= 512
	#define REDUCTION_BLOCK_DIM		128   // block size for reductions

	#define BLOCK_SIZE 				22  // N*N <= 512
	#if BLOCK_SIZE <= 4
		#define SHARED_MEM_RESULT_SIZE	16 		// [nearest power of 2 > BLOCK_SIZE]^2
	#elif BLOCK_SIZE <= 5
		#define SHARED_MEM_RESULT_SIZE	32
	#elif BLOCK_SIZE <= 8
		#define SHARED_MEM_RESULT_SIZE	64
	#elif BLOCK_SIZE <= 11
		#define SHARED_MEM_RESULT_SIZE	128
	#elif BLOCK_SIZE <= 16
		#define SHARED_MEM_RESULT_SIZE	256
	#elif BLOCK_SIZE <= 22						//sqrt 512 == 22.something
		#define SHARED_MEM_RESULT_SIZE	512
	#endif
#endif // USING_CUDA 1

#if GLVIS
	void GlutDisplay();
#endif

#include <iostream>

class Replica;

struct SimulationData
{
	long index;
	Replica *replica;
	int replicaCount;
	int GPUID;
	int threads;
	int streams;
	int MCsteps;
	int REsteps;
	int sampleFrequency;
	int sampleStartsAfter;
	float bound;
	int * waitingThreadCount;
	int * conformationsBound;
	FILE * fractionBound;
	FILE * boundConformations;
};

struct argdata
{
	int threads;
	int streams;
	int gpus;
	int MCsteps;
	int REsteps;
	int replicas;
	float bound;
	int sampleFrequency;
	int sampleStartsAfter;
	char prependageString[256];
	char file[256];
	char logfile[256];
	bool inputFile;
	int nonCrowders;
	float temperatureMin;
	float temperatureMax;
	int pid;
	char checkpointfilename[256];
	bool resume;
	int currentStep;
};

#endif /*DEFINITIONS_H_*/
