#include "Replica.h"

using namespace std;

#if USING_CUDA || INCLUDE_TIMERS
	#include <cutil.h>
#endif

Replica::Replica()
{
	temperature = 300.0f;
	label = -1;
	potential = 0.0f;
	moleculeCount = 0;
	moleculeArraySize = 0;
	residueCount = 0;
	maxMoleculeSize = 0;
	nonCrowderCount = 0;
	accept = 0;
	acceptA = 0;
	reject = 0;
	totalAcceptReject = 0;
	totalAccept = 0;
	boundSamples = 0;
	samplesSinceLastExchange = 0;
	totalBoundSamples = 0;
	totalSamples = 0;
	paircount = 0;
	acc_err = 0;
	//min_rotation = MIN_ROTATION;
	//min_translation = MIN_TRANSLATION;
	
	#if INCLUDE_TIMERS
	timersInit = false;
	#endif

}

#if INCLUDE_TIMERS
void Replica::initTimers()
{
	//timers for profiling the cuda functions
	replicaMCTimer = 0;
	replicaMCCounter = 0;

	replicaToGPUTimer = 0;
	replicaToGPUCounter = 0;
	replicaUpdateGPUTimer = 0;
	replicaUpdateGPUCounter = 0;
	replicaKernelTimer = 0;
	replicaKernelCounter = 0;
	replicaECUDATimer = 0;
	replicaECUDACounter = 0;
	replicaMoleculeUpdateTimer = 0;
	replicaMoleculeUpdateCounter = 0;
	replicaDeviceMCTimer = 0;
	replicaDeviceMCCounter = 0;
	initGPUMemoryTimer = 0;
	initGPUMemoryCounter = 0;
	replicaEHostCounter = 0;


	CUT_SAFE_CALL( cutCreateTimer(&replicaMCTimer) );
	CUT_SAFE_CALL( cutCreateTimer(&replicaToGPUTimer) );
	CUT_SAFE_CALL( cutCreateTimer(&replicaToHostTimer) );
	CUT_SAFE_CALL( cutCreateTimer(&replicaUpdateGPUTimer) );
	CUT_SAFE_CALL( cutCreateTimer(&replicaKernelTimer) );
	CUT_SAFE_CALL( cutCreateTimer(&replicaMoleculeUpdateTimer) );
	CUT_SAFE_CALL( cutCreateTimer(&replicaECUDATimer) );
	CUT_SAFE_CALL( cutCreateTimer(&replicaDeviceMCTimer) );
	CUT_SAFE_CALL( cutCreateTimer(&initGPUMemoryTimer) );
	CUT_SAFE_CALL( cutCreateTimer(&replicaEHostTimer) );
	timersInit = true;
}
#endif

Replica::Replica(const Replica& r)
{
	Replica();
	label = r.label;
	temperature = r.temperature;
	moleculeCount = r.moleculeCount;
	residueCount = r.residueCount;
	aminoAcids = r.aminoAcids;
	molecules = new Molecule[moleculeCount];
	maxMoleculeSize = r.maxMoleculeSize;
	nonCrowderCount = r.nonCrowderCount;
	for (size_t m=0;m<moleculeCount;m++)
	{
		molecules[m] = r.molecules[m];
		cout << "molecule copy: " << &molecules[m] << " <- " << &r.molecules[m]<< endl;
	}
	#if USING_CUDA
		blockSize = r.blockSize;
		sharedMemSize = r.sharedMemSize;
		replicaIsOnDevice = false;
	#endif

	paircount = r.paircount;
	acc_err = 0;
	min_translation = r.min_translation;
	min_rotation = r.min_rotation;
	temperatureBeta = r.temperatureBeta;
}

void Replica::setAminoAcidData(AminoAcids a)
{
	aminoAcids = a;
}

Replica::~Replica()
{
	//aminoAcids = AminoAcids();
	#if INCLUDE_TIMERS
	if (timersInit)
	{
		CUT_SAFE_CALL( cutDeleteTimer(replicaMCTimer) );
		CUT_SAFE_CALL( cutDeleteTimer(replicaToGPUTimer) );
		CUT_SAFE_CALL( cutDeleteTimer(replicaToHostTimer) );
		CUT_SAFE_CALL( cutDeleteTimer(replicaUpdateGPUTimer) );
		CUT_SAFE_CALL( cutDeleteTimer(replicaKernelTimer) );
		CUT_SAFE_CALL( cutDeleteTimer(replicaMoleculeUpdateTimer) );
		CUT_SAFE_CALL( cutDeleteTimer(replicaECUDATimer) );
		CUT_SAFE_CALL( cutDeleteTimer(replicaDeviceMCTimer) );
		CUT_SAFE_CALL( cutDeleteTimer(initGPUMemoryTimer) );
		CUT_SAFE_CALL( cutDeleteTimer(replicaEHostTimer) )
	}
	#endif
	if (moleculeCount > 0)
		delete [] molecules;
}

int Replica::countpairs()
{
	paircount = 0;
	for (size_t mI=0; mI<moleculeCount;mI++)
		for (size_t mJ=mI+1; mJ<moleculeCount;mJ++)
			for (size_t mi=0; mi<molecules[mI].residueCount; mi++)
				for (size_t mj=0; mj<molecules[mJ].residueCount; mj++)
					paircount++;
	return paircount;
}

void Replica::exchangeReplicas(Replica &r)
{
	swap(label,r.label);
	swap(temperature,r.temperature);

	// for statistics
	swap(totalAcceptReject,r.totalAcceptReject);
	swap(totalAccept,r.totalAccept);
	swap(boundSamples,r.boundSamples);
	swap(samplesSinceLastExchange,r.samplesSinceLastExchange);
	swap(totalBoundSamples,r.totalBoundSamples);
	swap(totalSamples,r.totalSamples);


}

void Replica::copy(const Replica &r)
{
	label = r.label;
	temperature = r.temperature;
	aminoAcids = r.aminoAcids;
	residueCount = r.residueCount;
	maxMoleculeSize = r.maxMoleculeSize;
	boundingValue = r.boundingValue;
	nonCrowderCount = r.nonCrowderCount;

	if (molecules != NULL && moleculeCount != r.moleculeCount )
	{
		try
		{
			delete [] molecules;
		}
		catch ( char * str )
		{
		   cout << "Exception raised: delete [] molecules in Replica::copy()" << endl;
		}
	}

	molecules = new Molecule[r.moleculeCount];
	moleculeCount = r.moleculeCount;

	//the next bit is important because it makes sure that residues are contiguous in memory => better cpu performance
	contiguousResiduesSize = int(32.0f*ceil(float(residueCount)/32.0f));
	contiguousResidues = new Residue[contiguousResiduesSize];
	// contiguousResidues[residueCount..contiguousResiduesSize-1].aminoAcidIndex == PADDER_IDENTIFIER
	for (int i=residueCount; i<contiguousResiduesSize;i++)
	{
		contiguousResidues[i].aminoAcidIndex = int(PADDER_IDENTIFIER);
	}


	int rescount = 0;
	for (size_t m=0;m<moleculeCount;m++)
	{
		memcpy(contiguousResidues+rescount,r.molecules[m].Residues,r.molecules[m].residueCount*sizeof(Residue));

		molecules[m].Residues = contiguousResidues+rescount;
		molecules[m].residueCount = r.molecules[m].residueCount;
		molecules[m].translationalStep = r.molecules[m].translationalStep;
		molecules[m].rotationalStep = r.molecules[m].rotationalStep;
		molecules[m].AminoAcidsData = r.molecules[m].AminoAcidsData;
		molecules[m].center = r.molecules[m].center;
		molecules[m].position = r.molecules[m].position;
		molecules[m].rotation = r.molecules[m].rotation;
		molecules[m].moleculeRoleIdentifier = r.molecules[m].moleculeRoleIdentifier;
		rescount += r.molecules[m].residueCount;
	}
	potential = r.potential;
	#if USING_CUDA
	blockSize = r.blockSize;
	sharedMemSize = r.sharedMemSize;
	#endif
}

void Replica::reserveConiguousMoleculeArray(int size)
{
	if (moleculeCount != 0)
	{
		cout << " ** Molecules already exist. Cannot reserve size" << endl;
		return;
	}
	moleculeArraySize = size;
	molecules = new Molecule[size];
}

int Replica::loadMolecule(const char* pdbfilename)
{
	if (moleculeCount+1 > moleculeArraySize) //need a new array
	{
		moleculeArraySize += 3;
		Molecule *new_molecules = new Molecule[moleculeArraySize];
		if (moleculeCount > 0)
			memcpy(new_molecules,molecules,sizeof(Molecule)*moleculeCount);

		try	{	delete [] molecules; }  // free existing array before forgetting about it
		catch ( char * str ) 	{   cout << "Exception raised: delete [] molecules failed" << endl;}

		molecules = new Molecule[moleculeArraySize];
		memcpy(molecules,new_molecules,sizeof(Molecule)*moleculeCount);
		delete [] new_molecules;
	}

	molecules[moleculeCount].AminoAcidsData = aminoAcids;
	molecules[moleculeCount].index = moleculeCount;
	molecules[moleculeCount].initFromPDB(pdbfilename);
	residueCount += molecules[moleculeCount].residueCount;
	maxMoleculeSize = max(maxMoleculeSize,molecules[moleculeCount].residueCount);
	moleculeCount++;
	return moleculeCount-1;
}

int Replica::loadMolecule(const char* pdbfilename, Vector3f position, Vector3double rotationAxis, double rotationAmount)
{
	int i = loadMolecule(pdbfilename);
	rotationAxis.normalizeInPlace();
	molecules[i].setPosition(position);
	molecules[i].rotateQ(rotationAxis,rotationAmount);
	return i;
}

void Replica::initRNGs()
{
	unsigned long long seed = time (NULL);
	srand(time(NULL)+(label+1)*(label+1));

	rng_moleculeSelection = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (rng_moleculeSelection,random());

	rng_rotate = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (rng_rotate,random());

	rng_rotateAmount = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (rng_rotateAmount,random());

	rng_translate = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (rng_translate,random());

	rng_translateAmount = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (rng_translateAmount,random());

	MCRng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set (MCRng,random());

	MCKbRng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set (MCKbRng,random());
}

void Replica::freeRNGs()
{
	gsl_rng_free (rng_moleculeSelection);
	gsl_rng_free (rng_rotateAmount);
	gsl_rng_free (rng_translateAmount);
	gsl_rng_free (rng_rotate);
	gsl_rng_free (rng_translate);
	gsl_rng_free (MCRng);
	gsl_rng_free (MCKbRng);
}

Vector3f Replica::createNormalisedRandomVector(gsl_rng * r)
{
	//	float denominator = gsl_rng_max(r) - gsl_rng_min(r);
	//	Vector3f x = 2.0f*(Vector3f(gsl_rng_get(r)/denominator-0.5f,gsl_rng_get(r)/denominator-0.5f,gsl_rng_get(r)/denominator-0.5f));
	//Vector3f x = (Vector3f(gsl_rng_uniform(r)-0.5,0.5-gsl_rng_uniform(r),gsl_rng_uniform(r)-0.5)).normalize();
	//return x;
	return (Vector3f(gsl_rng_uniform(r)-0.5,0.5-gsl_rng_uniform(r),gsl_rng_uniform(r)-0.5)).normalize();
}

Vector3double Replica::createNormalisedRandomVectord(gsl_rng * r)
{
	//double denominator = gsl_rng_max(r) - gsl_rng_min(r);
	Vector3double x(gsl_rng_uniform(r)-0.5,0.5-gsl_rng_uniform(r),gsl_rng_uniform(r)-0.5);
	x.normalizeInPlace();
	return x;
}

#define _translate	0
#define _rotate 	1

// here for profiling
inline void Replica::rotate(const int m, const double rotateStep)
{
	molecules[m].rotateQ(createNormalisedRandomVectord(rng_rotate),rotateStep);
}

// here for profiling
void Replica::translate(const int m, const float translateStep)
{
	// pick directions to rotate in, do by generating a random normalised vector of length equal to INITIAL_TRANSLATIONAL_STEP
	//Vector3f oldPosition = molecules[moleculeNo].center;
	Vector3f translateVector = translateStep * createNormalisedRandomVector(rng_translate);


	// bounding sphere conditions
	#if BOUNDING_METHOD == BOUNDING_SPHERE
		// if the translate falls inside the bounding radius
		if ((molecules[m].center + translateVector).sumSquares() < boundingValue*boundingValue)
		{
			molecules[m].translate(translateVector);
		}
		// if the change causes the molecule to fall outside the bounding radius
		else
		{
			#if OUTPUT_LEVEL >= PRINT_MC_MUTATIONS
				cout << "  - Reject: Replica "<< label << "/Molecule " << moleculeNo <<  " : falls outside bounding sphere" << endl;
			#endif
		}
	#elif BOUNDING_METHOD == PERIODIC_BOUNDRY 	// periodic boundary conditions
		Vector3f newPosition = molecules[m].center + translateVector;
		newPosition.x = fmod(newPosition.x+boundingValue,boundingValue);
		newPosition.y = fmod(newPosition.y+boundingValue,boundingValue);
		newPosition.z = fmod(newPosition.z+boundingValue,boundingValue);
		molecules[m].setPosition(newPosition);
	#endif
}



void Replica::MCSearch(int steps)
{
	float translateStep = INITIAL_TRANSLATIONAL_STEP;
	double rotateStep = INITIAL_ROTATIONAL_STEP;

	bool lastOperationWasRotate;

	float oldPotential = potential;

	Molecule savedMolecule;
	savedMolecule.reserveResidueSpace(maxMoleculeSize);

	for (int step=0;step<steps;step++)
	{
		#if OUTOUTPUT_LEVEL >= PRINT_MC_STEP_COUNT
		cout << "Step: " << step << endl;
		#endif
		//moleculeNo = ((int)(1000.0*_gsl_qrng_get (qrng_moleculeSelection)))%((int)moleculeCount);
		uint moleculeNo = (int) gsl_rng_uniform_int(rng_moleculeSelection,moleculeCount);

		uint mutationType = gsl_ran_bernoulli (MCRng,translate_rotate_bernoulli_bias);

		// save the current state so we can role back if it was not a good mutation.
		savedMolecule.saveBeforeStateChange(&molecules[moleculeNo]);

		switch (mutationType)
		{
			case _rotate:
			{
				//molecules[moleculeNo].rotateQ(createNormalisedRandomVectord(rng_rotate),rotateStep);
				rotateStep = temperatureBeta * gsl_rng_uniform_pos(rng_rotateAmount) * INITIAL_ROTATIONAL_STEP + INITIAL_ROTATIONAL_STEP/2.0f;
				rotate(moleculeNo, rotateStep);

				#if OUTPUT_LEVEL >= PRINT_MC_MUTATIONS
				cout << "    Rotate: Replica "<< label << "/Molecule " << moleculeNo << endl;
				#endif
				break;
			}

			case _translate:
			{
				translateStep = temperatureBeta * gsl_rng_uniform_pos(rng_translateAmount) * INITIAL_TRANSLATIONAL_STEP + INITIAL_TRANSLATIONAL_STEP/2.0f;
				translate(moleculeNo, translateStep);

				#if OUTPUT_LEVEL >= PRINT_MC_MUTATIONS
				cout << "    Translate: Replica "<< label << "/Molecule " << moleculeNo << endl;
				#endif

				break;
			}
			default:
				break;
		}

#if CUDA_E
		// copy host data to device. so we can do the calculations on it.
		MoleculeDataToDevice(moleculeNo);

		double newPotential(EonDevice());  // sequential cuda call

	#if PERFORM_GPU_AND_CPU_E
		float cpu_e(E());
		float err = abs(cpu_e-newPotential)/abs(cpu_e);
		//acc_err += err;
		printf("%24.20f %24.20f %24.20f\n",cpu_e,float(newPotential),err);
	#endif
#else // only CPU calls
		float newPotential = E();
#endif

		float delta = newPotential - oldPotential;

		// accept change if its better.
		if (delta < 0.0)
		{
			potential = newPotential;
			oldPotential = potential;
			accept++;
			#if OUTPUT_LEVEL >= PRINT_MC_MUTATIONS
			cout << "  * Replace: Replica "<< label << "/Molecule " << moleculeNo <<  " : delta E = " << delta << " E = " << potential << endl;
			#endif
		}
		// accept change if it meets the boltzmann criteria, must be (kJ/mol)/(RT), delta is in kcal/mol @ 294K
		else if (gsl_rng_uniform(MCKbRng)<exp(-(delta*4184.0f)/(Rgas*temperature)))
		{
			potential = newPotential;
			oldPotential = potential;
			acceptA++;
			#if OUTPUT_LEVEL >= PRINT_MC_MUTATIONS
			cout << "  **Replace: Replica "<< label << "/Molecule " << moleculeNo <<  " : delta E = " << delta << " U < " << exp(-delta * 4.184f/(Rgas*temperature)) << " E = " << potential << endl;
			#endif
		}
		//reject
		else
		{
			#if OUTPUT_LEVEL >= PRINT_MC_MUTATIONS
				cout << "  - Reject: Replica "<< label << "/Molecule " << moleculeNo <<  " : delta E = " << delta << " E = " << oldPotential << endl;
			#endif

			reject++;
			molecules[moleculeNo].undoStateChange(&savedMolecule);

			#if CUDA_E
					MoleculeDataToDevice(moleculeNo); // you have to update the device again because the copy will be inconsistent
			#endif
			potential = oldPotential;
		}
	}
	delete [] savedMolecule.Residues;

}

bool Replica::savePDB(const char *filename) // saves multiple pdb files per replica
{
	char filenameExt[256];
	char tmp[64];
	for (size_t i=0;i<moleculeCount;i++)
	{
		strcpy (filenameExt,filename);
		sprintf (tmp,"%02d",int(i));
		strcat (filenameExt,tmp);
		strcat (filenameExt,".pdb");
		molecules[i].saveAsPDB(filenameExt);
	}
	return true;
}

void Replica::saveAsSinglePDB(const char *filename)
{
	FILE * output;
	output = fopen (filename,"w");
	fprintf(output,"REMARK %s \n",filename);
	fprintf(output,"REMARK potential: %0.10f \n",float(potential));
	fprintf(output,"REMARK temperature: %5.1f \n",temperature);

	for (size_t i=0;i<moleculeCount;i++)
	{
		fprintf(output,"REMARK Molecule: %d\n",int(i));
		fprintf(output,"REMARK Rotation relative to input Q(w,x,y,z): %f %f %f %f\n",molecules[i].rotation.w,molecules[i].rotation.x,molecules[i].rotation.y,molecules[i].rotation.z);
		fprintf(output,"REMARK Centriod position P(x,y,z): %f %f %f\n",molecules[i].center.x,molecules[i].center.y,molecules[i].center.z);
	}

	char chainId = 'A';
	int itemcount = 0;
	int lastSeqNo = 0;
	for (size_t m=0;m<moleculeCount;m++)
	{
		size_t i=0;
		while (i<molecules[m].residueCount)
		{
			itemcount++;
			fprintf(output,"ATOM  %5d %4s%C%3s %C%4d%C  %8.3f%8.3f%8.3f%6.2f%6.2f\n",itemcount,"CA",' ',aminoAcids.get(molecules[m].Residues[i].aminoAcidIndex).getSNAME(),chainId,molecules[m].Residues[i].resSeq,' ',molecules[m].Residues[i].position.x,molecules[m].Residues[i].position.y,molecules[m].Residues[i].position.z,1.0f,1.0f);
			lastSeqNo = molecules[m].Residues[i].resSeq;
			i++;
		}
		fprintf(output,"TER   %5d      %3s %C%4d \n",itemcount,aminoAcids.get(molecules[m].Residues[i-1].aminoAcidIndex).getSNAME(),chainId,lastSeqNo);
		chainId++;
		fflush(output);
	}
	fprintf(output,"END \n");
	fflush(output);
	fclose(output);
}
// simulation evaluations

inline float Replica::phi(Residue& i, Residue& j)
{
	// eps prevents division by zero
	float r = EPS + sqrt( (i.position.x-j.position.x)*(i.position.x-j.position.x)+(i.position.y-j.position.y)*(i.position.y-j.position.y)+(i.position.z-j.position.z)*(i.position.z-j.position.z));

	// first do the lennard jones pair interactions
	float LJ = 0.0f;
	float DH = 0.0f;
	#if REPULSIVE_CROWDING
	if(i.aminoAcidIndex==CROWDER_IDENTIFIER || j.aminoAcidIndex==CROWDER_IDENTIFIER)
	{
		//float ctmp = 6.0f/r;
		//LJ = ctmp*ctmp*ctmp*ctmp*ctmp*ctmp;
		LJ = 46656.0f / powf(r,6);
	}
	else
	{
	#endif
		float Eij = lambda*(aminoAcids.LJpotentials[i.aminoAcidIndex][j.aminoAcidIndex] - e0);

		// sigmaij is the average atomic radius determined by the van der waal radius in kim2008
		float sigmaij = 0.5f * ( i.vanderWaalRadius + j.vanderWaalRadius );
		//float r0 = sigmaij*1.122462048309372981433533049679f; //pow(2.0,(1.0/6.0));
		float sigT = sigmaij / r ;
		float LJtmp = sigT*sigT*sigT*sigT*sigT*sigT;
		if (Eij<0.0f)  // attractive pairs
		{
			LJ = 4.0f*abs(Eij)*( LJtmp*(LJtmp-1.0f));
		}
		else 	// repulsive pairs
		{
			//if (r<r0)
			if (r<sigmaij*1.122462048309372981433533049679f)
			{
				LJ = 4.0f*Eij*( LJtmp*(LJtmp-1.0f)) + 2.0f*Eij;
			}
			else // r >= r0
			{
				LJ = -4.0f*Eij*( LJtmp*(LJtmp-1.0f));
			}
		}
		// do the debye huckel long range intereactions
		float DH_constant_component =  1.602176487f * 1.602176487f * DH_CONVERSION_FACTOR;
		DH = i.electrostaticCharge * j.electrostaticCharge * expf(-r/Xi) / r * DH_constant_component;

	#if REPULSIVE_CROWDING
	}
	#endif
	return (LJ * LJ_CONVERSION_FACTOR + DH); // convert to kcal/mol
}

inline float crowderPairPotential(const float r)
{
	return powf(6.0f/r,6);
}

// evaluate the conformations  energy
double Replica::E()
{
#if INCLUDE_TIMERS
	replicaEHostCounter++;
	CUT_SAFE_CALL( cutStartTimer(replicaEHostTimer) );
#endif

	double epotential = 0.0f;
    double LJAccumulator = 0.0f;
	double DHAccumulator = 0.0f;
	double DH_constant_component =  DH_CONVERSION_FACTOR * 1.602176487f * 1.602176487f ;

	#if COMPENSATE_KERNEL_SUM
	double c_lj(0.0f);
	double c_dh(0.0f);

	/*
	function kahanSum(input)
	 var potential = 0
	 var c = 0
	 for i = 0 to blockdim-1
	  y = Pij - c
	  t = potential + y
	  c = (t - potential) - y
	  potential = t
	 next i
	return sum
	*/
	#endif

	#define iRes molecules[mI].Residues[mi]
	#define jRes molecules[mJ].Residues[mj]
	// between each residue do the LJ and DH potentials
	// do summation for LJ and DH potentials
	for (size_t mI=0; mI<moleculeCount;mI++)
	{
		#if REPULSIVE_CROWDING
		if (molecules[mI].moleculeRoleIdentifier == CROWDER_IDENTIFIER)
		{
			for (size_t mJ=mI+1; mJ<moleculeCount;mJ++)
			{
				for (size_t mi=0; mi<molecules[mI].residueCount; mi++)
				{
					for (size_t mj=0; mj<molecules[mJ].residueCount; mj++)
					{
						//r = EPS + sqrtf( (iRes.position.x-jRes.position.x)*(iRes.position.x-jRes.position.x)+(iRes.position.y-jRes.position.y)*(iRes.position.y-jRes.position.y)+(iRes.position.z-jRes.position.z)*(iRes.position.z-jRes.position.z));
						//LJAccumulator += (46656.0 / pow(r,6));
						//LJAccumulator += crowderPairPotential(scalar_distance(molecules[mI].Residues[mi].position,molecules[mJ].Residues[mj].position)+ EPS);
						#if COMPENSATE_KERNEL_SUM
							double 	y(crowderPairPotential(scalar_distance(molecules[mI].Residues[mi].position,molecules[mJ].Residues[mj].position)+ EPS) - c_lj);
							double t(LJAccumulator + y);
							c_lj = (t-LJAccumulator)-y;
							LJAccumulator = t;
						#else
							LJAccumulator += crowderPairPotential(scalar_distance(molecules[mI].Residues[mi].position,molecules[mJ].Residues[mj].position)+ EPS);
						#endif

					}
				}
			}
		}
		else
		#endif
		for (size_t mJ=mI+1; mJ<moleculeCount;mJ++)
		{

			#if REPULSIVE_CROWDING
			if (molecules[mJ].moleculeRoleIdentifier == CROWDER_IDENTIFIER)
			{
				for (size_t mi=0; mi<molecules[mI].residueCount; mi++)
				{
					for (size_t mj=0; mj<molecules[mJ].residueCount; mj++)
					{
						#if COMPENSATE_KERNEL_SUM
							double 	y(crowderPairPotential(scalar_distance(molecules[mI].Residues[mi].position,molecules[mJ].Residues[mj].position)+ EPS) - c_lj);
							double t(LJAccumulator + y);
							c_lj = (t-LJAccumulator)-y;
							LJAccumulator = t;
						#else
							LJAccumulator += crowderPairPotential(scalar_distance(molecules[mI].Residues[mi].position,molecules[mJ].Residues[mj].position)+ EPS);
						#endif
					}
				}
			}
			else
			#endif

			for (size_t mi=0; mi<molecules[mI].residueCount; mi++)
			{
				for (size_t mj=0; mj<molecules[mJ].residueCount; mj++)
				{

					double r(scalar_distance(molecules[mI].Residues[mi].position,molecules[mJ].Residues[mj].position) + EPS);
					#if COMPENSATE_KERNEL_SUM
						double y = (iRes.electrostaticCharge * jRes.electrostaticCharge * expf(-r/Xi) / r) - c_dh;
						double t = DHAccumulator + y;
						c_dh = (t-DHAccumulator)-y;
						DHAccumulator = t;
					#else
						DHAccumulator += (iRes.electrostaticCharge * jRes.electrostaticCharge * expf(-r/Xi) / r);
					#endif

					// first do the lennard jones pair interactions, put first to overlap the fetch woth the calcualtion of r
					float Eij(lambda*(aminoAcids.LJpotentials[iRes.aminoAcidIndex][jRes.aminoAcidIndex] - e0));


					// sigmaij is the average atomic radius determined by the van der waal radius in kim2008
					float sigmaij(0.5f * ( iRes.vanderWaalRadius + jRes.vanderWaalRadius ));
					//float r0 = powf(2.0f,(1.0f/6.0f));
					//float sigT = sigmaij / r ;
					double LJtmp(powf( sigmaij / r,6.0f)); //sigT*sigT*sigT*sigT*sigT*sigT;

					double LJ(-4.0f*Eij*LJtmp*(LJtmp-1.0f));
					//LJAccumulator += -4.0f*Eij*LJtmp*(LJtmp-1.0f);
					if (Eij>0.0f && r<sigmaij*1.122462048309372981433533049679f)  // attractive pairs
					{
						LJ = -LJ + 2.0f*Eij;
					}

					#if COMPENSATE_KERNEL_SUM
					y = LJ - c_lj;
					t = LJAccumulator + y;
					c_lj = (t-LJAccumulator)-y;
					LJAccumulator = t;
					#else
					LJAccumulator += LJ;
					#endif

				}
			}
		}
	}
	epotential = ( LJAccumulator * LJ_CONVERSION_FACTOR + DHAccumulator * DH_constant_component) * KBTConversionFactor;

#if INCLUDE_TIMERS
	CUT_SAFE_CALL( cutStopTimer(replicaEHostTimer) );
#endif
	return epotential;
}
/*
inline float Replica::fullPotentialLoop(const uint y, const uint x, const uint cacheTile, double * LJacc, double * DHacc )
{
	for (uint xi=x;xi<x+cacheTile,xi < residueCount; xi++)
		if (xi > y && contiguousResidues[y].moleculeID != contiguousResidues[xi].moleculeID)
		{
			#if REPULSIVE_CROWDING
			if(contiguousResidues[xi].aminoAcidIndex = CROWDER_IDENTIFIER)
			{
				LJacc[0] += crowderPairPotential(scalar_distance(contiguousResidues[y].position,contiguousResidues[xi].position) + EPS);
			}
			else
			{
			#endif

				float r(scalar_distance(contiguousResidues[y].position,contiguousResidues[xi].position) + EPS);
				DHacc[0] += (contiguousResidues[y].electrostaticCharge * contiguousResidues[xi].electrostaticCharge * expf(-r/Xi) / r);

				// first do the lennard jones pair interactions
				float Eij(lambda*(aminoAcids.LJpotentials[contiguousResidues[y].aminoAcidIndex][contiguousResidues[xi].aminoAcidIndex] - e0));
				// sigmaij is the average atomic radius determined by the van der waal radius in kim2008
				float sigmaij(0.5f * ( contiguousResidues[y].vanderWaalRadius + contiguousResidues[xi].vanderWaalRadius ));
				//float r0 = powf(2.0f,(1.0f/6.0f));
				//float sigT = sigmaij / r ;
				float LJtmp(powf( sigmaij / r,6.0f)); //sigT*sigT*sigT*sigT*sigT*sigT;

				float LJ(-4.0f*Eij*LJtmp*(LJtmp-1.0f));
				//LJAccumulator += -4.0f*Eij*LJtmp*(LJtmp-1.0f);
				if (Eij>0.0f && r<sigmaij*1.122462048309372981433533049679f)  // attractive pairs
				{
					LJ = -LJ + 2.0f*Eij;
				}
				LJacc[0] += LJ;

			#if REPULSIVE_CROWDING
			}
			#endif
		}
}

inline float Replica::crowderLoop(const uint y, const uint x, const uint cacheTile, double * LJacc)
{
	for (uint xi=x;xi<x+cacheTile,xi < residueCount;xi++)
		if (contiguousResidues[y].moleculeID != contiguousResidues[xi].moleculeID)
			(*LJacc) += crowderPairPotential(scalar_distance(contiguousResidues[y].position,contiguousResidues[xi].position) + EPS);
}*/

float Replica::E_sse()
{
	#if INCLUDE_TIMERS
		replicaEHostCounter++;
		CUT_SAFE_CALL( cutStartTimer(replicaEHostTimer) );
	#endif

	float epotential = 0.0f;

	double LJAccumulator = 0.0f;
	double DHAccumulator = 0.0f;
	float DH_constant_component =  DH_CONVERSION_FACTOR * 1.602176487f * 1.602176487f ;

	//float r;

	#define iRes molecules[mI].Residues[mi]
	#define jRes molecules[mJ].Residues[mj]
	// between each residue do the LJ and DH potentials
	// do summation for LJ and DH potentials
	for (size_t mI=0; mI<moleculeCount;mI++)
	{
		#if REPULSIVE_CROWDING
		if (molecules[mI].moleculeRoleIdentifier == CROWDER_IDENTIFIER)
		{
			for (size_t mJ=mI+1; mJ<moleculeCount;mJ++)
			{
				for (size_t mi=0; mi<molecules[mI].residueCount; mi++)
				{
					for (size_t mj=0; mj<molecules[mJ].residueCount; mj++)
					{
						//r = EPS + sqrtf( (iRes.position.x-jRes.position.x)*(iRes.position.x-jRes.position.x)+(iRes.position.y-jRes.position.y)*(iRes.position.y-jRes.position.y)+(iRes.position.z-jRes.position.z)*(iRes.position.z-jRes.position.z));
						//LJAccumulator += (46656.0 / pow(r,6));
						LJAccumulator += crowderPairPotential(scalar_distance(molecules[mI].Residues[mi].position,molecules[mJ].Residues[mj].position)+ EPS);
					}
				}
			}
		}
		else
		#endif
		for (size_t mJ=mI+1; mJ<moleculeCount;mJ++)
		{

			#if REPULSIVE_CROWDING
			if (molecules[mJ].moleculeRoleIdentifier == CROWDER_IDENTIFIER)
			{
				for (size_t mi=0; mi<molecules[mI].residueCount; mi++)
				{
					for (size_t mj=0; mj<molecules[mJ].residueCount; mj++)
					{
						//LJAccumulator += (46656.0 / pow(r,6));
						LJAccumulator += crowderPairPotential(scalar_distance(molecules[mI].Residues[mi].position,molecules[mJ].Residues[mj].position)+ EPS);
					}
				}
			}
			else
			#endif

			for (size_t mi=0; mi<molecules[mI].residueCount; mi++)
			{
				for (size_t mj=0; mj<molecules[mJ].residueCount; mj++)
				{

					float r(scalar_distance(molecules[mI].Residues[mi].position,molecules[mJ].Residues[mj].position) + EPS);
					DHAccumulator += (iRes.electrostaticCharge * jRes.electrostaticCharge * expf(-r/Xi) / r);

					// first do the lennard jones pair interactions, put first to overlap the fetch woth the calcualtion of r
					float Eij(lambda*(aminoAcids.LJpotentials[iRes.aminoAcidIndex][jRes.aminoAcidIndex] - e0));

					// sigmaij is the average atomic radius determined by the van der waal radius in kim2008
					float sigmaij(0.5f * ( iRes.vanderWaalRadius + jRes.vanderWaalRadius ));
					//float r0 = powf(2.0f,(1.0f/6.0f));
					//float sigT = sigmaij / r ;
					float LJtmp(powf( sigmaij / r,6.0f)); //sigT*sigT*sigT*sigT*sigT*sigT;

					float LJ(-4.0f*Eij*LJtmp*(LJtmp-1.0f));
					//LJAccumulator += -4.0f*Eij*LJtmp*(LJtmp-1.0f);
					if (Eij>0.0f && r<sigmaij*1.122462048309372981433533049679f)  // attractive pairs
					{
						LJ = -LJ + 2.0f*Eij;
					}
					LJAccumulator += LJ;
				}
			}
		}
	}
	epotential = ( LJAccumulator * LJ_CONVERSION_FACTOR + DHAccumulator * DH_constant_component) * KBTConversionFactor;

#if INCLUDE_TIMERS
	CUT_SAFE_CALL( cutStopTimer(replicaEHostTimer) );
#endif
	return epotential;
}

float Replica::E(Molecule *a,Molecule *b)
{
	float epotential = 0.0f;
	float LJAccumulator = 0.0f;
	float DHAccumulator = 0.0f;
	float DH_constant_component =  1.602176487f * 1.602176487f * DH_CONVERSION_FACTOR;
	float r;

	#define aRes a->Residues[mi]
	#define bRes b->Residues[mj]

	for (size_t mi=0; mi<a->residueCount; mi++)
	{
		for (size_t mj=0; mj<b->residueCount; mj++)
		{
			r = EPS + sqrt( (aRes.position.x-bRes.position.x)*(aRes.position.x-bRes.position.x)+(aRes.position.y-bRes.position.y)*(aRes.position.y-bRes.position.y)+(aRes.position.z-bRes.position.z)*(aRes.position.z-bRes.position.z));
			DHAccumulator += aRes.electrostaticCharge * bRes.electrostaticCharge * expf(-r/Xi) / r;

			// first do the lennard jones pair interactions
			float Eij = lambda*(aminoAcids.LJpotentials[aRes.aminoAcidIndex][bRes.aminoAcidIndex] - e0);
			// sigmaij is the average atomic radius determined by the van der waals radius in kim2008
			float sigmaij = 0.5f * ( aRes.vanderWaalRadius + bRes.vanderWaalRadius );
			float LJtmp = powf( sigmaij / r,6.0f); //sigT*sigT*sigT*sigT*sigT*sigT;
			float LJ = -4.0f*Eij*LJtmp*(LJtmp-1.0f);
			if (Eij>0.0f && r<sigmaij*1.122462048309372981433533049679f)  // attractive pairs
			{
				LJ = -LJ + 2.0f*Eij;
			}
			LJAccumulator += LJ;
		}
	}
	epotential = ( DHAccumulator * DH_constant_component + LJAccumulator * LJ_CONVERSION_FACTOR ) * KBTConversionFactor;
	return epotential;
}

float Replica::geometricDistance(Molecule *a,Molecule *b, float cutoff)
{
	float minimum = (a->center - b->center).magnitude();
	float tmp;
	for (size_t mi=0; mi<a->residueCount; mi++)
	{
		for (size_t mj=0; mj<b->residueCount; mj++)
		{
			tmp = Vector3f(a->Residues[mi].position - b->Residues[mj].position).magnitude();
			if (tmp < minimum)
				minimum = tmp;
		}
	}
	return minimum;
}

/*
float length(const float4 a, const float4 b)
{
	return sqrtf((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z));
}

// evaluate the conformations  energy
float Replica::Eopt()
{
#if INCLUDE_TIMERS
	replicaEHostCounter++;
	CUT_SAFE_CALL( cutStartTimer(replicaEHostTimer) );
#endif

	float epotential(0.0f);

	float DHAccumulator= 0.0f;
	float LJAccumulator= 0.0f;

	for (size_t tx=0; tx<residueCount;tx++)
	{
		for (size_t ty=tx+1; ty<residueCount;ty++)
		{
			float4 rIp = host_float4_residuePositions[tx];
			float4 rJp = host_float4_residuePositions[ty];

			if (rIp.w == rJp.w)		// same molecule
			{
				// residueIp.w == residueJp.w means that they are the same molecule
				// if either is -1 then its a padding residue and must be ignored
			}
			else
			{

				float4 rIm = host_float4_residueMeta[tx];
				float4 rJm = host_float4_residueMeta[ty];

				float r = length(rIp,rJp) + EPS;  // add eps so that r is never 0, can happen in a collision


				//if there are crowders the LJ is replaced by the repulsive formula.
				float LJ(0);   // needs to be here for scoping issues
				float DH(0);

				#if REPULSIVE_CROWDING
				if (rIm.w == CROWDER_IDENTIFIER || rJm.w == CROWDER_IDENTIFIER) // repulsive crowder interaction
				{
					LJ = (6.0f*6.0f*6.0f*6.0f*6.0f*6.0f)/(r*r*r*r*r*r);  // powf(6.0f/r,6);
				}
				else  // normal LJ interaction
				{
				#endif
					// do the debye huckel long range intereactions
					DH = rIm.y * rJm.y * 1.602176487f * 1.602176487f * expf(-r/Xi) / r;
					float sigmaij = (rIm.z + rJm.z ) * 0.5f;
					float Eij = lambda*(aminoAcids.LJpotentials[int(rIm.x)][int(rJm.x)] - e0);


					// sigmaij is the average atomic radius determined by the van der waals radius in kim2008
					//float r0 = sigmaij*1.122462048309372981433533049679f; //pow(2.0,(1.0/6.0));
					// the following 2 lines are needed for preserving numerical accuracy on a cpu. not needed on GPU
					//float LJtmp = powf(sigmaij/r, 6.0f);
					float LJtmp = (sigmaij*sigmaij*sigmaij*sigmaij*sigmaij*sigmaij)/(r*r*r*r*r*r);

					LJ = -4.0f*Eij*LJtmp*(LJtmp-1.0f);
					if (Eij>0.0f)  // attractive pairs
					{
						//if (r<r0)
						if (r<sigmaij*1.122462048309372981433533049679f)
						{
							LJ = -LJ + 2.0f*Eij;
						}
					}
				#if REPULSIVE_CROWDING
				}  // end conditional branch for LJ or repulsive short-range energy
				#endif
				//epotential += LJ * LJ_CONVERSION_FACTOR + DH; // convert to kcal/mol
				DHAccumulator += DH;
				LJAccumulator += LJ;
			}
		}
	}
	//epotential *= KBTConversionFactor;
	epotential = (DHAccumulator * DH_CONVERSION_FACTOR + LJAccumulator * LJ_CONVERSION_FACTOR) * KBTConversionFactor;
	//cout << "opt:  " << epotential << endl;
	//cout << DHAccumulator * DH_CONVERSION_FACTOR * KBTConversionFactor << " " <<  LJAccumulator * LJ_CONVERSION_FACTOR * KBTConversionFactor << endl;

	/*
	float epotential = 0.0f;

	float LJtot = 0.0f;
	float DHtot = 0.0f;

	//float LJAccumulator = 0.0f;
	//float DHAccumulator = 0.0f;
	float DH_constant_component =  1.602176487f * 1.602176487f * DH_CONVERSION_FACTOR;

	float r;

	#define	 CPU_STRIDE 16

	int gridDim = contiguousResiduesSize/CPU_STRIDE;

	//Residue ResArrX[CPU_STRIDE];
	//Residue ResArrY[CPU_STRIDE];
	//memcpy (&ResArrX,contiguousResidues+gx*CPU_STRIDE,CPU_STRIDE);
	/memcpy (&ResArrY,contiguousResidues+gy*CPU_STRIDE,CPU_STRIDE);

	for (size_t gx=0; gx<gridDim;gx++)
	{
		for (size_t gy=gx; gy<gridDim;gy++)
		{

			float LJAccumulator = 0.0f;
			float DHAccumulator = 0.0f;

			for (int x=0;x<CPU_STRIDE;x++)   // go through the first set
			{
				int gxx = gx*CPU_STRIDE+x;

				bool xcrowd = ((contiguousResidues+gxx)->aminoAcidIndex == CROWDER_IDENTIFIER);

				if ((contiguousResidues+gxx)->aminoAcidIndex != PADDER_IDENTIFIER)  // its a padder -> ignore
				{
					for (int y=0; y<CPU_STRIDE;y++)
					{
						int gyy = gy*CPU_STRIDE+y;

						if ((contiguousResidues+gyy)->aminoAcidIndex != PADDER_IDENTIFIER)  // not a padder or crowder
						{
							if ((contiguousResidues+gxx)->moleculeID!=(contiguousResidues+gyy)->moleculeID) // if they are different molecules
							{
								#if REPULSIVE_CROWDING
								if (xcrowd || (contiguousResidues+gyy)->aminoAcidIndex == CROWDER_IDENTIFIER) // one is a crowder
								{
									r = EPS + sqrt( ((contiguousResidues+gxx)->position.x-(contiguousResidues+gyy)->position.x)*((contiguousResidues+gxx)->position.x-(contiguousResidues+gyy)->position.x)+((contiguousResidues+gxx)->position.y-(contiguousResidues+gyy)->position.y)*((contiguousResidues+gxx)->position.y-(contiguousResidues+gyy)->position.y)+((contiguousResidues+gxx)->position.z-(contiguousResidues+gyy)->position.z)*((contiguousResidues+gxx)->position.z-(contiguousResidues+gyy)->position.z));
									LJAccumulator += powf(6.0f/r,6.0f);
									// only use repulsive force
									//DHAccumulator += iRes.electrostaticCharge * jRes.electrostaticCharge * expf(-r/Xi) / r;
								}
								else
								#endif //REPULSIVE_CROWDING
								{
									r = EPS + sqrt( ((contiguousResidues+gxx)->position.x-(contiguousResidues+gyy)->position.x)*((contiguousResidues+gxx)->position.x-(contiguousResidues+gyy)->position.x)+((contiguousResidues+gxx)->position.y-(contiguousResidues+gyy)->position.y)*((contiguousResidues+gxx)->position.y-(contiguousResidues+gyy)->position.y)+((contiguousResidues+gxx)->position.z-(contiguousResidues+gyy)->position.z)*((contiguousResidues+gxx)->position.z-(contiguousResidues+gyy)->position.z));
									DHAccumulator += (contiguousResidues+gxx)->electrostaticCharge * (contiguousResidues+gyy)->electrostaticCharge * expf(-r/Xi) / r;

									// first do the lennard jones pair interactions
									float Eij = lambda*(aminoAcids.LJpotentials[(contiguousResidues+gxx)->aminoAcidIndex][(contiguousResidues+gyy)->aminoAcidIndex] - e0);
									// sigmaij is the average atomic radius determined by the van der waal radius in kim2008
									float sigmaij = 0.5f * ( (contiguousResidues+gxx)->vanderWaalRadius + (contiguousResidues+gyy)->vanderWaalRadius );
									//float r0 = powf(2.0f,(1.0f/6.0f));
									//float sigT = sigmaij / r ;
									float LJtmp = powf( sigmaij / r,6.0f); //sigT*sigT*sigT*sigT*sigT*sigT;

									float LJ = -4.0f*Eij*LJtmp*(LJtmp-1.0f);
									//LJAccumulator += -4.0f*Eij*LJtmp*(LJtmp-1.0f);
									if (Eij>0.0f && r<sigmaij*1.122462048309372981433533049679f)  // attractive pairs
									{
										LJ = -LJ + 2.0f*Eij;
									}
									LJAccumulator += LJ;
								}
							}
						}
					}
				}

			}
			if (gx==gy)
			{
				LJtot += 0.5f*LJAccumulator;
				DHtot += 0.5f*DHAccumulator;
			}
			else
			{
				LJtot += LJAccumulator;
				DHtot += DHAccumulator;
			}
		}
	}
	epotential = ( DHtot * DH_constant_component + LJtot * LJ_CONVERSION_FACTOR ) * KBTConversionFactor;

	cout << DHtot * DH_constant_component * KBTConversionFactor << " " <<  LJtot * LJ_CONVERSION_FACTOR * KBTConversionFactor << endl;*/
/*#if INCLUDE_TIMERS
	CUT_SAFE_CALL( cutStopTimer(replicaEHostTimer) );
#endif
	return epotential;
}
*/
// evaluate conformations energy with optimizations
float Replica::Eorig()
{
#if INCLUDE_TIMERS
	replicaEHostCounter++;
	CUT_SAFE_CALL( cutStartTimer(replicaEHostTimer) );
#endif

	float epotential = 0.0f;

	// between each residue do the LJ and DH potentials
	// do summation for LJ and DH potentials
	for (size_t mI=0; mI<moleculeCount;mI++)
	{
		for (size_t mJ=mI+1; mJ<moleculeCount;mJ++)
		{
			for (size_t mi=0; mi<molecules[mI].residueCount; mi++)
			{
				for (size_t mj=0; mj<molecules[mJ].residueCount; mj++)
				{
					if (mI!=mJ) //mj == mi && ) // if this is the same molecule, the same residue cannot be compared to itself, do nothing
					{
						epotential += phi(molecules[mI].Residues[mi],molecules[mJ].Residues[mj]);
                    }
				}
			}
		}
	}
#if INCLUDE_TIMERS
	CUT_SAFE_CALL( cutStopTimer(replicaEHostTimer) );
#endif
	return epotential * KBTConversionFactor;
}

void Replica::printTimers()
{
	#if INCLUDE_TIMERS

	float ave_cputimer=cutGetTimerValue(replicaEHostTimer)/(replicaEHostCounter==0 ? 1: float(replicaEHostCounter));
	float ave_gputimer=cutGetTimerValue(replicaECUDATimer)/(replicaECUDACounter==0 ? 1: float(replicaECUDACounter));

	cout << "Timer count and values (ms)" << endl;
	cout << "Count\tTotal(ms)   \tAverage(ms)  \t Action" << endl;

	printf("%5d\t%12.6f\t%12.6f\t Replica to GPU (init & transfer)\n"	,replicaToGPUCounter,			cutGetTimerValue(replicaToGPUTimer),			cutGetTimerValue(replicaToGPUTimer)/(replicaToGPUCounter==0 ? 1: float(replicaToGPUCounter)));
	printf("%5d\t%12.6f\t%12.6f\t Replica Update on GPU (transter)\n"	,replicaUpdateGPUCounter,		cutGetTimerValue(replicaUpdateGPUTimer),		cutGetTimerValue(replicaUpdateGPUTimer)/(replicaUpdateGPUCounter==0 ? 1: float(replicaUpdateGPUCounter)));
	printf("%5d\t%12.6f\t%12.6f\t Kernel Timer (computation)\n"			,replicaECUDACounter,			cutGetTimerValue(replicaECUDATimer),			cutGetTimerValue(replicaECUDATimer)/(replicaECUDACounter==0 ? 1: float(replicaECUDACounter)));
	printf("%5d\t%12.6f\t%12.6f\t Host Timer (computation)\n"			,replicaEHostCounter,			cutGetTimerValue(replicaEHostTimer),			cutGetTimerValue(replicaEHostTimer)/(replicaEHostCounter==0 ? 1: float(replicaEHostCounter)));
	printf("%5d\t%12.6f\t%12.6f\t Update Molecule on GPU (transfer)\n"	,replicaMoleculeUpdateCounter,	cutGetTimerValue(replicaMoleculeUpdateTimer),	cutGetTimerValue(replicaMoleculeUpdateTimer)/(replicaMoleculeUpdateCounter==0 ? 1: float(replicaMoleculeUpdateCounter)));
	printf("%5d\t%12.6f\t%12.6f\t Mutation of GPU (computation) \n"		,replicaDeviceMCCounter,		cutGetTimerValue(replicaDeviceMCTimer),			cutGetTimerValue(replicaDeviceMCTimer)/(replicaDeviceMCCounter==0 ? 1: float(replicaDeviceMCCounter)));
	printf("%5d\t%12.6f\t%12.6f\t GPU Memory Initialisation (malloc)\n"	,initGPUMemoryCounter,			cutGetTimerValue(initGPUMemoryTimer),			cutGetTimerValue(initGPUMemoryTimer)/(initGPUMemoryCounter==0 ? 1: float(initGPUMemoryCounter)));
	printf("Kernel Speedup:  %0.1fx\n", ave_cputimer/ave_gputimer);

	#else

	cout << "Timers disabled: set INCLUDE_TIMERS 1 and USING_CUDA 1 in definitions.h" << endl;

	#endif
}

void Replica::setBoundingValue(float value)
{
	boundingValue = value;
}

#if USING_CUDA
// all functions following this line are dependent on CUDA

void Replica::initialiseDeviceLJPotentials()
{
	copyLJPotentialDataToDevice (device_LJPotentials, &aminoAcids);
}
/*
#if LJ_LOOKUP_METHOD == TEXTURE_MEM
void bindLJtoTexture(int device, float *ljp)
{
	bindLJTexture(device, ljp, sizeof(float)*AA_COUNT*AA_COUNT); // typically 4*20*20
}

textureReference* bindLJTexture2D(float *ljp)
{
	cudaGetTextureReference(&texRefPtr, "LJTexture");
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
	cudaBindTexture2D(0, texRefPtr, ljp, &channelDesc, AA_COUNT, AA_COUNT, AA_COUNT*sizeof(float));
	return texRefPtr;
}

void unbindLJTexture2D (textureReference* texRefPtr)
{
	cudaUnbindTexture(LJTexture2D);
	delete texRefPtr;
}

void FreeLJfromTexture(int device)
{
	freeLJTexture(device);
}
#endif*/

int Replica::getBlockSize()
{
    return blockSize;
}

void Replica::setBlockSize(int blockSize)
{
    this->blockSize = blockSize;

    // also set the shared memory required per block as the 2 are linked
    sharedMemSize = sizeof(float);
	#if POSITIONDATA_MEMORY == SHARED_MEM
		sharedMemSize += sizeof(float4);
	#endif
	#if METADATA_MEMORY == SHARED_MEM
		sharedMemSize += sizeof(float4);
	#endif
	sharedMemSize *= blockSize;
}

int Replica::getDataSetSize()
{
    return dataSetSize;
}

void Replica::setDataSetSize(int dataSetSize)
{
    this->dataSetSize = dataSetSize;
}

int Replica::getPaddedSize()
{
    return paddedSize;
}

void Replica::setPaddedSize(int paddedSize)
{
    this->paddedSize = paddedSize;
}

int Replica::getGridSize()
{
    return gridSize;
}

void Replica::setGridSize(int gridSize)
{
    this->gridSize = gridSize;
}

int Replica::getResultSize()
{
    return resultSize;
}

void Replica::setResultSize(int resultSize)
{
    this->resultSize = resultSize;
}


#if CUDA_STREAMS
void Replica::ReserveSumSpace()
{
	if (!replicaIsOnDevice)
	{
		cout << "! Error: replica not on device implies paddedSize==0; cannot perform Replica::ReserveSumSpace()" << endl;
	}
	// result stored on the device
	// gridSize can be arbitrary
	gridSize = paddedSize/blockSize;
	resultSize = gridSize;
	//for a parallel sum each grid must have one cell in the array of results from all the threads
	cudaMallocHost((void**)&kernelResult, sizeof(float)*resultSize*resultSize);
	cudaMalloc((void **)&device_kernelResult,sizeof(float)*resultSize*resultSize);

	cudaError_t err;
	if ((err = cudaGetLastError()) != cudaSuccess)
	  printf("CUDA error: %s\n", cudaGetErrorString(err));

	memset(kernelResult,0,sizeof(float)*resultSize*resultSize);
	//CUDA_memcpy_to_device_async(device_kernelResult,kernelResult,sizeof(float)*resultSize*resultSize,cudaStream);
	//CUDA_memcpy_to_device(device_kernelResult,kernelResult,sizeof(float)*resultSize*resultSize);
	cudaMemset(device_kernelResult, 0, sizeof(float)*resultSize*resultSize);
}

void Replica::FreeSumSpace()
{
	cudaFreeHost(kernelResult);
	cudaFree((void**)&device_kernelResult);
}

float Replica::SumGridResults()
{
	cudaStreamSynchronize(cudaStream);
	float potentialSum = 0.0f;
	float c(0.0f); //kahanSum
	for (int i=0;i<resultSize*resultSize;i++)
	{
			float y(kernelResult[i] - c);
			float t((potentialSum) + y);
			c = (t-(potentialSum)) -y;
			potentialSum = t;
	}


	//for (int i=0;i<resultSize*resultSize;i++)
	//{
	//	potentialSum += kernelResult[i];
	//}

	cudaMemset(device_kernelResult, 0, sizeof(float)*resultSize*resultSize);
	return potentialSum;
}
#endif

#if CUDA_MC
void Replica::MCSearchOnDevice(int steps)
{
	float translateStep = INITIAL_TRANSLATIONAL_STEP;
	float rotateStep = INITIAL_ROTATIONAL_STEP;

	Vector3f translateVector;
	Vector3f rotateVector;

	float oldPotential = potential;

	for (int step=0;step<steps;step++)
	{
		//unsigned int type(_gsl_ran_bernoulli (MCRng,translate_rotate_bernoulli_bias));
		uint moleculeNo = ((int)(1000.0*_gsl_qrng_get (qrng_moleculeSelection)))%((int)moleculeCount);

		switch (_gsl_ran_bernoulli (MCRng,translate_rotate_bernoulli_bias))
		{
			case _rotate:
				rotateVector = createNormalisedRandomVector(rng_rotate);
				rotateOnDevice(moleculeNo,rotateVector,rotateStep);
				break;

			case _translate:
				// pick directions to rotate in, do by generating a random normalised vector of length equal to INITIAL_TRANSLATIONAL_STEP
				translateVector = translateStep * createNormalisedRandomVector(rng_translate);

				// if the translate falls inside the bounding radius
				if ((molecules[moleculeNo].center + translateVector).sumSquares() < BOUNDING_RADIUS*BOUNDING_RADIUS)
				{
					translateOnDevice(moleculeNo,translateVector);
				}
				// if the change causes the molecule to fall outside the bounding radius
				else
				{
#ifdef VERY_VERBOSE
					cout << "  - Reject: Replica "<< label << "/Molecule " << moleculeNo <<  " : falls outside bounding sphere" << endl;
#endif
				}
				break;

			default:
				break;
		}

		EonDeviceAsync();
		float newPotential = SumGridResults();
		oldPotential = potential;  // save the current potential value
#ifdef 	CUDA_STREAMS
		// we have to sync here to be sure that the correct value is in newPotential, the rest is async
		cudaStreamSynchronize(cudaStream);
#endif

		float delta = newPotential - oldPotential;
		// if the change is bad then discard it.
		if (delta < 0)
		{
			potential = newPotential;
#ifdef VERY_VERBOSE
			cout << "  * Replace: Replica "<< label << "/Molecule " << moleculeNo <<  " : delta E = " << delta << endl;
#endif
		}
		else if (gsl_rng_uniform(MCKbRng)>exp(-delta*4.184f/(Rgas*temperature)))
		{
			potential = newPotential;
#ifdef VERY_VERBOSE
			cout << "  **Replace: Replica "<< label << "/Molecule " << moleculeNo <<  " : delta E = " << delta << " U < " << exp(-delta*4.184f/(Rgas*temperature)) << endl;
#endif
		}
		else
		{
#ifdef VERY_VERBOSE
			cout << "  - Reject: Replica "<< label << "/Molecule " << moleculeNo <<  " : delta E = " << delta << endl;
#endif
			cudaRollbackMutation();  // async rollback so that it queues up the next mutation.
		}
	}
	cout << " == End MC search. Replica:" << label << endl;
}
#endif

#if CUDA_STREAMS

// 1/3 of the above function, does the mutation on the gpu asynchronously
void Replica::MCSearchMutate()
{
	float translateStep = INITIAL_TRANSLATIONAL_STEP;
	double rotateStep = INITIAL_ROTATIONAL_STEP;

	oldPotential = potential;

	uint moleculeNo = (int) gsl_rng_uniform_int(rng_moleculeSelection,moleculeCount);
	lastMutationIndex = moleculeNo;

	uint mutationType = gsl_ran_bernoulli (MCRng,translate_rotate_bernoulli_bias);

	// save the current state so we can role back if it was not a good mutation.
	savedMolecule.saveBeforeStateChange(&molecules[moleculeNo]);

	switch (mutationType)
	{
		case _rotate:
		{
			Vector3double rotateVector = createNormalisedRandomVectord(rng_rotate);
			molecules[moleculeNo].rotateQ(rotateVector,rotateStep);

			#if OUTPUT_LEVEL >= PRINT_MC_MUTATIONS
			cout << "    Rotate: Replica "<< label << "/Molecule " << moleculeNo << endl;
			#endif
		}
		break;

		case _translate:
		{
			// pick directions to rotate in, do by generating a random normalised vector of length equal to INITIAL_TRANSLATIONAL_STEP
			Vector3f translateVector = translateStep * createNormalisedRandomVector(rng_translate);

			#if OUTPUT_LEVEL >= PRINT_MC_MUTATIONS
				cout << "    Translate: Replica "<< label << "/Molecule " << moleculeNo << endl;
			#endif

			// bounding sphere conditions
			#if BOUNDING_METHOD == BOUNDING_SPHERE
				// if the translate falls inside the bounding radius
				if ((molecules[moleculeNo].center + translateVector).sumSquares() < BOUNDING_RADIUS*BOUNDING_RADIUS)
				{
					molecules[moleculeNo].translate(translateVector);
				}
				// if the change causes the molecule to fall outside the bounding radius
				else
				{
					#if OUTPUT_LEVEL >= PRINT_MC_MUTATIONS
						cout << "  - Reject: Replica "<< label << "/Molecule " << moleculeNo <<  " : falls outside bounding sphere" << endl;
					#endif
				}
			#elif BOUNDING_METHOD == PERIODIC_BOUNDRY 	// periodic boundary conditions
				Vector3f newPosition = molecules[moleculeNo].center + translateVector;
				newPosition.x = fmod(newPosition.x+BOUNDING_VALUE,BOUNDING_VALUE);
				newPosition.y = fmod(newPosition.y+BOUNDING_VALUE,BOUNDING_VALUE);
				newPosition.z = fmod(newPosition.z+BOUNDING_VALUE,BOUNDING_VALUE);
				molecules[moleculeNo].setPosition(newPosition);
			#endif
		}
		break;
		default:
				break;
	}
	MoleculeDataToDevice(moleculeNo);
	newPotential = 0.0f;
}

void Replica::MCSearchEvaluate()
{
	EonDeviceAsync();
}

void Replica::MCSearchAcceptReject()
{
	cudaStreamSynchronize(cudaStream);  // sync, newPotential needs to have been returned
	newPotential = SumGridResults();

	float delta = (newPotential - oldPotential);  // needs to be in K_bT
	// if the change is bad then discard it.
	if (delta < 0)
	{
		potential = newPotential;
		oldPotential = potential;
		#if OUTPUT_LEVEL >= PRINT_MC_MUTATIONS
		cout << "  * Replace: Replica "<< label << "/Molecule " << lastMutationIndex <<  " : delta E = " << delta << " E = " << potential << endl;
		#endif
	}
	else if (gsl_rng_uniform(MCKbRng)<exp(-delta*4.184f/(Rgas*temperature)))
	{
		potential = newPotential;
		oldPotential = potential;
		#if OUTPUT_LEVEL >= PRINT_MC_MUTATIONS
		cout << "  **Replace: Replica "<< label << "/Molecule " << lastMutationIndex <<  " : delta E = " << delta << " U < " << exp(-delta*4.184f/(Rgas*temperature)) << " E = " << potential << endl;
		#endif
	}
	else
	{
		#if OUTPUT_LEVEL >= PRINT_MC_MUTATIONS
			cout << "  - Reject: Replica "<< label << "/Molecule " << lastMutationIndex <<  " : delta E = " << delta << " E = " << oldPotential << endl;
		#endif
		molecules[lastMutationIndex].undoStateChange(&savedMolecule);
		potential = oldPotential;
		MoleculeDataToDevice(lastMutationIndex); // you have to update the device again because the copy will be inconsistent
	}
}
#endif  // streams

// copy the replica to the device
void Replica::ReplicaDataToDevice()
{
#if INCLUDE_TIMERS
	replicaToGPUCounter++;
	initGPUMemoryCounter++;
	CUT_SAFE_CALL(cutStartTimer(initGPUMemoryTimer));
#endif

	// allocate total size needed.

	paddedSize = int(ceil(float(residueCount)/float(blockSize)))*blockSize; // reserve a blocksize multiple because it allows for efficient summation

	dataSetSize = paddedSize;





	#if CUDA_STREAMS
	cudaMallocHost((void**)&host_float4_residuePositions,sizeof(float4)*paddedSize);
	#else
	host_float4_residuePositions = new float4[paddedSize];
	#endif

	//host_float4_residuePositions = new float4[paddedSize];
	host_float4_residueMeta = new float4[paddedSize];
	host_moleculeStartPositions = new int[moleculeCount];

	cudaMalloc((void**)&device_float4_residuePositions,sizeof(float4)*paddedSize);
	cudaMalloc((void**)&device_float4_residueMeta,sizeof(float4)*paddedSize);
	cudaMalloc((void**)&device_moleculeStartPositions,sizeof(int)*moleculeCount);
	cudaMalloc((void**)&device_residueCount,sizeof(int));
	cudaMalloc((void**)&device_moleculeCount,sizeof(int));

#if INCLUDE_TIMERS
	CUT_SAFE_CALL(cutStopTimer(initGPUMemoryTimer));
	CUT_SAFE_CALL(cutStartTimer(replicaToGPUTimer));
#endif
	// keep a counter to calculate the offset from the beginning of the residues array
	int arrayIndex = 0;
	// pack the residues into an array and copy it to the correct place on the GPU
	for (int m=0; m<moleculeCount;m++)
	{
		host_moleculeStartPositions[m] = arrayIndex;   // start position of the first molecule in the host/gpu array

		float mf(m);

		for (int rc=0; rc<molecules[m].residueCount;rc++)
		{
			host_float4_residuePositions[arrayIndex].x = molecules[m].Residues[rc].position.x;
			host_float4_residuePositions[arrayIndex].y = molecules[m].Residues[rc].position.y;
			host_float4_residuePositions[arrayIndex].z = molecules[m].Residues[rc].position.z;
			host_float4_residuePositions[arrayIndex].w = mf;  // residue belongs to this molecule id
			host_float4_residueMeta[arrayIndex].x = molecules[m].Residues[rc].aminoAcidIndex;
			host_float4_residueMeta[arrayIndex].y = molecules[m].Residues[rc].electrostaticCharge;
			host_float4_residueMeta[arrayIndex].z = molecules[m].Residues[rc].vanderWaalRadius;
			host_float4_residueMeta[arrayIndex].w = molecules[m].moleculeRoleIdentifier;

			arrayIndex++;
		}
	}

	while ( arrayIndex < paddedSize )
	{
		host_float4_residuePositions[arrayIndex].x = 0;//float(rand()); // stops division by zero
		host_float4_residuePositions[arrayIndex].y = 0.0f;
		host_float4_residuePositions[arrayIndex].z = 0.0f;
		host_float4_residuePositions[arrayIndex].w = PADDER_IDENTIFIER;  // residue belongs to this molecule id

		//h_id[arrayIndex] = PADDER_IDENTIFIER;

		host_float4_residueMeta[arrayIndex].x = PADDER_AATYPE;  // amino acid index
		host_float4_residueMeta[arrayIndex].y = 0.0f;
		host_float4_residueMeta[arrayIndex].z = 0.0f;
		host_float4_residueMeta[arrayIndex].w = 0.0f;
		arrayIndex++;
	}


// cant stream as the amount of stuff is larger than pagable memory
	#if 	CUDA_STREAMS
	cudaMemcpyAsync(device_float4_residuePositions, host_float4_residuePositions, sizeof(float4)*paddedSize, cudaMemcpyHostToDevice,cudaStream);
	cudaStreamSynchronize(cudaStream);
	#else
	// copy all the above to the device
	cudaMemcpy(device_float4_residuePositions, host_float4_residuePositions, sizeof(float4)*paddedSize, cudaMemcpyHostToDevice);
	#endif
	cudaMemcpy(device_float4_residueMeta, host_float4_residueMeta, sizeof(float4)*paddedSize, cudaMemcpyHostToDevice);
	cudaMemcpy(device_moleculeStartPositions, host_moleculeStartPositions, sizeof(int)*moleculeCount, cudaMemcpyHostToDevice);
	cudaMemcpy(device_residueCount, &residueCount, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(device_moleculeCount, &moleculeCount, sizeof(int), cudaMemcpyHostToDevice);


	#if METADATA_MEMORY == TEXTURE_MEM
	bindMetaDataToTexture(device_float4_residueMeta,sizeof(float4)*paddedSize);
	#endif
	#if POSITIONDATA_MEMORY == TEXTURE_MEM
	bindPositionDataToTexture(device_float4_residuePositions,sizeof(float4)*paddedSize);
	#endif


#if INCLUDE_TIMERS
	CUT_SAFE_CALL( cutStopTimer(replicaToGPUTimer) );
	//cout << "Data to GPU Initialisation and Transfer time: " << (cutGetTimerValue(initGPUMemoryTimer)+cutGetTimerValue(replicaToGPUTimer)) << "ms"<< endl;
#endif

	replicaIsOnDevice = true;

	return;
}

void Replica::setLJpotentials(float *ljp)
{
	device_LJPotentials = ljp;
}

void Replica::freeLJpotentials()
{
	cudaFree(device_LJPotentials);
}

void Replica::ReplicaDataToHost()
{

#if INCLUDE_TIMERS
	replicaToHostCounter++;
	CUT_SAFE_CALL(cutStartTimer(replicaToHostTimer));
#endif
	int residueCount = 0;

	for (int m=0; m<moleculeCount;m++)
	{
		residueCount += molecules[m].residueCount;
	}


//#ifdef 	CUDA_STREAMS ==== doesn't work because host code needs it to be synchronous for now
//	CUDA_memcpy_to_host_async(device_float4_residuePositions, host_float4_residuePositions, sizeof(float4)*dataSetSize, cudaStream);
//	CUDA_memcpy_to_host_async(device_float4_residueMeta, host_float4_residueMeta, sizeof(float4)*dataSetSize, cudaStream);
//#else
	// copy all the above to the device
	cudaMemcpy(host_float4_residuePositions, device_float4_residuePositions, sizeof(float4)*paddedSize, cudaMemcpyDeviceToHost);
	//CUDA_memcpy_to_host(host_float4_residuePositions, device_float4_residuePositions, sizeof(float4)*paddedSize);

	// I dont think this next line is neccessary because it would not have changed on the device
	//CUDA_memcpy_to_host(device_float4_residueMeta, host_float4_residueMeta, sizeof(float4)*paddedSize);
//#endif


#if INCLUDE_TIMERS
	CUT_SAFE_CALL( cutStopTimer(replicaToHostTimer) );
#endif
	// position of the molecule in the gpu array
	int arrayIndex = 0;
	// unpack the residues from the host array
	for (int m=0; m<moleculeCount;m++)
	{
		arrayIndex = host_moleculeStartPositions[m];   // start position of the first molecule in the host/gpu array
		for (int offset=0; offset<molecules[m].residueCount;offset++)
		{
			molecules[m].Residues[offset].position.x = host_float4_residuePositions[arrayIndex+offset].x;
			molecules[m].Residues[offset].position.y = host_float4_residuePositions[arrayIndex+offset].y;
			molecules[m].Residues[offset].position.z = host_float4_residuePositions[arrayIndex+offset].z;
		}
	}

#if INCLUDE_TIMERS
	CUT_SAFE_CALL( cutStopTimer(replicaToHostTimer) );
#endif
	return;
}


void Replica::UpdateDeviceData()
{
	if (!replicaIsOnDevice)
	{
		cout << "ERROR: Replica::UpdateDeviceData() called without initialising device data"<< endl;
		ReplicaDataToDevice();
	}
#if INCLUDE_TIMERS
	replicaUpdateGPUCounter++;
	CUT_SAFE_CALL( cutStartTimer(replicaUpdateGPUTimer) );
#endif
	// works only if the replica is already on the device
	// keep a counter to calculate the offset from the beginning of the residues array
	int arrayIndex = 0;
	// pack the residues into an array and copy it to the correct place on the GPU
	for (int m=0; m<moleculeCount;m++)
	{
		host_moleculeStartPositions[m] = arrayIndex;   // start position of the first molecule in the host/gpu array
		for (int rc=0; rc<molecules[m].residueCount;rc++)
		{
			host_float4_residuePositions[arrayIndex].x = molecules[m].Residues[rc].position.x;
			host_float4_residuePositions[arrayIndex].y = molecules[m].Residues[rc].position.y;
			host_float4_residuePositions[arrayIndex].z = molecules[m].Residues[rc].position.z;
			host_float4_residuePositions[arrayIndex].w = float(m);  // residue belongs to this molecule id
		//	host_float4_residueMeta[arrayIndex].x = molecules[m].Residues[rc].aminoAcidIndex;
		//	host_float4_residueMeta[arrayIndex].y = molecules[m].Residues[rc].electrostaticCharge;
		//	host_float4_residueMeta[arrayIndex].z = molecules[m].Residues[rc].vanderWaalRadius;
		//	host_float4_residueMeta[arrayIndex].w = temperature;
			arrayIndex++;
		}
	}

	// copy all the above to the device
#if 	CUDA_STREAMS
	cudaMemcpyAsync(device_float4_residuePositions, host_float4_residuePositions, sizeof(float4)*paddedSize, cudaMemcpyHostToDevice,cudaStream);
	//CUDA_memcpy_to_device_async(device_float4_residuePositions, host_float4_residuePositions, sizeof(float4)*paddedSize,cudaStream);

#else
	cudaMemcpy(device_float4_residuePositions, host_float4_residuePositions, sizeof(float4)*paddedSize, cudaMemcpyHostToDevice);
	//CUDA_memcpy_to_device(device_float4_residuePositions, host_float4_residuePositions, sizeof(float4)*paddedSize);

#endif
#if INCLUDE_TIMERS
	CUT_SAFE_CALL( cutStopTimer(replicaUpdateGPUTimer) );
#endif
	return;
}

// update a molecule on the device
void Replica::MoleculeDataToDevice(int MoleculeID)
{
	if (!replicaIsOnDevice)
	{
		cout << "ERROR: Replica::MoleculeDataToDevice("<< MoleculeID << ") called without initialising device data."<< endl;
		ReplicaDataToDevice();
	}
#if INCLUDE_TIMERS
	replicaMoleculeUpdateCounter++;
	CUT_SAFE_CALL( cutStartTimer(replicaMoleculeUpdateTimer) );
#endif

	// TODO: change to use float4 for everything, will eliminate this copy

	/*
	1. get molecule on device's residues pointer
	2. make residue array
	3. copy to device at location found in 1
	*/

	int residueIndex = host_moleculeStartPositions[MoleculeID];
	int memoryPosition = residueIndex;
	int moleculeSize = molecules[MoleculeID].residueCount;
	for (int rc=0; rc<molecules[MoleculeID].residueCount;rc++)
	{
		host_float4_residuePositions[residueIndex].x = molecules[MoleculeID].Residues[rc].position.x;
		host_float4_residuePositions[residueIndex].y = molecules[MoleculeID].Residues[rc].position.y;
		host_float4_residuePositions[residueIndex].z = molecules[MoleculeID].Residues[rc].position.z;
		host_float4_residuePositions[residueIndex].w = float(MoleculeID);  // residue belongs to this molecule id
		residueIndex++;
	}

	// copy to the device
#if CUDA_STREAMS
	cudaMemcpyAsync(&device_float4_residuePositions[memoryPosition], &host_float4_residuePositions[memoryPosition], sizeof(float4)*moleculeSize,cudaMemcpyHostToDevice,cudaStream);
	//CUDA_memcpy_to_device_async(&device_float4_residuePositions[memoryPosition], &host_float4_residuePositions[memoryPosition], sizeof(float4)*moleculeSize,cudaStream);
	//residueMeta will already be on the device
#else
	cudaMemcpy(&device_float4_residuePositions[memoryPosition], &host_float4_residuePositions[memoryPosition], sizeof(float4)*moleculeSize, cudaMemcpyHostToDevice);
	//CUDA_memcpy_to_device(&device_float4_residuePositions[memoryPosition], &host_float4_residuePositions[memoryPosition], sizeof(float4)*moleculeSize);
	//residueMeta will already be on the device
#endif
#if INCLUDE_TIMERS
	CUT_SAFE_CALL( cutStopTimer(replicaMoleculeUpdateTimer) );
#endif
	return;
}

void Replica::EonDeviceAsync()
{
#if INCLUDE_TIMERS
	replicaECUDACounter++;
	CUT_SAFE_CALL( cutStartTimer(replicaECUDATimer) );
#endif
#if CUDA_STREAMS
	// compute potential parts
	CUDA_EonDevice_async(device_float4_residuePositions, device_float4_residueMeta, device_residueCount, device_moleculeStartPositions, device_moleculeCount, device_LJPotentials, device_kernelResult, resultSize, blockSize, dataSetSize, sharedMemSize, cudaStream);
	//CUDA_EonDevice(device_float4_residuePositions, device_float4_residueMeta, device_residueCount, device_moleculeStartPositions, device_moleculeCount, device_LJPotentials, device_kernelResult,blockSize,dataSetSize);

	// write back to CPU
	//CUDA_memcpy_to_host_async(kernelResult,device_kernelResult,resultSize*resultSize*sizeof(float),cudaStream);
	cudaMemcpyAsync(kernelResult,device_kernelResult,resultSize*resultSize*sizeof(float),cudaMemcpyDeviceToHost,cudaStream);

#else
	cout << " ! Replica::EonDeviceAsync() can only be run using streams" << endl;
#endif
#if INCLUDE_TIMERS
	CUT_SAFE_CALL( cutStopTimer(replicaECUDATimer) );
#endif
}

double Replica::EonDevice()
{
#if INCLUDE_TIMERS
	replicaECUDACounter++;
	CUT_SAFE_CALL( cutStartTimer(replicaECUDATimer) );
#endif

	double result(0.0);

	CUDA_EonDevice(device_float4_residuePositions, device_float4_residueMeta, device_residueCount, device_moleculeStartPositions, device_moleculeCount, device_LJPotentials, &result,blockSize,dataSetSize,sharedMemSize);

	#if INCLUDE_TIMERS
	CUT_SAFE_CALL( cutStopTimer(replicaECUDATimer) );
#endif
	return result;
}

void Replica::FreeDevice()
{
	replicaIsOnDevice = false;  // no longer on device

	#if CUDA_STREAMS
		cudaFreeHost(host_float4_residuePositions);
	#else
		delete [] host_float4_residuePositions;
	#endif
	#if METADATA_MEMORY == TEXTURE_MEM
		freeMetaDataTexture();
	#endif
	#if POSITIONDATA_MEMORY == TEXTURE_MEM
		freePositionDataTexture();
	#endif

	delete [] host_float4_residueMeta;
	delete [] host_moleculeStartPositions;

	cudaFree(device_float4_residuePositions);
	cudaFree(device_float4_residueMeta);
	cudaFree(device_moleculeStartPositions);
	cudaFree(device_residueCount);
	cudaFree(device_moleculeCount);

#if CUDA_MC
	cudaFree(device_moleculeCenters);
	cudaFree(device_moleculeLenghts);
	cudaFree(device_translationVector);
	cudaFree(device_reverseTranslationVector);
	cudaFree(device_rotationVector);
	cudaFree(device_reverseRotationVector);
	cudaFreeHost(host_translationVector);
	cudaFreeHost(host_reverseTranslationVector);
	cudaFreeHost(host_rotationVector);
	cudaFreeHost(host_reverseRotationVector);
#endif
}

#if CUDA_MC
// rotate the molecule on device about a vector and amount
bool Replica::rotateOnDevice(int moleculeID, Vector3f rvector, float amount)
{
	//device_rotationVector;  // vector(x,y,z)|amount(w)

	host_rotationVector->x = rvector.x;
	host_rotationVector->y = rvector.y;
	host_rotationVector->z = rvector.z;
	host_rotationVector->w = amount;
	CUDA_memcpy_to_device_async(device_rotationVector,&host_rotationVector,sizeof(float4),cudaStream);

	CUDA_rotateMolecule(device_float4_residuePositions, &device_moleculeStartPositions[moleculeID], &device_moleculeLenghts[moleculeID], molecules[moleculeID].length, device_rotationVector, &device_moleculeCenters[moleculeID], cudaStream);

	// create the undo vector while its busy
	host_reverseRotationVector->x = rvector.x;
	host_reverseRotationVector->y = rvector.y;
	host_reverseRotationVector->z = rvector.z;
	host_reverseRotationVector->w = -amount;
	CUDA_memcpy_to_device_async(device_reverseTranslationVector,&host_reverseTranslationVector,sizeof(float4),cudaStream);


	return true;
}

bool Replica::translateOnDevice(int moleculeID, Vector3f translation)
{
	// send the vector to the device
	host_translationVector->x = translation.x;
	host_translationVector->y = translation.y;
	host_translationVector->z = translation.z;
	CUDA_memcpy_to_device_async(device_translationVector,&host_translationVector,sizeof(float4),cudaStream);

//  CUDA_translateMolecule (float4 *residuePositions, int *startPosition, int *moleculeLength, int moleculeLength, int *moleculeId, float4* translation, cudaStream_t stream)

	CUDA_translateMolecule(device_float4_residuePositions, &device_moleculeStartPositions[moleculeID], &device_moleculeLenghts[moleculeID], molecules[moleculeID].length, device_translationVector, &device_moleculeCenters[moleculeID], cudaStream);

	// create the undo vector while its busy
	host_reverseTranslationVector->x = -translation.x;
	host_reverseTranslationVector->y = -translation.y;
	host_reverseTranslationVector->z = -translation.z;
	CUDA_memcpy_to_device_async(device_reverseTranslationVector,&host_reverseTranslationVector,sizeof(float4),cudaStream);

	return true;
}

void Replica::cudaRollbackMutation()
{
	if (lastMutationWasATranslate)
	{
		CUDA_translateMolecule(device_float4_residuePositions, &device_moleculeStartPositions[lastMutatedMolecule], &device_moleculeLenghts[lastMutatedMolecule], molecules[lastMutatedMolecule].length, device_reverseTranslationVector,&device_moleculeCenters[lastMutatedMolecule], cudaStream);
	}
	else
	{
		CUDA_rotateMolecule(device_float4_residuePositions, &device_moleculeStartPositions[lastMutatedMolecule], &device_moleculeLenghts[lastMutatedMolecule], molecules[lastMutatedMolecule].length, device_reverseRotationVector,&device_moleculeCenters[lastMutatedMolecule], cudaStream);
	}
}
#endif // CUDA MC
#endif  // CUDA

bool Replica::sample(SimulationData *data, int current_step, float boundEnergyThreshHold, pthread_mutex_t *writeFileMutex)
{
	bool isBound = false;
	float nonCrowderPotential(0.0f);
	samplesSinceLastExchange++;


	// todo: change to a binding distance matrix and check if all are less than the cut off in the case of 3 or more proteins
	//float bindingDistance = data->bound*1.733f; // just larger than the max distance in the bounding box
	//float tmp_d;

	// in the case of all molecules being of interest then just return the count
	if (moleculeCount == nonCrowderCount)
	{
		nonCrowderPotential = potential;
	}
	else
	{
		//TODO: this must be changed to a cuda call because it is slow for large sizes
		//need to make a new version of E to accept non crowders only
		// ie: 	EonDeviceNC(&nonCrowderPotential);

		// if there are crowders and molecules of interest the only use the energy of the interesting ones
		for (int i=0;i<nonCrowderCount;i++)
		{
			for (int j=i+1;j<nonCrowderCount;j++)
			{
				if (moleculeCount != nonCrowderCount)
				{
					nonCrowderPotential += E(&molecules[i],&molecules[j]);
				}
			}
		}
	}

	if (nonCrowderPotential < boundEnergyThreshHold)
	{
		isBound = true;
		boundSamples++;
		pthread_mutex_lock(writeFileMutex);
		fprintf(data->boundConformations, "%d; %0.5f (%0.5f); %0.1f\n", current_step,nonCrowderPotential,potential,temperature);

		for (int a=0;a<moleculeCount;a++)
		{
			if (a<nonCrowderCount)
				fprintf(data->boundConformations,"0 %2d: %f %f %f %f %f %f %f\n", a, molecules[a].rotation.w,molecules[a].rotation.x,molecules[a].rotation.y,molecules[a].rotation.z,molecules[a].center.x,molecules[a].center.y,molecules[a].center.z);
			else
				fprintf(data->boundConformations,"1 %2d: %f %f %f %f %f %f %f\n", a, molecules[a].rotation.w,molecules[a].rotation.x,molecules[a].rotation.y,molecules[a].rotation.z,molecules[a].center.x,molecules[a].center.y,molecules[a].center.z);
		}


		pthread_mutex_unlock(writeFileMutex);

	}
	return isBound;
}

void Replica::save(char* filename)
{
	cout << " *** INSERT SAVECODE CODE HERE *** " << endl;
}
