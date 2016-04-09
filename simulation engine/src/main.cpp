#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <cmath>
#include <pthread.h>
#include <cuda.h>
#include <cutil.h>  // CUDA c util package
#include <cuda_runtime_api.h>
#include <cutil_inline.h>
#include <cutil_inline_runtime.h>
#include <sys/mman.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "definitions.h"
#include "AminoAcid.h"
#include "TorsionalLookupMatrix.h"
#include "Replica.h"
#include "cudaExterns.h"
#include "vector3f.h"
#include "Quaternion.h"
#include <map>
#include <unistd.h>  // include for profiling, TAU cannot find getpid() as it is an externed call

//all cuda things
#if USING_CUDA

#include "cudaExterns.h"
int cuda_blockSize = TILE_DIM;
bool auto_blockdim = true;

#endif

using namespace std;

float e = 2.71828182845904523536f;
size_t lowestEnergy;
Replica replica[REPLICA_COUNT];   // container for all replicas in the simulation
vector<uint> nonCrowderMolecules; // contains a list of each molecule/protein we want, lets us quickly determine which ones to save


int __nsleep(const struct timespec *req, struct timespec *rem)
{
	struct timespec temp_rem;
	if(nanosleep(req,rem)==-1)
		__nsleep(rem,&temp_rem);
	else
		return 1;
}

int msleep(unsigned long milisec)
{
	struct timespec req={0},rem={0};
	time_t sec=(int)(milisec/1000);
	milisec=milisec-(sec*1000);
	req.tv_sec=sec;
	req.tv_nsec=milisec*1000000L;
	__nsleep(&req,&rem);
	return 1;
}

gsl_rng  * REMCRng;			// random numbers for the REMC method to swap labels
bool exitCondition = false;
bool viewConditions = false;
bool skipsimulation = false;
int mcstepcount = 0;
int threads = THREAD_COUNT;   // 2 threads per core i think
int streams = STREAM_COUNT;

void REMCSimulation(Replica *initialReplica, argdata *parameters);
void REMCSimulationResume(Replica *initialReplica, argdata *parameters);
double DRMS(Replica* a, Replica* b);
bool getArgs(argdata * d, int argc, char **argv);
void printHelp(bool badArg);
void loadArgsFromFile(argdata * parameters, Replica *initialReplica);
void *MCthreadableFunction(void *arg);
void initSamplingFiles (argdata *args, FILE * fractionBoundFile, FILE * boundConformationsFile, FILE * acceptanceRatioFile);
void flushSamplingFiles (FILE * fractionBoundFile, FILE * boundConformationsFile, FILE * acceptanceRatioFile);
void closeSamplingFiles (argdata *args, FILE * fractionBoundFile, FILE * boundConformationsFile, FILE * acceptanceRatioFile);
void saveSimulation(const Replica *replicaState, argdata *parameters);

// All the opengl stuff
#if GLVIS
Replica *GLreplica;
#include "openglvis.h"
#endif

void saveSimulation(const Replica *replicaState, argdata *parameters)
{
	FILE * output;
	char filename[256];
	sprintf(filename,"checkpoints/%s_%d",parameters->prependageString,parameters->pid);
	output = fopen (filename,"w");
	fprintf(output,"#Saved simulation state file, prependage string, original pid\n");
	fprintf(output,"%s\n",parameters->file);
	fprintf(output,"%s_%d\n",parameters->prependageString,parameters->pid);
	fprintf(output,"%d\n",parameters->pid);
	fprintf(output,"# threads, streams, gpus, MCsteps, RESteps, sampleFrequency, sampleStartsAfter,nonCrowders, temperatureMin, temperatureMax, currentStep\n");
	fprintf(output,"%d %d %d %d %d %d %d %d %f %f %d\n",parameters->threads,parameters->streams,parameters->gpus,parameters->MCsteps,parameters->REsteps,parameters->sampleFrequency,parameters->sampleStartsAfter,parameters->nonCrowders,parameters->temperatureMin,parameters->temperatureMax,parameters->currentStep);
	fprintf(output,"#molecules\n");
	fprintf(output,"%d\n",replicaState[0].moleculeCount);
	fprintf(output,"#replicas\n");
	fprintf(output,"%d\n",parameters->replicas);
	for (int r=0;r<parameters->replicas;r++)
	{
		fprintf(output,"#replica %d, molecules: x,y,z translation. w,x,y,z rotation quaternion\n",r);
		for (int m=0;m<replicaState[0].moleculeCount;m++)
		{
			fprintf(output,"%.7f %.7f %.7f %.7f %.7f %.7f %.7f\n",replicaState[r].molecules[m].position.x,replicaState[r].molecules[m].position.y,replicaState[r].molecules[m].position.z,
					replicaState[r].molecules[m].rotation.w, replicaState[r].molecules[m].rotation.x, replicaState[r].molecules[m].rotation.y, replicaState[r].molecules[m].rotation.z);
		}
		fprintf(output,"#replica %d, stats: accept, acceptA, reject, totalAcceptReject, totalAccept, boundSamples, samplesSinceLastExchange, totalBoundSamples, totalSamples\n",r);
		fprintf(output,"%d %d %d %d %d %d %d %d %d\n", replicaState[r].accept, replicaState[r].acceptA, replicaState[r].reject, replicaState[r].totalAcceptReject, replicaState[r].totalAccept, replicaState[r].boundSamples, replicaState[r].samplesSinceLastExchange, replicaState[r].totalBoundSamples, replicaState[r].totalSamples);
	}
	fclose(output);
}

void loadSimulation(Replica *initialReplica, Replica *replicas, argdata *parameters, char *filename)
{
	FILE * input;
	input = fopen (filename,"w");
	int fscanf_return;
	argdata savedParameters;

	fscanf_return = fscanf(input,"#Saved simulation state file and prependage\n","string");
	fscanf_return = fscanf(input,"%s\n",savedParameters.file);
	fscanf_return = fscanf(input,"%s\n",savedParameters.prependageString);
	fscanf_return = fscanf(input,"# threads, streams, gpus, MCsteps, RESteps, sampleFrequency, sampleStartsAfter,nonCrowders, temperatureMin, temperatureMax, currentStep\n");
	fscanf_return = fscanf(input,"%d %d %d %d %d %d %d %d %f %f %d\n",&savedParameters.threads,&savedParameters.streams,&savedParameters.gpus,&savedParameters.MCsteps,&savedParameters.REsteps,&savedParameters.sampleFrequency,&savedParameters.sampleStartsAfter,&savedParameters.nonCrowders,&savedParameters.temperatureMin,&savedParameters.temperatureMax,&savedParameters.currentStep);

	loadArgsFromFile(parameters,initialReplica);

	fscanf_return = fscanf(input,"#molecules\n");
	int ignoreInt(0);
	fscanf_return = fscanf(input,"%d\n",&ignoreInt);
	for (int i=0;i<initialReplica->moleculeCount;i++)
	{
		fscanf_return = fscanf(input,"\n");
	}
	fscanf_return = fscanf(input,"#replicas\n");
	fscanf_return = fscanf(input,"%d\n",&parameters->replicas);
	for (int r=0;r<parameters->replicas;r++)
	{
		fscanf_return = fscanf(input,"#replica %d, molecules: x,y,z translation. w,x,y,z rotation quaternion\n",&ignoreInt);
		for (int m=0;m<replicas[0].moleculeCount;m++)
		{
			Vector3f translation;
			Quaternionf rotationf;
			fscanf_return = fscanf(input,"%f %f %f %f %f %f %f\n",&translation.x,&translation.y,&translation.z,&rotationf.w,&rotationf.x,&rotationf.y,&rotationf.z);
			Quaternion rotation = rotationf.makeQuaternion();
			rotation.normalize();
			replicas[r].molecules[m].setPosition(translation);
			replicas[r].molecules[m].setRotation(rotation);

		}
		fscanf_return = fscanf(input,"#replica %d, stats: accept, acceptA, reject, totalAcceptReject, totalAccept, boundSamples, samplesSinceLastExchange, totalBoundSamples, totalSamples\n",&ignoreInt);
		fscanf_return = fscanf(input,"%d %d %d %d %d %d %d %d %d\n", &replicas[r].accept, &replicas[r].acceptA, &replicas[r].reject, &replicas[r].totalAcceptReject, &replicas[r].totalAccept, &replicas[r].boundSamples, &replicas[r].samplesSinceLastExchange, &replicas[r].totalBoundSamples, &replicas[r].totalSamples);
	}
	fclose(input);
}



void printHelp(bool badArg)
{
	if (badArg)
		cout << " * Bad arguments" << endl;
	cout << "Usage: REMCDockingXXX -f <filename> [-c] [-h] [-v] [-0] [-t x] [-s x] [-g x] [-m x ] [-e x] [-r x] [-o x] [-b x] [-t0 x] [-t1 x] [-bx x]"<< endl;
	cout << "\t-h|--help: show this dialog" << endl;
	cout << "\t-c: Check, perform 10 reference potential sums, gpu vs cpu" << endl;
	cout << "\t-f|--file <file>: Input config file" << endl;
	cout << "\t-v|--view:        use the open GL preview of this configuration, performs no simulation" << endl;
	cout << "\t-q|--nosim:       Do everything except the simulation (for use with -v)" << endl;
	cout << "The following values override those in the config file" << endl;
	cout << "\t-t|--threads x:  The number of CPU/pthreads to use" << endl;
	cout << "\t-s|--streams x:  The number of CUDA Streams" << endl;
	cout << "\t-g|--gpus x:     The number of GPUS to use" << endl;
	cout << "\t-m|--mcsteps x:  How many MC Steps to perform per replica " << endl;
	cout << "\t-e|--resteps x:  How many replica exhanges to perform" << endl;
	cout << "\t-r|--replicas x: How many replicas" << endl;
	cout << "\t-o|--output x:   The output prefix for files created by the simulation" << endl;
	cout << "\t-b|--boundry x:  The bounding box edge length" << endl;
	cout << "\t-t0|--tmax x: 	The temperature of the highest replica" << endl;
	cout << "\t-t1|--tmin x:    The temperature of the lowest replica" << endl;
	cout << "\t-t1|--tmin x:    The temperature of the lowest replica" << endl;
	cout << "\t-bx|--blockdim x: Number of threads per CUDA block" << endl;


	exit(0);
}

bool getArgs(argdata * d, int argc, char **argv)
{
	d->resume = false;
	d->nonCrowders = 0;

	if (!d->inputFile)
	{
		d->gpus = 1;
		d->threads = THREAD_COUNT;
		d->streams = STREAM_COUNT;
		d->MCsteps = MC_STEPS;
		d->REsteps = REMC_STEPS;
		d->sampleFrequency = SAMPLE_FREQ;
		d->sampleStartsAfter = STEPS_BEFORE_SAMPLE;
		d->inputFile = false;
		d->replicas = REPLICA_COUNT;
		d->bound = BOUNDING_VALUE;
		d->temperatureMin = LOWTEMP;
		d->temperatureMax = HIGHTEMP;
	}

	memset(&d->logfile,0,256);
	sprintf(d->logfile,"output/%d_logfile",d->pid);

	if (argc <= 1)
	{
		cout << "No arguments." << endl;
		printHelp(false);
		return true;  // use default values
	}

	int i = 1;
	while (i<argc)
	{
		// threads to be used
		if (strcmp(argv[i],"-v")==0 || strcmp(argv[i],"--view")==0)
		{
			#ifndef	EnableOPENGL
				cout << "!!! This build does not support OpenGL and the -v option" << endl;
			#else
			viewConditions = true;
			#endif
			i++;
		}

		else if (strcmp(argv[i],"-q")==0 || strcmp(argv[i],"--nosim")==0)
		{
			skipsimulation = true;
			i++;
		}

#if USING_CUDA
		else if (strcmp(argv[i],"-bx")==0 || strcmp(argv[i],"--blockdim")==0)
		{
			if (i+1>=argc)
				printHelp(true);
			cuda_blockSize = atoi(argv[i+1]);
			cout << "Block size changed to: " << cuda_blockSize << endl;
			i+=2;
			auto_blockdim = false;
		}
#endif
		// threads to be used
		else if (strcmp(argv[i],"-t")==0 || strcmp(argv[i],"--threads")==0)
		{
			if (i+1>=argc)
				printHelp(true);

			d->threads = atoi(argv[i+1]);
			i+=2;
		}

		// streams to be used
		else if (strcmp(argv[i],"-s")==0 || strcmp(argv[i],"--streams")==0)
		{
			if (i+1>=argc)
				printHelp(true);
			d->streams = atoi(argv[i+1]);
			i+=2;
		}

		// gpus to be used
		else if (strcmp(argv[i],"-g")==0 || strcmp(argv[i],"--gpus")==0)
		{
			if (i+1>=argc)
				printHelp(true);

			d->gpus = atoi(argv[i+1]);
			i+=2;
		}

		//input file for molecules
		else if (strcmp(argv[i],"-f")==0 || strcmp(argv[i],"--file")==0)
		{
			if (i+1>=argc)
				printHelp(true);
			strcpy(d->file,argv[i+1]);
			d->inputFile = true;
			i+=2;
		}

		else if (strcmp(argv[i],"-m")==0 || strcmp(argv[i],"--mcsteps")==0)
		{
			if (i+1>=argc)
				printHelp(true);

			d->MCsteps = atoi(argv[i+1]);
			i+=2;
		}

		else if (strcmp(argv[i],"-e")==0 || strcmp(argv[i],"--resteps")==0)
		{
			if (i+1>=argc)
				printHelp(true);
			d->REsteps = atoi(argv[i+1]);
			i+=2;
		}

		else if (strcmp(argv[i],"-r")==0 || strcmp(argv[i],"--replicas")==0)
		{
			if (i+1>=argc)
				printHelp(true);
			d->replicas = atoi(argv[i+1]);
			i+=2;
		}

		else if (strcmp(argv[i],"-o")==0 || strcmp(argv[i],"--output")==0)
		{
			if (i+1>=argc)
				printHelp(true);
			strcpy(d->logfile,argv[i+1]);
			i+=2;
		}
		else if (strcmp(argv[i],"-b")==0 || strcmp(argv[i],"--boundry")==0)
		{
			if (i+1>=argc)
				printHelp(true);
			d->bound = atof(argv[i+1]);
			i+=2;
		}
		else if (strcmp(argv[i],"-t0")==0 || strcmp(argv[i],"--tmin")==0)
		{
			if (i+1>=argc)
				printHelp(true);
			d->temperatureMin = atof(argv[i+1]);
			i+=2;
		}
		else if (strcmp(argv[i],"-t1")==0 || strcmp(argv[i],"--tmax")==0)
		{
			if (i+1>=argc)
				printHelp(true);
			d->temperatureMax = atof(argv[i+1]);
			i+=2;
		}
		else if (strcmp(argv[i],"-h")==0 || strcmp(argv[i],"--help")==0)
		{
			printHelp(false);
		}
		else if (strcmp(argv[i],"-z")==0 || strcmp(argv[i],"--resume")==0)
		{
			#if CHECKPOINTING
			if (i+1>=argc)
				printHelp(true);
			strcpy(d->checkpointfilename,argv[i+1]);
			FILE *checkpointfile;
			checkpointfile = fopen (d->checkpointfilename,"r");
			int result = fscanf(checkpointfile,"%s\n",d->file);
			result += fscanf(checkpointfile,"%d\n",&d->pid);
			fclose(checkpointfile);
			if (result < 2)
			{
				cout << "Checkpoint load failed." << endl;
				exit(0);
			}
			d->resume = true;
			d->inputFile = true;
			i+=2;
			#else
			cout << "!!! Checkpointing is not enabled in this build, cannot resume" << endl;
			i++;
			#endif
		}
		else
		{
			i++;
		}
	}
	return false;
}

void loadArgsFromFile(argdata * parameters, Replica *initialReplica)
{
	/* use an input file
	 * files start after the "files" word in a file.
	 * format is: px py pz rx ry rz r filename
	 * meaning: position 3 floats x,y,z, rotation axis floats rx,ry,rz, and amount r float, filename is the rest of the line
	 *
	 * file can define parameters too.
	 * one per line args:
	 * threads 2  (2 threads)
	 * streams 4  (4 streams) etc....
	 * files
	 * p 1 4 7 r 1 0 0 0 pdb data/conf1/1a.pdb
	 * crowders
	 * p 1 8 9 r 1 9 2 0.5 pdb data/conf1/x.pdb
	 *
	 */
	if (parameters->inputFile)
	{
		ifstream input(parameters->file);
		if (!input.good())
		{
			cout << "Failed to open file: " << parameters->file << "\n";
			exit(0);
		}

		char line[512] = {0};
		char saveline[512];

		bool parameterSection = true; // parameters come first
		bool moleculeSection = false;
		bool crowdersSection = false;
		int lc = 0;
		// for each line in the file
		while (!input.eof())
		{
			input.getline(line,512);
			lc++;
			strcpy(saveline,line);
			char *token = strtok(line," ");

			if (token==NULL || line[0] == '#' || strlen(token)==0)
			{
				// comment or blank line
			}
			else if (strcmp(token,"files")==0)
			{
				parameterSection = false;
				moleculeSection = true;
			}
			else if (strcmp(token,"crowders")==0)
			{
				moleculeSection = false;
				crowdersSection = true;
			}
			else if (parameterSection)  // parameters
			{
				//gpus,streams,threads,mcsteps,resteps,replicas
				if (strcmp(token,"gpus")==0)
				{
					token = strtok(NULL," ");
					parameters->gpus = atoi(token);
				}
				if (strcmp(token,"streams")==0)
				{
					token = strtok(NULL," ");
					parameters->streams = atoi(token);
				}
				if (strcmp(token,"threads")==0)
				{
					token = strtok(NULL," ");
					parameters->threads = atoi(token);
				}
				if (strcmp(token,"mcsteps")==0)
				{
					token = strtok(NULL," ");
					parameters->MCsteps = atoi(token);
				}
				if (strcmp(token,"resteps")==0)
				{
					token = strtok(NULL," ");
					parameters->REsteps = atoi(token);
				}
				if (strcmp(token,"temperaturemax")==0)
				{
					token = strtok(NULL," ");
					parameters->temperatureMax = atof(token);
				}
				if (strcmp(token,"temperaturemin")==0)
				{
					token = strtok(NULL," ");
					parameters->temperatureMin = atof(token);
				}
				if (strcmp(token,"boundry")==0)
				{
					token = strtok(NULL," ");
					parameters->bound = atof(token);
				}
				if (strcmp(token,"replicas")==0)
				{
					token = strtok(NULL," ");
					parameters->replicas = atoi(token);
				}
				if (strcmp(token,"samplefrequency")==0)
				{
					token = strtok(NULL," ");
					parameters->sampleFrequency = atoi(token);
				}
				if (strcmp(token,"sampleafter")==0)
				{
					token = strtok(NULL," ");
					parameters->sampleStartsAfter = atoi(token);
				}
				if (strcmp(token,"prefix")==0)
				{
					token = strtok(NULL," ");
					strcpy(parameters->prependageString,"");
					strcpy(parameters->prependageString,token);
					#if REPULSIVE_CROWDING
						strcat(parameters->prependageString,"_repulsive");
					#else
						strcat(parameters->prependageString,"_full");
					#endif
				}
			}
			else if (moleculeSection||crowdersSection) // files
			{
				float px,py,pz,rx,ry,rz,ra;
				char *pdbfilename = new char [256];
				bool translate = true;
				int result = 0;
				if (saveline[0]=='t')
				{
					result = sscanf(saveline,"t(%f,%f,%f) r(%f,%f,%f,%f) %s", &px,&py,&pz,&rx,&ry,&rz,&ra,pdbfilename);
				}
				if (saveline[0]=='p')
				{
					result = sscanf(saveline,"p(%f,%f,%f) r(%f,%f,%f,%f) %s", &px,&py,&pz,&rx,&ry,&rz,&ra,pdbfilename);
					translate = false;
				}
				if (result<8)
				{
					cout << "Failed to parse line " << lc << ": " << saveline << endl;
				}
				else
				{
					int moleculeId = initialReplica->loadMolecule(pdbfilename);
					if (translate)
					{
						initialReplica->molecules[moleculeId].translate(Vector3f(px,py,pz));
					}
					else
					{
						initialReplica->molecules[moleculeId].setPosition(Vector3f(px,py,pz));

					}

					if (ra >= 0.000)
					{
						Vector3double v = Vector3double(rx,ry,rz);
						v.normalizeInPlace();
						initialReplica->molecules[moleculeId].rotateQ(v,ra);
					}

					if(crowdersSection)
					{
						initialReplica->molecules[moleculeId].setMoleculeRoleIdentifier(CROWDER_IDENTIFIER);
					}
					else
					{
						parameters->nonCrowders++;
					}
				}

				delete [] pdbfilename;
			}
		}
		input.close();
	}
	if (parameters->bound <= 0)
	{
		cout << "! Bounding value too small, setting equal to " << BOUNDING_VALUE << endl;
	}
	initialReplica->setBoundingValue(parameters->bound);

	if (parameters->temperatureMax < parameters->temperatureMin)
	{
		cout << "! Maximum temperature < minimum temperature, swapping " << parameters->temperatureMax << " <-> " << parameters->temperatureMin << endl;
		float tmp = parameters->temperatureMax;
		parameters->temperatureMax = parameters->temperatureMin;
		parameters->temperatureMin = tmp;
	}

	if (parameters->threads > parameters->replicas)
	{
		cout << "! Too many threads, setting equal to " << parameters->replicas << endl;
		parameters->threads = parameters->replicas;
	}
	if (parameters->threads > parameters->streams)
	{
		parameters->streams = threads;
		if (parameters->streams > 16*parameters->gpus)
		{
			parameters->streams = 16*parameters->gpus;
			if (parameters->streams > parameters->replicas)
				parameters->streams = parameters->replicas;
			cout << "! Too many streams, setting equal to " << parameters->gpus << endl;
		}
	}

	initialReplica->nonCrowderCount = parameters->nonCrowders;
	char * fileindex = new char[256];
	sprintf(fileindex,"output/%s_%d_fileindex",parameters->prependageString,parameters->pid);
	FILE * fileindexf = fopen (fileindex,"w");
	fprintf(fileindexf,"index molecule_file_path crowder(Y/N)\n");
	delete [] fileindex;

	#if OUTPUT_LEVEL > 0

	cout << "Argument data from file:" << endl;
	cout << "-------------------------------------" << endl;
	cout << "threads " << parameters->threads << endl;
	cout << "streams " << parameters->streams << endl;
	cout << "GPUs " << parameters->gpus << endl;
	cout << "mc steps " << parameters->MCsteps << endl;
	cout << "re steps " << parameters->REsteps << endl;
	cout << "replicas " << parameters->replicas << endl;
	cout << "sampling frequency (mc steps) " << parameters->sampleFrequency << endl;
	cout << "sampling starts after (mc steps) " << parameters->sampleStartsAfter << endl;
	cout << "bounding box size " << parameters->bound << endl;
	cout << "non-crowder molecules " << parameters->nonCrowders << endl;
	cout << "maximum temperature " << parameters->temperatureMax << endl;
	cout << "minimum temperature " << parameters->temperatureMin << endl;

	cout << "Loaded: "<< endl;
	cout << "-------------------------------------------------------------"<< endl;
	for (int z=0; z<parameters->nonCrowders;z++)
	{
		printf("%2d %s centered @ (%0.3f,%0.3f,%0.3f)\n",z,initialReplica->molecules[z].filename,initialReplica->molecules[z].center.x,initialReplica->molecules[z].center.y,initialReplica->molecules[z].center.z);
		fprintf(fileindexf,"%2d %s N\n",z,initialReplica->molecules[z].filename,initialReplica->molecules[z].center.x,initialReplica->molecules[z].center.y,initialReplica->molecules[z].center.z);
	}
	for (int z=parameters->nonCrowders; z<initialReplica->moleculeCount;z++)
	{
		printf("%2d %s crowder centered @ (%0.3f,%0.3f,%0.3f)\n",z,initialReplica->molecules[z].filename,initialReplica->molecules[z].center.x,initialReplica->molecules[z].center.y,initialReplica->molecules[z].center.z);
		fprintf(fileindexf,"%2d %s Y\n",z,initialReplica->molecules[z].filename);
	}
	cout << "-------------------------------------------------------------"<< endl;
	#endif
	fclose(fileindexf);
}

void initSamplingFiles (argdata * args, FILE ** fractionBoundFile, FILE ** boundConformationsFile, FILE ** acceptanceRatioFile,  FILE ** exchangeFrequencyFile)
{
	char * sfraftionBound = new char[256];
	sprintf(sfraftionBound,"output/%s_%d_fractionBound",args->prependageString,args->pid);
	*fractionBoundFile = fopen (sfraftionBound,"at"); // attempt append a file of the same name
	if (!*fractionBoundFile) *fractionBoundFile = fopen(sfraftionBound, "wt"); // create if that failed, and open and append failed (eg: file permissions r--)
	if (!*fractionBoundFile)
	{
	   printf("Cannot open/create file: %s to record fraction bound.\n",sfraftionBound);
	   return;
	}
	fprintf(*fractionBoundFile,"# Fraction Bound\n#Iteration InstantaniousAve CumulativeAve\n");
	delete [] sfraftionBound;

	char * sboundConformationsFile = new char[256];
	sprintf(sboundConformationsFile,"output/%s_%d_boundconformations",args->prependageString,args->pid);
	*boundConformationsFile = fopen (sboundConformationsFile,"a+");
	if (!*boundConformationsFile) *boundConformationsFile = fopen(sboundConformationsFile, "wt");
	if (!*boundConformationsFile)
	{
	   printf("Cannot open/create file: %s to record bound conformations.\n",sboundConformationsFile);
	   return;
	}
	fprintf(*boundConformationsFile,"# iteration; complex free energy (replica free energy); temperature;\n# molecule: rotation(w,x,y,z) position(x,y,z) 1 line per other molecule in the bound state\n");
	delete [] sboundConformationsFile;

	char * acceptanceRatioName = new char[256];
	sprintf(acceptanceRatioName,"output/%s_%d_acceptance_ratios",args->prependageString,args->pid);
	*acceptanceRatioFile = fopen (acceptanceRatioName,"a+");
	if (!*acceptanceRatioFile) *acceptanceRatioFile = fopen(acceptanceRatioName, "wt");
	if (!*acceptanceRatioFile)
	{
	   printf("Cannot open/create file: %s to record acceptance ratio.\n",acceptanceRatioName);
	   return;
	}
	fprintf(*acceptanceRatioFile,"# Iteration AcceptanceRatio\n");
	delete [] acceptanceRatioName;


	char * exchangeFrequancyName = new char[256];
	sprintf(exchangeFrequancyName,"output/%s_%d_exchange_freq",args->prependageString,args->pid);
	*exchangeFrequencyFile = fopen (exchangeFrequancyName,"a+");
	if (!*exchangeFrequencyFile) *exchangeFrequencyFile = fopen(exchangeFrequancyName, "wt");
	if (!*exchangeFrequencyFile)
	{
	   printf("Cannot open/create file: %s to record exchanges.\n",exchangeFrequancyName);
	   return;
	}
	fprintf(*exchangeFrequencyFile,"# Iteration continous_ratio, inst_ratio\n");
	delete [] exchangeFrequancyName;

	fflush(*acceptanceRatioFile);
	fflush(*fractionBoundFile);
	fflush(*boundConformationsFile);
	fflush(*exchangeFrequencyFile);

	return;
}

void flushSamplingFiles (FILE * fractionBoundFile, FILE * boundConformationsFile, FILE * acceptanceRatioFile)
{
	fflush(fractionBoundFile);
	fflush(acceptanceRatioFile);
	fflush(boundConformationsFile);
}

void closeSamplingFiles (argdata *args, FILE * fractionBoundFile, FILE * boundConformationsFile, FILE * acceptanceRatioFile, FILE * exchangeFrequencyFile)
{
	fclose(fractionBoundFile);
	fclose(boundConformationsFile);
	fclose(acceptanceRatioFile);
	fclose(exchangeFrequencyFile);
	#if OUTPUT_LEVEL > 0

	cout << ">>> Sampling files closed." << endl;
	printf("    - output/%s_%d_fractionBound\n",args->prependageString,args->pid);
	printf("    - output/%s_%d_boundconformations\n",args->prependageString,args->pid);
	printf("    - output/%s_%d_acceptance_ratios\n",args->prependageString,args->pid);
	#endif
}


void initLattice(Replica *initialReplica)
{
	/*initialReplica->reserveConiguousMoleculeArray(2);
	float interMolecularSpacing = BOUNDING_VALUE / 2 - 1;
	nonCrowderMolecules.clear();
	//initialReplica->loadMolecule("data/conf1/1a.pdb");
	//initialReplica->loadMolecule("data/conf1/1b.pdb");

	#if BOUNDING_SPHERE

	initialReplica->loadMolecule("data/application/ubq.pdb",Vector3f(-20,-20,-20),Vector3f(0,1,0),3.142);
	initialReplica->loadMolecule("data/application/uim.pdb",Vector3f(20,20,20),Vector3f(1,0,0),-3.142);

	#else

	initialReplica->loadMolecule("data/application/ubq.pdb",Vector3f(45,45,45),Vector3double(0,1,0),3.142);
	nonCrowderMolecules.push_back(uint(initialReplica->moleculeCount));
	initialReplica->loadMolecule("data/application/uim.pdb",Vector3f(60,60,60),Vector3double(1,0,0),-3.142);
	nonCrowderMolecules.push_back(uint(initialReplica->moleculeCount));


	gsl_rng *init_rotate = gsl_rng_alloc (gsl_rng_mt19937);

	//create a lattice of protein crowders
	for (int x=0;x<3;x++)
	{
		for (int y=0;y<2;y++)
		{
			for (int z=0;z<2;z++)
			{
				if(!(x==1 && y==1 && z==1))
				{
					initialReplica->loadMolecule("data/application/cspA.pdb",Vector3f(x,y,z)*interMolecularSpacing,initialReplica->createNormalisedRandomVectord(init_rotate),1.0);
					initialReplica->molecules[initialReplica->moleculeCount-1].setMoleculeRoleIdentifier(CROWDER_IDENTIFIER);
				}
			}
		}
	}
	gsl_rng_free(init_rotate);

	#endif

	cout << "Loaded: "<< endl;
	cout << "-------------------------------------------------------------"<< endl;
	for (int z=0; z<initialReplica->moleculeCount;z++)
		cout << initialReplica->molecules[z].filename << " " << initialReplica->molecules[z].moleculeRoleIdentifier << endl;
	cout << initialReplica->residueCount << " residues" << endl;
	cout << "-------------------------------------------------------------"<< endl;
	 */
	cout << "EXIT: Program requires a configuration file." << endl;
	cout << "use: REMCDockingGPU -f filename" << endl;
	exit(0);
}

float isBound(Molecule *A, Molecule *B)
{
	for (size_t residueI=0; residueI<A->residueCount;residueI++)
	{
		for (size_t residueJ=0; residueJ<B->residueCount;residueJ++)
		{
			float distance = (A->Residues[residueI].position - B->Residues[residueJ].position).magnitude();

			// if the experiment pair counts as a native contact
			if (distance < 8.0)
			{
				return 1.0f;
			}
		}
	}
	return 0.0f;
}

// calculate the root mean square between replicas
double DRMS(Replica* exp, Replica* sim)
{
	double drms = 0;
	double N = 0;
	double native_contacts = 0;
	double experimental_contacts = 0;
	double binding_contacts = 0;

	for (size_t moleculeI=0; moleculeI<exp->moleculeCount;moleculeI++)
	{
		for (size_t moleculeJ=moleculeI; moleculeJ<exp->moleculeCount;moleculeJ++)
		{
			if (moleculeI!= moleculeJ)
			{
				for (size_t residueI=0; residueI<exp->molecules[moleculeI].residueCount;residueI++)
				{
					for (size_t residueJ=0; residueJ<exp->molecules[moleculeJ].residueCount;residueJ++)
					{
						N++;
						drms += abs((exp->molecules[moleculeI].Residues[residueI].position - exp->molecules[moleculeJ].Residues[residueJ].position).magnitude() -
								(sim->molecules[moleculeI].Residues[residueI].position - sim->molecules[moleculeJ].Residues[residueJ].position).magnitude());

						// get native contacts. etc for quality stuff.

						// if the experiment pair counts as a native contact
						if ((exp->molecules[moleculeI].Residues[residueI].position - exp->molecules[moleculeJ].Residues[residueJ].position).magnitude() < 12*Angstrom)
						{
							native_contacts++;
							if ((sim->molecules[moleculeI].Residues[residueI].position - sim->molecules[moleculeJ].Residues[residueJ].position).magnitude() < 15*Angstrom)
							{
								experimental_contacts++;
							}
						}
						if ((sim->molecules[moleculeI].Residues[residueI].position - sim->molecules[moleculeJ].Residues[residueJ].position).magnitude() < 12*Angstrom)
						{
							binding_contacts++;
						}
					}
				}
			}
		}
	}
	drms /= N;
	cout << "DRMS (angstrom): " << drms/Angstrom << endl;
	cout << "Quality of interfacial residue data [0:1]: " << experimental_contacts/native_contacts << endl;
	cout << "Binding contacts discovered: " << binding_contacts << endl;
	return drms;
}

double DRMS(Molecule* a, Molecule* b)
{
	double drms = 0;
	double N = 0;

	for (size_t residueA=0; residueA<a->residueCount;residueA++)
	{
		for (size_t residueB=0; residueB<b->residueCount;residueB++)
		{
			N++;
			drms += abs((a->Residues[residueA].position - b->Residues[residueB].position).magnitude());
		}
	}
	drms /= N;
	return drms;
}

// required objects for synchronisation
pthread_mutex_t waitingThreadMutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t waitingCounterMutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t reMutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t writeFileMutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t waitingThreadCond;
pthread_cond_t waitingReplicaExchangeCond;
int waitingThreadCount;

void *MCthreadableFunction(void *arg)
{

	SimulationData *data = (SimulationData *) arg;
	long threadIndex = data->index;
	Replica *replica = data->replica;



	printf ("--- Monte-Carlo thread %d running. ---.\n",int(threadIndex+1));

	#if GLVIS
	GLreplica = &replica[0];
	#endif
	// determine the number of replicas the thread will run
	int tReplicas = int(ceil(float(data->replicaCount)/float(data->threads)));
	//to make sure we have all replicas execute, check that tReplicas*(i+1) < REPLICA_COUNT
	int replicasInThisThread = tReplicas;
	if  (tReplicas*(threadIndex+1)>data->replicaCount)
		replicasInThisThread = tReplicas - (tReplicas*data->threads)%data->replicaCount;

	if (replicasInThisThread-1+threadIndex*tReplicas > data->replicaCount)
	{
		cout << "TOO MANY THREADS. Either increase the number of replicas or decrease the thread count. " << endl;
		exit(0);
	}

	#ifdef VERY_VERBOSE
	pthread_mutex_lock( &replicaMutex );
	for (int tx=0;tx<replicasInThisThread;tx++)
	{
		cout << "+ Thread " << threadIndex << " running replica " << tx+threadIndex*tReplicas << endl;
	}
	pthread_mutex_unlock( &replicaMutex );
	#endif
	#if USING_CUDA

	//initialise cuda for use in this thread
	cuInit(0);
	cutilCheckMsg("Failed to initialise CUDA runtime.");

	cudaSetDevice(data->GPUID);
	cutilCheckMsg("Failed to pick device for the CUDA runtime.");

	// copy the LJpotentials to gpu memory in this thread context so it can access it
	float *deviceLJpTmp;
	cudaMalloc((void**)&deviceLJpTmp,sizeof(float)*AA_COUNT*AA_COUNT);
	cutilCheckMsg("Failed to allocate contact potential memory on the GPU");

	CUDA_setBoxDimension(replica->boundingValue);
	cutilCheckMsg("Failed to copy box dimensions to GPU");

	copyLJPotentialDataToDevice(deviceLJpTmp, &replica[threadIndex*tReplicas].aminoAcids);
	cutilCheckMsg("Failed to code contact potentials to device.");


	// create streams for each subsequent MCSearch
	// ensure the replica can find the lookup table
	// initialise data on the device and copy the initial batch

	#if CUDA_STREAMS

	cudaStream_t streams[16];   // reserve 16 stream slots but create only as many as needed
	// the stream/replica ration must be a whole number otherwise there will be lots of waste, ie dormant streams etc
	int streamsAvailiablePerThread = data->streams/data->threads;
	int replicasPerStream = int(ceil(float(replicasInThisThread)/float(streamsAvailiablePerThread)));
	int streamsPerThread  = replicasInThisThread/replicasPerStream;

	for (int is=0;is<streamsPerThread;is++)
	{
		cudaStreamCreate(&streams[is]);
	}
	#endif

	#if LJ_LOOKUP_METHOD == TEXTURE_MEM
		bindLJTexture(deviceLJpTmp);
	#endif

		for (int tx=0;tx<replicasInThisThread;tx++)
	{
		replica[tx+threadIndex*tReplicas].setLJpotentials(deviceLJpTmp);
		replica[tx+threadIndex*tReplicas].ReplicaDataToDevice();


		//printf("Initial potential value for replica %d (compute thread %d): %20.10f\n",int(tx+threadIndex*tReplicas) ,int(threadIndex), replica[tx+threadIndex*tReplicas].EonDevice());
		#if CUDA_STREAMS
			data->replica[tx+threadIndex*tReplicas].cudaStream = streams[tx%streamsPerThread];	// use rotation to assign replicas to streams
			replica[tx+threadIndex*tReplicas].ReserveSumSpace();								// reserve a space for the potential summation to be stored
			replica[tx+threadIndex*tReplicas].savedMolecule.reserveResidueSpace(replica[tx+threadIndex*tReplicas].maxMoleculeSize);
		#endif

	}



	#endif

	int mcstepsPerRE = data->MCsteps/data->REsteps;	// number of MC steps to do at a time

	// run the MC loop
	#if !CUDA_MC && !CUDA_STREAMS	// if no MC on device
	int mcstep = 0;
	int samples = 0;
	while (mcstep<data->MCsteps)
	{
		//cout << "Starting mc loop at " << mcstep << " steps" << endl;
		// do all MC steps before RE
		for (int tx=0;tx<replicasInThisThread;tx++)
		{
			// to sample at the correct rate split this into another set of loops
			for (int s=0;s<mcstepsPerRE/data->sampleFrequency;s++)
			{
				replica[tx+threadIndex*tReplicas].MCSearch(data->sampleFrequency);

				if (mcstep+s*data->sampleFrequency >= data->sampleStartsAfter)
				{
					// cout << "Sampling at step:" << mcstep+(s+1)*data->sampleFrequency << endl;
					// if E < -1.1844 kcal/mol then its bound
					replica[tx+threadIndex*tReplicas].sample(data,mcstep+s*data->sampleFrequency,BOUND_ENERGY_VALUE,&writeFileMutex);
				}
			}
		}
		mcstep += mcstepsPerRE;
		//cout << "thread " << threadIndex << " waits for mutex " << endl;

		// do replica exchange
		pthread_mutex_lock(&waitingCounterMutex);     			  // lock the counter
		//cout << "thread " << threadIndex << " gets mutex " << endl;

		data->waitingThreadCount[0] = data->waitingThreadCount[0]+1;
		if (data->waitingThreadCount[0]==data->threads)					// if all threads in waiting state
		{
			pthread_cond_signal(&waitingReplicaExchangeCond);
			//cout << "thread " << threadIndex << " signals parent" << endl;
		}
		if (mcstep<data->MCsteps)						// wait if another MC loop must be done
		{
			//cout << "thread " << threadIndex << " waits after " << mcstep << endl;

			pthread_cond_wait(&waitingThreadCond,&waitingCounterMutex);	// wait for replica exchange to finish.
			// NB! unlock the mutex locked upon the condition above
			// being met because this will unblock all threads such
			// that they will continue concurrently, if not unlocked
			// other threads will run sequentially
			//cout << "thread " << threadIndex << " releases mutex after wait" << endl;

		//	pthread_mutex_unlock(&waitingCounterMutex);
		}
		pthread_mutex_unlock(&waitingCounterMutex);


	} // continue MC
	#endif


	#if CUDA_STREAMS	// use streams
	int mcstep = 0;
	while (mcstep<data->MCsteps)
	{
		// to sample at the correct rate split this into another set of loops
		for (int mcx=0;mcx<mcstepsPerRE;mcx++) // at each mc step
		{
			for (int index=0;index<replicasInThisThread;index+=replicasPerStream)
			{
				for (int rps=0;rps<replicasPerStream;rps++)
				{
					// batch replicas such that no stream is shared per batch
					replica[threadIndex*tReplicas+index+rps].MCSearchMutate();
					replica[threadIndex*tReplicas+index+rps].MCSearchEvaluate();

				}
				for (int rps=0;rps<replicasPerStream;rps++)
				{
					replica[threadIndex*tReplicas+index+rps].MCSearchAcceptReject();
				}

				for (int rps=0;rps<replicasPerStream;rps++)
				{
					if (mcstep+mcx > data->sampleStartsAfter && mcstep+mcx % data->sampleStartsAfter == 0)
					{
						replica[threadIndex*tReplicas+index+rps].sample(data,mcstep+mcx*data->sampleFrequency,BOUND_ENERGY_VALUE,&writeFileMutex);
					}
				}
			}
			mcstep++;
		}

		// do replica exchange
		pthread_mutex_lock(&waitingCounterMutex);     			  // lock the counter
		data->waitingThreadCount[0] = data->waitingThreadCount[0]+1;
		if (data->waitingThreadCount[0]==data->threads)
		{
			pthread_cond_signal(&waitingReplicaExchangeCond);
		}
		if (mcstep<data->MCsteps)											// wait if another MC loop must be done
		{
			pthread_cond_wait(&waitingThreadCond,&waitingCounterMutex);		// wait for replica exchange to finish.
		}
		pthread_mutex_unlock(&waitingCounterMutex);

	} // continue MC
	#endif



#if USING_CUDA
	// sync all streams and free gpu memory
	for (int tx=0;tx<replicasInThisThread;tx++)
	{
		#if CUDA_STREAMS
		replica[tx+threadIndex*tReplicas].FreeSumSpace();
		delete [] replica[tx+threadIndex*tReplicas].savedMolecule.Residues;
		#endif
		replica[tx+threadIndex*tReplicas].FreeDevice();
	}

	#if CUDA_STREAMS
	for (int i=0;i<streamsPerThread;i++)
	{
		cudaStreamDestroy(streams[i]);
	}
	#endif

	// free the memory used by the LJ table
	#if LJ_LOOKUP_METHOD == TEXTURE_MEM
		unbindLJTexture();
	#endif
	cudaFree(deviceLJpTmp);
#endif //using cuda
	/*pthread_mutex_lock(&waitingCounterMutex);     			  	// lock the counter
	data->waitingThreadCount[0] = data->waitingThreadCount[0]+1;		// wait one last time for all threads to complete before joining them.
	if (data->waitingThreadCount[0]==data->threads)
	{
		pthread_cond_signal(&waitingReplicaExchangeCond);
	}
	pthread_mutex_unlock(&waitingCounterMutex);	  			//unlock the counter
*/
	printf (" >>> Monte-Carlo thread %d exited.\n",int(threadIndex));
	return 0;
}

// Multithreaded simulation
void REMCSimulation(Replica *initialReplica, argdata *parameters)
{
	if (!parameters->resume)
		parameters->currentStep = 0;
	int checkpointFrequency = CHECKPOINTFREQUENCY; // in resteps steps

	int threadCount = parameters->threads;
	int availiableStreams = parameters->streams;
	int availiableGpus = parameters->gpus;
	int waitingThreads = 0;
	int conformationsBound = 0;
	FILE *boundConformationsFile;
	FILE *fractionBoundFile;
	FILE *acceptanceRatioFile;
	FILE *exchangeFrequencyFile;


	#if USING_CUDA
	cudaGetDeviceCount(&availiableGpus);
	if (availiableGpus>1) cout << availiableGpus << " GPUs availiable" << endl;
	//cutilCheckMsg("");
	#endif

	if (initialReplica->moleculeCount == 0)  // make sure something is loaded
	{
		cout << "-!- No molecules loaded. -!-" << endl;
		cout << "-!- Aborting REMCSimulation. -!- " << endl;
		return;
	}

	#if INCLUDE_TIMERS
	uint RELoopTimer;
	uint MCLoopTimer;
	cutCreateTimer(&RELoopTimer);
	cutCreateTimer(&MCLoopTimer);
	#endif

	//char *startEs = new char[512];
	//char *currentEs = new char[512];
	//char *tmps = new char[64];

	float lowestPairwisePotential = 1000000000;

	REMCRng = gsl_rng_alloc (gsl_rng_mt19937);
	//strcpy (startEs,"");

	// copy the data into each replica
	// copy initial replica to other replicas and init the rngs for each
	// we can't use pthreads and CUDA at the moment, but we are going to use streams

	//basically we change the simulation to spawn N threads N = number of CPU cores

	//initialise the replicas data
	int _300kReplica = 0;

	double geometricTemperature = pow(double(parameters->temperatureMax/parameters->temperatureMin),double(1.0/double(parameters->replicas-1)));
	double geometricTranslate = pow(double(MAX_TRANSLATION/MIN_TRANSLATION),double(1.0/double(parameters->replicas-1)));
	double geometricRotation = pow(double(MAX_ROTATION/MIN_ROTATION),double(1.0/double(parameters->replicas-1)));

	for (size_t i=0;i<parameters->replicas;i++)
	{
		replica[i].copy(*initialReplica);
		#if INCLUDE_TIMERS
			replica[i].initTimers();
		#endif
		replica[i].setLabel(i);

		// changed to geometric sequence
		replica[i].setTemperature(parameters->temperatureMin * pow(geometricTemperature,int(i)));
		replica[i].setTemperatureBeta(replica[i].temperature/300.0f);

		replica[i].setRotateStep(MIN_ROTATION * pow(geometricRotation,int(i)));
		replica[i].setTranslateStep(MIN_TRANSLATION * pow(geometricTranslate,int(i)));
		printf ("Replica %d %.3f %.3f %.3f\n",int(i),replica[i].temperature,replica[i].translateStep,replica[i].rotateStep);


		// note which replica is the 300K replica for sampling
		// if there is no replica at 300K then choose the closest one.
		if ( abs(replica[i].temperature-300.0f) < abs(replica[_300kReplica].temperature-300.0f))
			_300kReplica = i;
		replica[i].initRNGs();
		replica[i].threadCount = parameters->threads;
	}

	/*#if CHECKPOINTING
	if (parameters->resume)
	{
		FILE *checkpointfile;
		checkpointfile = fopen (parameters->checkpointfilename,"r");
		fscanf(checkpointfile,"%s\n",parameters->file);
		fscanf(checkpointfile,"%d\n",&parameters->pid);
		fscanf(checkpointfile,"%s\n",parameters->prependageString);
		fscanf(checkpointfile,"%d\n",&parameters->currentStep);

		for (int i=0;i<parameters->replicas;i++)
		{
			// fb ar
			fscanf(checkpointfile,"%d %d\n",&replica[i].totalBoundSamples,&replica[i].totalSamples);
			fscanf(checkpointfile,"%d %d %d\n",&replica[i].acceptA,&replica[i].accept,&replica[i].reject);
			for (int m=0;m<replicaState[i].moleculeCount;m++)
			{
				Quaternion qi;
				Vector3f vi;
				fscanf(checkpointfile,"    %f %f %f %f %f %f %f\n",&vi.x,&vi.y,&vi.z, &qi.w, &qi.x, &qi.y, &qi.z);
				qi.normalize();
				replica[i].molecules[m].rotation = qi.rotate(replica[i].molecules[m].rotation);
				replica[i].molecules[m].setPosition(vi);
			}
		}

		fclose(checkpointfile);
	}
	#endif*/

	initSamplingFiles(parameters,&fractionBoundFile,&boundConformationsFile,&acceptanceRatioFile,&exchangeFrequencyFile);


	// round the ~300K replica to 300K to ensure there is at least one T/T_i == 1
	//replica[_300kReplica].temperature = 300.0f;

	// create a lookup map for replica temperatures to allow for sampling temperatures and in place exchanges
	// samples are by temperature not replica index, hence temperatures get shuffled about in replica[],
	// so use replicaTemperaturesmap[] to track where a specific temperature is

	typedef map<float, int> KRmap;
	KRmap replicaTemperatureMap;
	typedef map<int,float> RKmap;
	RKmap temperatureMap;


	// build a map of lookups for temperature lookups while sampling
	for (size_t i=0;i<parameters->replicas;i++)
	{
		// lookup of where replica with temperature x is
		replicaTemperatureMap[replica[i].temperature] = i;
		// lookup of which replica has temperature x
		temperatureMap[i] = replica[i].temperature;
	}

	fprintf(fractionBoundFile,"Iteration:  ");
	fprintf(acceptanceRatioFile,"Iteration:  ");

	for(KRmap::const_iterator iterator = replicaTemperatureMap.begin(); iterator != replicaTemperatureMap.end(); ++iterator)
	{
		float temperature = iterator->first;
		fprintf(fractionBoundFile,  "%0.1fKi %0.1fKc ",temperature,temperature);
		fprintf(acceptanceRatioFile,"%0.1fKi %0.1fKc ",temperature,temperature);
	}

	fprintf(fractionBoundFile,"	\n");
	fprintf(acceptanceRatioFile,"\n");


	/*#if OUTPUT_LEVEL > 1
	for (int i=0;i<parameters->replicas;i++)
	{
		sprintf (tmps,"%3d: %1.1e",i,replica[i].potential);
		strcat(startEs,tmps);
	}
	cout << "start energies (kcal/mol): " << startEs << endl;
	#endif*/

	pthread_t *thread = new pthread_t[threadCount];;
	//For portability, explicitly create threads in a joinable state
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	pthread_cond_init (&waitingThreadCond, NULL);
	pthread_cond_init (&waitingReplicaExchangeCond, NULL);
	pthread_mutex_init(&waitingCounterMutex, NULL);
	pthread_mutex_init(&waitingThreadMutex, NULL);
	pthread_mutex_init(&reMutex, NULL);

	SimulationData *data = new SimulationData[threadCount];

	cout << cudaGetErrorString(cudaGetLastError()) << endl;


	#if INCLUDE_TIMERS
		CUT_SAFE_CALL( cutStartTimer(RELoopTimer) );
	#endif

	//if there are multiple gpus then associate here
	if (parameters->gpus > availiableGpus)
		parameters->gpus = availiableGpus;

	int gpuid = 0; // the first device, increment per thread


	int mcstepsPerRE = parameters->MCsteps/parameters->REsteps;	// number of MC steps to do at a time
	int totalMCsteps = parameters->MCsteps;

	// if samping does not divide into steps in this loop
	if (mcstepsPerRE % parameters->sampleFrequency > 0)// && mcstepsPerRE > data->sampleFrequency)
	{
		if (mcstepsPerRE > parameters->sampleFrequency)
		{
			// fit sampling frequency such that its a divisor or equal to the mc steps we do before RE
			while (mcstepsPerRE % parameters->sampleFrequency != 0)
				parameters->sampleFrequency++;
			cout << "-!- CHANGED: sample frequency changed as it does divide into the replica exchange frequency. -!-" << endl;
			cout << "-!- CHANGED: sample frequency <- "<< parameters->sampleFrequency <<". -!-" << endl;

		}
		if (mcstepsPerRE < parameters->sampleFrequency)
		{
			parameters->sampleFrequency = mcstepsPerRE;
			cout << "-!- CHANGED: sample frequency too long for MC loop. -!-" << endl;
			cout << "-!- CHANGED: sample frequency <- "<< mcstepsPerRE <<". -!-" << endl;
		}
	}

	#if CHECKPOINTING
	if (checkpointFrequency % mcstepsPerRE != 0) //if it doesnt occur during a RE
	{
		if (checkpointFrequency < mcstepsPerRE)
		{
			checkpointFrequency = mcstepsPerRE;
		}
		else
		{
			checkpointFrequency = ceil( float(checkpointFrequency) / float(mcstepsPerRE)) * mcstepsPerRE; // make it so
		}
	}
	int stepsUntilCheckpoint = checkpointFrequency;
	printf ("--- The simulation will checkpoint every %d monte carlo steps ---.\n",checkpointFrequency);
	#endif

	printf ("--- Staring simulation ---.\n");


	pthread_mutex_lock(&waitingCounterMutex);

	//TODO: assign streams as a function of the number of GPUs and threads

	//spawn N threads to do the Monte-Carlo mutations
	for (int i=0;i<threadCount;i++)
	{
		data[i].replica = replica;
		data[i].replicaCount = parameters->replicas;
		data[i].index = i;
		data[i].threads = threadCount;
		data[i].streams = availiableStreams;///threadCount; //determines the number of streams available
		data[i].MCsteps = parameters->MCsteps;
		data[i].REsteps = parameters->REsteps;
		data[i].sampleFrequency = parameters->sampleFrequency;
		data[i].sampleStartsAfter = parameters->sampleStartsAfter;
		data[i].bound = parameters->bound;
		data[i].streams = parameters->streams;
		data[i].waitingThreadCount = &waitingThreads;
		data[i].conformationsBound = &conformationsBound;
		data[i].fractionBound = fractionBoundFile;
		data[i].boundConformations = boundConformationsFile;

		// assign gpus in rotation per thread, t0 = gpu0, t1 = gpu1 etc
		// % #gpus so they share if threads > gpus
		// will perform best if threads:gpus = 1:1
		data[i].GPUID = i % parameters->gpus;
		cout << "Assign thread " << i << " to GPU " << data[i].GPUID << endl;

		//start the MC loops
		pthread_create(&thread[i], &attr, MCthreadableFunction, (void*)&data[i]);
		//MCthreadableFunction(&data[i]);
	}

	int exchanges = 0;  // number of exchanges performed
	int tests = 0; 		// number of replica exchange tests
	int totalExchanges = 0; // total number of exchanges
	int totalTests = 0; // total number of exchanges

	int offset = 0;  // RE offset
	int steps = 0;	 // # re steps performed
	while (++steps < parameters->REsteps)  // until enough steps taken
	{
		#if INCLUDE_TIMERS
		CUT_SAFE_CALL( cutStartTimer(MCLoopTimer) );
		#endif

		pthread_cond_wait(&waitingReplicaExchangeCond,&waitingCounterMutex);

		#if INCLUDE_TIMERS
		CUT_SAFE_CALL( cutStopTimer(MCLoopTimer) );
		#endif

		//if the sampling has started
		if (parameters->MCsteps/parameters->REsteps*(steps-1) >= parameters->sampleStartsAfter)
		{
			// do fraction bound and acceptance ratio

			fprintf(fractionBoundFile,"%9d: ",parameters->MCsteps/parameters->REsteps*steps);
			fprintf(acceptanceRatioFile,"%9d: ",parameters->MCsteps/parameters->REsteps*steps);

			for(RKmap::const_iterator iterator = temperatureMap.begin(); iterator != temperatureMap.end(); ++iterator)
			{
				int index = replicaTemperatureMap[iterator->second];

				// fraction bound
				float fractionBound = float(replica[index].boundSamples)/max(1.0f,float(replica[index].samplesSinceLastExchange));
				replica[index].totalBoundSamples += replica[index].boundSamples;
				replica[index].totalSamples += replica[index].samplesSinceLastExchange;
				replica[index].boundSamples = 0;
				replica[index].samplesSinceLastExchange = 0;
				float accumulativeFractionBound = float(replica[index].totalBoundSamples)/float(replica[index].totalSamples);
				fprintf(fractionBoundFile,"| %6.4f %6.4f ",fractionBound,accumulativeFractionBound);

				// acceptance rate
				float acceptanceRatio =  float(replica[index].acceptA + replica[index].accept)/float(replica[index].acceptA + replica[index].accept + replica[index].reject);
				replica[index].totalAccept += replica[index].acceptA+replica[index].accept ;
				replica[index].totalAcceptReject += replica[index].acceptA+replica[index].accept+replica[index].reject ;
				fprintf(acceptanceRatioFile,"| %6.4f %6.4f ",acceptanceRatio , float(replica[index].totalAccept)/float(replica[index].totalAcceptReject));
				replica[index].acceptA=0 ;
				replica[index].accept=0 ;
				replica[index].reject=0 ;
			}

			fprintf(fractionBoundFile,"\n");
			fflush(fractionBoundFile);
			fprintf(acceptanceRatioFile,"\n");
			fflush(acceptanceRatioFile);
		}

		int i = offset;  // replica exchange offset
		while (i < parameters->replicas-1) // loop through the replicas in pairs
		{
			int j = i+1; // replica to swap with current replica

			// i and j represent temperatures, hence we need to map temperature to position
			// z = temperatureMap[x] -> temperature of position i
			// replicaTemperatureMap[z] -> where the replica with temperature z actually is

			int t_i = replicaTemperatureMap[temperatureMap[i]];
			int t_j = replicaTemperatureMap[temperatureMap[j]];


			double delta = (1.0/replica[t_i].temperature - 1.0/replica[t_j].temperature)*(replica[t_i].potential - replica[t_j].potential)*(4184.0f/Rgas);
			if (gsl_rng_uniform(REMCRng) < min(1.0,exp(delta)))
			{
				replica[t_i].exchangeReplicas(replica[t_j]);

				replicaTemperatureMap[replica[t_i].temperature] = t_i;
				replicaTemperatureMap[replica[t_j].temperature] = t_j;

				exchanges++; // sampling
			}
			i += 2; // compare the next two neighbours
			tests++;// tests performed samplng
		}
		offset = 1-offset; // switch the swap neighbours in the replica exchange

		parameters->currentStep = mcstepsPerRE*steps;

		#if CHECKPOINTING
		stepsUntilCheckpoint -= mcstepsPerRE; // decrement the counter

		if (stepsUntilCheckpoint<=0) //if enough steps have passed to checkpoint
		{
			stepsUntilCheckpoint = checkpointFrequency; // reset the counter
			// if it is a checkpoint place
			saveSimulation(replica,parameters);
		}
		#endif

		waitingThreads = 0;
		pthread_cond_broadcast(&waitingThreadCond);

		#if GLVIS
		GLreplica = &replica[_300kReplica];
		GlutDisplay();
		#endif
		#if OUTPUT_LEVEL > 0
		cout << "Replica Exchange step " << steps << " of " << parameters->REsteps << " complete" << endl;
		cout.flush();
		#endif

		totalExchanges += exchanges;
		totalTests += tests;

		fprintf(exchangeFrequencyFile,"%10d %6.4f %6.4f\n",steps, float(totalExchanges)/float(totalTests), float(exchanges)/float(tests));
		tests = 0;
		exchanges = 0;

	}

	//printf ("Replica Exchange step %d skipped as it has no effect\n",steps);

	#if INCLUDE_TIMERS
	CUT_SAFE_CALL( cutStartTimer(MCLoopTimer) );
	#endif


	pthread_cond_broadcast(&waitingThreadCond);
	pthread_mutex_unlock(&waitingCounterMutex);  // release the mutex so MC threads can continue.


	printf ("--- Replica Exchanges Complete.---\n");
	printf ("--- Waiting for threads to exit. ---\n");

	// join the threads that have finished
	for (int i=0;i<threadCount;i++)
		pthread_join(thread[i],NULL);

	#if INCLUDE_TIMERS
		CUT_SAFE_CALL( cutStopTimer(MCLoopTimer) );
	#endif

	printf ("--- All threads complete.---\n");
	fprintf(fractionBoundFile,"%9d: ",parameters->MCsteps);
	fprintf(acceptanceRatioFile,"%9d: ",parameters->MCsteps);

	//final fraction bound
	for(RKmap::const_iterator iterator = temperatureMap.begin(); iterator != temperatureMap.end(); ++iterator)
	{
		int index = replicaTemperatureMap[iterator->second];

		// fraction bound
		float fractionBound = float(replica[index].boundSamples)/float(replica[index].samplesSinceLastExchange);
		replica[index].totalBoundSamples += replica[index].boundSamples;
		replica[index].totalSamples += replica[index].samplesSinceLastExchange;
		replica[index].boundSamples = 0;
		replica[index].samplesSinceLastExchange = 0;
		float accumulativeFractionBound = float(replica[index].totalBoundSamples)/float(replica[index].totalSamples);
		fprintf(fractionBoundFile,"| %6.4f %6.4f  ",fractionBound,accumulativeFractionBound);

		// acceptance rate
		float acceptanceRatio =  float(replica[index].acceptA)/float(replica[index].acceptA + replica[index].accept + replica[index].reject);
		fprintf(acceptanceRatioFile,"| %6.4f %6.4f ",acceptanceRatio , float(replica[index].totalAccept)/float(replica[index].totalAcceptReject));


	}
	fprintf(fractionBoundFile,"\n");
	fflush(fractionBoundFile);
	fprintf(acceptanceRatioFile,"\n");
	fflush(acceptanceRatioFile);

	printf ("--- Simulation finished.---.\n\n");


	#if INCLUDE_TIMERS
		CUT_SAFE_CALL( cutStopTimer(RELoopTimer) );
		printf("Simulation Timers\n");
		printf("MC Loop:   Tot  %10.5f ms  Ave %10.5fms (%d steps, %d replicas, %d threads, %d streams)\n"  	,cutGetTimerValue(MCLoopTimer),cutGetTimerValue(MCLoopTimer)/float(parameters->REsteps),parameters->MCsteps,parameters->replicas,threadCount,availiableStreams);
		printf("Simulation:     %10.5f ms  (%d exchanges)\n"	,cutGetTimerValue(RELoopTimer),parameters->REsteps);
		cutDeleteTimer(RELoopTimer);
		cutDeleteTimer(MCLoopTimer);
	#endif

	pthread_mutex_lock(&writeFileMutex);
	closeSamplingFiles(parameters,fractionBoundFile,boundConformationsFile,acceptanceRatioFile,exchangeFrequencyFile);
	pthread_mutex_unlock(&writeFileMutex);


	// Clean up/
	pthread_attr_destroy(&attr);
	pthread_mutex_destroy(&waitingCounterMutex);
	pthread_mutex_destroy(&waitingThreadMutex);
	pthread_cond_destroy(&waitingThreadCond);
	pthread_cond_destroy(&waitingReplicaExchangeCond);
	gsl_rng_free(REMCRng);
	delete [] data;
	delete [] thread;
	return;
}



int main(int argc, char **argv)
{
	int sysreturn;
	sysreturn = system("mkdir -p checkpoints");
	sysreturn = system("mkdir -p output");

	cout << "Compiled with:" << endl;
	#ifdef EnableOPENGL
		cout << "  OpenGL support" << endl;
	#endif
	#if EnableCUDA
		cout << "  CUDA support" << endl;
	#endif
	#if CUDA_STREAMS
		cout << "  Asynchronous GPU calls (CUDA capability 1.1+ required)" << endl;
	#endif

	#if COMPENSATE_KERNEL_SUM
		cout << "  Kahan summation in kernels" << endl;
	#endif

	// get global options for the simulation
	argdata parameters;
	parameters.pid = int(getpid());
	parameters.inputFile = false;
	memset(parameters.prependageString,0,256);


	nonCrowderMolecules.clear();
	bool use_defaults = getArgs(&parameters,argc,argv);
	#if CHECKPOINTING
	if (!parameters.resume)
	{
		sprintf(parameters.checkpointfilename,"checkpoints/checkpoint.%s_%d",parameters.prependageString,parameters.pid);
	}
	else
	{
		cout << "The simulation will resume from the checkpoint file: " << parameters.checkpointfilename << endl;
		cout << "Original input: " << parameters.file << endl;
	}
	#endif

	cout.precision(8);

	AminoAcids aminoAcidData;

	#if OUTPUT_LEVEL > 0
	cout << "Loading amino acid data. " << AMINOACIDDATASOURCE;
	#endif
	aminoAcidData.loadAminoAcidData(AMINOACIDDATASOURCE);
	#if OUTPUT_LEVEL > 0
	cout << "Loading pair lookup table. " << LJPDSOURCE << endl;
	#endif
	aminoAcidData.loadLJPotentialData(LJPDSOURCE);

#ifdef displayLJpotentials
	printPotentials(&aminoAcidData);
#endif

	#if USING_CUDA
	cout << "Initialising CUDA. "<< endl ;
	cuInit(0);
	cudaInfo();

	#if OUTPUT_LEVEL > 0
	cout << "CUDA parameters and options for this run:\n" ;
	cout << "-----------------------------------------\n" ;
	cout << "Tile Size " << TILE_DIM << endl;
	cout << "LJ lookup memory type: ";
	#if LJ_LOOKUP_METHOD == SHARED_MEM
	cout << "Shared" << endl;
	#elif LJ_LOOKUP_METHOD == CONST_MEM
	cout << "Constant" << endl;
	#elif LJ_LOOKUP_METHOD == GLOBAL_MEM
	cout << "Global" << endl;
	#elif LJ_LOOKUP_METHOD == TEXTURE_MEM
	cout << "Texture" << endl;
	#endif
	cout << "-----------------------------------------\n" ;
	#endif

	if (strcmp(argv[argc-1],"-c")==0)
	{
		cout << "performing check..." << endl;
	  // get values for the example conformations
		Replica exampleReplicas[10];
		char *egnames = new char[60];
		float * ljp_t;
		cudaMalloc((void**)&ljp_t,LJArraySize);
		cutilCheckMsg("Failed to cudaMalloc");
		copyLJPotentialDataToDevice(ljp_t,&aminoAcidData);
		// set box dimensions
		float testboxdim(118.4f);
		CUDA_setBoxDimension(testboxdim);

		#if LJ_LOOKUP_METHOD == TEXTURE_MEM
			bindLJTexture(ljp_t);
			//bindLJTexture2D(ljp_t);
		#endif
		//cout << "Reference Conformation energies:" << endl;
		//printf ("U           LJ         DH\n");
		for (int i = 1;i<=10;i++)
		{
			exampleReplicas[i-1].aminoAcids = aminoAcidData;
			exampleReplicas[i-1].label = i;
			exampleReplicas[i-1].setBoundingValue(testboxdim);
			sprintf(egnames,"data/conf%d/1a.pdb",i);
			exampleReplicas[i-1].loadMolecule(egnames);
			sprintf(egnames,"data/conf%d/1b.pdb",i);
			exampleReplicas[i-1].loadMolecule(egnames);
			exampleReplicas[i-1].initTimers();
			uint cpu_E_timer;
			exampleReplicas[i-1].E();

			//printf ("%10.6f %10.6f %10.6f\n",exampleReplicas[i-1].potential,exampleReplicas[i-1].E_LJ,exampleReplicas[i-1].E_DH);
			//show that cuda returns the same energys for each
			// make sure the replicas know where to find the data, as the use it for calling cuda kernels
			exampleReplicas[i-1].setDeviceLJPotentials(ljp_t);
			exampleReplicas[i-1].setBlockSize(TILE_DIM);
			exampleReplicas[i-1].ReplicaDataToDevice();
			///cout << "CPU execution time: " << cutGetTimerValue(cpu_E_timer)/100.0f << "ms"<<  endl;
			uint gpu_E_timer;
			double r = exampleReplicas[i-1].EonDevice();;

			double cpu_ee = exampleReplicas[i-1].E();
			printf ("CPU,GPU,err: %24.20f %24.20f %24.20f\n",cpu_ee,r,abs(cpu_ee-r));
					//exampleReplicas[i-1].printTimers();
			exampleReplicas[i-1].FreeDevice();

		}
		cudaFree(ljp_t);
		#if LJ_LOOKUP_METHOD == TEXTURE_MEM
			unbindLJTexture();
		#endif
		cout.flush();
		exit(0);
	} // end of check
	#endif

	Replica initialReplica;
	initialReplica.setAminoAcidData(aminoAcidData);
	initialReplica.setBoundingValue(BOUNDING_VALUE);
	initialReplica.reserveConiguousMoleculeArray(30);
	//if there is an input file load its contents here
	if (parameters.inputFile)
	{
		loadArgsFromFile(&parameters,&initialReplica);
		getArgs(&parameters,argc,argv); // second pass to override any variables if doing performance tests
	}
	else
	{
		parameters.nonCrowders = 2;
		initLattice(&initialReplica);
	}

	cout << "Loaded: " << initialReplica.residueCount << " residues in " << initialReplica.moleculeCount << " molecules:\n";

	for (int i=0;i<initialReplica.moleculeCount;i++)
	{
		printf("%2d: %3d %12.3f A^3 %s\n",i , initialReplica.molecules[i].residueCount, initialReplica.molecules[i].getVolume(),initialReplica.molecules[i].filename);
	}

	#if INCLUDE_TIMERS
		initialReplica.initTimers();
	#endif

	Replica tmp;

	float p = initialReplica.E();
	initialReplica.potential = p;
	cout << initialReplica.countpairs() << " residue interaction pairs." << endl;
	printf("CPU initial energy value: %26.23f kcal/mol\n", initialReplica.potential);


	#if USING_CUDA
	if (auto_blockdim)
		(initialReplica.residueCount < 1024) ? cuda_blockSize = 32 : cuda_blockSize = 64;

	// set box size
	initialReplica.setBlockSize(cuda_blockSize);
	initialReplica.setBoundingValue(parameters.bound);
	// set box dimensions
	CUDA_setBoxDimension(parameters.bound);


	#if CUDA_STREAMS
		cudaStreamCreate(&initialReplica.cudaStream);
	#endif
	initialReplica.ReplicaDataToDevice();

	float * ljp_tmp;
	cudaMalloc((void**)&ljp_tmp,LJArraySize);
	copyLJPotentialDataToDevice(ljp_tmp,&aminoAcidData);
	initialReplica.setDeviceLJPotentials(ljp_tmp);


	#if LJ_LOOKUP_METHOD == TEXTURE_MEM
		bindLJTexture(initialReplica.device_LJPotentials);
	#endif


	double GPUE = initialReplica.EonDevice();
	initialReplica.potential = GPUE;

	printf("GPU initial energy value: %26.23f kcal/mol\n", GPUE);
	printf("Absolute error:           %26.23f kcal/mol\n", abs(GPUE-p));
	printf("Relative error:           %26.23f kcal/mol\n", abs(p-GPUE)/abs(p));



	#if CUDA_STREAMS
		cudaStreamSynchronize(initialReplica.cudaStream);
		cudaStreamDestroy(initialReplica.cudaStream);
	#endif

	if (abs(p-GPUE)/abs(p)>0.01)
	{
		cout << "!\n! WARNING: Inconsistency between GPU and CPU result, check simulation \n!" << endl;
	}
	initialReplica.FreeDevice();
	cutilCheckMsg("Error freeing device");
	#if LJ_LOOKUP_METHOD == TEXTURE_MEM
		unbindLJTexture();
		cutilCheckMsg("Error freeing texture");
	#endif
	cudaFree(initialReplica.device_LJPotentials);
	cudaFree(ljp_tmp);

	//cutilCheckMsg("Last Cuda Error");

	#endif
/*#if INCLUDE_TIMERS
	initialReplica.printTimers();
#endif*/

	// Collision tests
/*	GlutDisplay();
	sleep(5);

    replica[0].molecules[0].translate(Vector3f(-25,0,0));
    replica[0].molecules[1].translate(Vector3f(25,0,0));

	for (int i = 0; i < 100; i++)
    {
        replica[0].molecules[0].translate(Vector3f(1,0,0));
        GlutDisplay();
	sleep(1);

        float separationDistance = (replica[0].molecules[1].center - replica[0].molecules[0].center).magnitude();
        replica[0].E();
        cout << separationDistance << " " << replica[0].potential << endl;
    }

    // End: Collision tests
*/

	if (initialReplica.nonCrowderCount < initialReplica.moleculeCount)
	{
		#if REPULSIVE_CROWDING
			cout << "-!- Crowding modelled using: u(r) = (6/r)^ 6." << endl;
		#else
			cout << "-!- Crowding is modelled using the full potential calculations." << endl;
		#endif
	}

	cout.flush();

	if (viewConditions)
	{
		#if GLVIS
			glutInit(&argc, argv);
			camera.setPosition(Vector3f(-15,15,15),Vector3f(1,-1,-1),Vector3f(0,1,0));
			char windowName[64] = {"REMC Protein Docker"};
			GlutInit(WIDTH,HEIGHT,windowName);
			gl_replicaCount = parameters.replicas;
			gl_boundingValue = int(parameters.bound);
			GLreplica = &initialReplica;
		#endif
	}

	if (!skipsimulation)
	{
		#if OUTPUT_LEVEL > 0
		cout << "Beginning simulation" << endl;
		#endif

		// need at most REPLICA_COUNT threads
		cout << "Output files will be prefixed by " << parameters.prependageString << "_" << parameters.pid << endl;
		REMCSimulation(&initialReplica,&parameters);
		#if OUTPUT_LEVEL > 0
		//printf("%d %0.20f\n",replica[0].countpairs(),replica[0].acc_err/float(parameters.MCsteps));
		cout << endl << "Simulation done" << endl;
		#endif

		for (size_t i=0;i<parameters.replicas;i++)
		{
			#if INCLUDE_TIMERS
			cout << "Replica " << i << " timers" << endl;
			replica[i].printTimers();
			#endif
		}
	}

	if (viewConditions)
	{
		#if GLVIS
		GLreplica = &initialReplica;
		cout << "Entering free gl viewing mode." << endl;
		glutMainLoop();
		#endif
	}

	if (!skipsimulation)
	{
		for (size_t i=0;i<parameters.replicas;i++)
		{
			replica[i].freeRNGs();
			delete [] replica[i].contiguousResidues;
		}
	}

	cout << "Finished." << endl;
	return 0;
}
