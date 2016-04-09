#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "definitions.h"
#include "AminoAcid.h"
#include "Replica.h"
#include "drms.h"
#include <unistd.h>
#include <sys/stat.h>


using namespace std;

int main(int argc, char **argv)
{
	cout << "usage: pdbgen <indexfile> <boxdim> [0/1: reverse drms cal] [1: outputcrowders crowders]" << endl;


	AminoAcids aminoAcidData;

	cout << "Loading amino acid data. ";
	aminoAcidData.loadAminoAcidData(AMINOACIDDATASOURCE);
	cout << "Done." << endl;
	cout << "Loading pair lookup table. " << endl;
	aminoAcidData.loadLJPotentialData(LJPDSOURCE);
	cout << "Done." << endl;

	vector<Replica> Conformations;
	Replica initalReplica;
	initalReplica.setAminoAcidData(aminoAcidData);
	initalReplica.reserveConiguousMoleculeArray(30);

	Replica reference;
	reference.setAminoAcidData(aminoAcidData);

	char *indexfilename = new char[256];
	char *savefilename = new char[256];
	char *datafilename = new char[256];
	char *drmsEfilename = new char[256];

	float boundary(atof(argv[2]));	

	cout << "File: " << argv[1] << endl;
	cout << "Bounding box size: " << boundary << "A**3" << endl; 
	
	bool printCrowders(false);

	bool reversedrms=false;
	if (argc > 3)
	{
		if (argv[3][0]=='1')
			reversedrms = true;
	}

	if (argc > 4)
	{
		if (argv[4][0]=='1')
			printCrowders = true;
	}

	strcpy(indexfilename,argv[1]);

	char * strstrresult = strstr(indexfilename,"fileindex");
	if (strstrresult==NULL)
	{
		cout << "You must provide and index file as the input" << endl;
		exit(0);
	}

	ifstream indexfile(indexfilename);
	if (!indexfile.good())
	{
		cout << "Failed to open file: " << indexfilename << endl;
		exit(0);
	}

	char line[512] = {0};

	int loadcount = 0;
	float max_d = 0;

	while (!indexfile.eof())
	{
		indexfile.getline(line,512);

		int index;
		char crowder;
		int result = sscanf(line,"%d %s %c",&index,datafilename,&crowder);
		if (result==3)
		{
			if (crowder=='N' || (crowder=='Y' && printCrowders)) 
			{
				cout << "load : " << datafilename << endl;
				initalReplica.loadMolecule(datafilename);
				reference.loadMolecule(datafilename);
				if (loadcount > 0)
				{
					Vector3f c = initalReplica.molecules[initalReplica.moleculeCount-1].position;
					for(int i=0;i<initalReplica.molecules[initalReplica.moleculeCount-1].residueCount;i++)
					{
						Vector3f r = initalReplica.molecules[initalReplica.moleculeCount-1].Residues[i].position;
						max_d = max((c-r).magnitude(),max_d);
					}
				}
				loadcount++;
			}
		}
	}
	cout << "finished loading" << endl;
	indexfile.close();

	memset(datafilename,0,256);
	strncpy(datafilename,indexfilename,strlen(indexfilename)-strlen("fileindex"));
	strcat(datafilename,"boundconformations");


	memset(drmsEfilename,0,256);
	strncpy(drmsEfilename,indexfilename,strlen(indexfilename)-strlen("fileindex"));
	strcat(drmsEfilename,"drms.dat");



	int iteration;
	float energy;   // system energy
	float temperature;
	float energyC;  // complex energy
	int samples = 0;

	ifstream datafile(datafilename);
	if (!datafile.good())
	{
		cout << "Failed to open file: " << datafilename << endl;
		exit(0);
	}
	else
	{
		cout << "loading data: "<< datafilename << endl;
	}


	//ofstream drmslog(drmsEfilename);
	//drmslog.precision(7);

	//first 2 lines contain nothing
	datafile.getline(line,512);
	datafile.getline(line,512);

	// next lines are in sets of the number of non crowders
	// 501000; -9.18593; 5.3; 297.6
	// 0 0: -0.258606 -0.929792 -0.256388 -0.053619 7.369041 103.941734 91.936638
	// 0 1: -0.055370 -0.464538 -0.429429 0.772482 7.673889 89.995583 96.419502
	char savename[256];

	int instanceNo = 0;
	while (!datafile.eof())
	{
		datafile.getline(line,512);
		//cout << "line: "<< line << endl;

		//int result = sscanf(line,"%d; %f; %f; %f;",&iteration,&energy,&distance,&temperature);
		int result = sscanf(line,"%d; %f (%f); %f;",&iteration,&energy,&energyC,&temperature);
		//int result = sscanf(line,"%d; %f; %f;",&iteration,&energy,&temperature);
		initalReplica.temperature = temperature;
		initalReplica.potential = energy;

		if (abs(temperature-300.0f)<1.0f && iteration>10000000) // if the instance is at 300K and more than 1m iterations are done
		{
			Quaternion rotation0;
			Vector3f vrotation0;
			Vector3f position0;

			Replica instance;
			instance.copy(initalReplica);
			
			for (int i=0;i<instance.moleculeCount;i++)
			{
				datafile.getline(line,512);
				Quaternion rotationi;
				Vector3f positioni;
				Vector3f vrotationi;
				int mno;
				int crowderyn;
				//cout << line << endl;
				char * token = strtok(line, ":");
				sscanf(token,"%d %d",&crowderyn,&mno); // new format
				
				if (crowderyn==1 && !printCrowders)
					continue; // ignore crowders
 
				rotationi.w = atof(strtok(NULL, " "));
				rotationi.x = atof(strtok(NULL, " "));
				rotationi.y = atof(strtok(NULL, " "));
				rotationi.z = atof(strtok(NULL, " "));
				positioni.x = atof(strtok(NULL, " "));
				positioni.y = atof(strtok(NULL, " "));
				positioni.z = atof(strtok(NULL, " "));

				rotationi.normalize();
				instance.molecules[i].setRotation(rotationi);

				if (crowderyn==0) // only two in the simulation
				{
					Vector3f v = positioni - instance.molecules[0].position;
					
					positioni.x -= boundary * round(v.x/boundary);
					positioni.y -= boundary * round(v.y/boundary);
					positioni.z -= boundary * round(v.z/boundary);
				}
				instance.molecules[i].setPosition(positioni);
				
				
			}
			memset(savename,0,256);
			instanceNo++;
			sprintf(savename, "pdbgen_results/sample_%d_%0.1fK_%5.2f.pdb",instanceNo,temperature,energy);
			//if (energy < -70.0)
			instance.saveAsSinglePDB(savename);
			//double drms(DRMS(&instance,&reference,reversedrms));
			//drmslog << energy << " " << drms << endl;
			cout << "Saved frame " << instanceNo << ": "  << savename << endl;


		}
	}
	datafile.close();
	//drmslog.close();
	cout << "Exiting..." << endl;
	return 0;
}
