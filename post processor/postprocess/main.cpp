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
#include "drms.h"
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

AminoAcids aminoAcidData;

// calculate the root mean square between replicas
/*double DRMS(Replica* sim, Replica* exp)
{
	double drms = 0.0;
	double N = 0.0;
	double native_contacts = 0;
	double experimental_contacts = 0;
	double binding_contacts = 0;


	for (int moleculeI=0; moleculeI<sim->moleculeCount;moleculeI++)
	{
		for (int moleculeJ=moleculeI+1; moleculeJ<sim->moleculeCount;moleculeJ++)
		{
			for (int residueI=0; residueI<sim->molecules[0].residueCount;residueI++)
			{
				for (int residueJ=0; residueJ<sim->molecules[1].residueCount;residueJ++)
				{
						N++;
						drms += abs(
								(sim->molecules[moleculeI].Residues[residueI].position - sim->molecules[moleculeJ].Residues[residueJ].position).magnitude() -
								(exp->molecules[moleculeI].Residues[residueI].position - exp->molecules[moleculeJ].Residues[residueJ].position).magnitude()
								);

						//if (residueI%50==0&&residueJ%50==0)
						//cout << residueI << ":" << residueJ << " " << (sim->molecules[0].Residues[residueI].position - sim->molecules[1].Residues[residueJ].position).magnitude() -
						//		(exp->molecules[0].Residues[residueI].position - exp->molecules[1].Residues[residueJ].position).magnitude() << endl;

						// get native contacts. etc for quality stuff.

						// if the experiment pair counts as a native contact
						//if ((exp->molecules[moleculeI].Residues[residueI].position - exp->molecules[moleculeJ].Residues[residueJ].position).magnitude() < 12*Angstrom)
						//{
						//	native_contacts++;
						//	if ((sim->molecules[moleculeI].Residues[residueI].position - sim->molecules[moleculeJ].Residues[residueJ].position).magnitude() < 15*Angstrom)
						//	{
						//		experimental_contacts++;
						//	}
						//}
						//if ((sim->molecules[moleculeI].Residues[residueI].position - sim->molecules[moleculeJ].Residues[residueJ].position).magnitude() < 12*Angstrom)
						//{
						//	binding_contacts++;
						//}

				}
			}
		}
	}
	drms /= N;
	//cout << "DRMS (angstrom): " << drms << endl;
	//cout << "Quality of interfacial residue data [0:1]: " << experimental_contacts/native_contacts << endl;
	//cout << "Binding contacts discovered: " << binding_contacts << endl;
	return drms;
}*/

bool extractMolecule(Molecule *molecule, ifstream &trajectory, int length, int index)
{
	bool success = false;
	char line[512] = {' '};
	char tmpLine[512] = {' '};
	int atomsRead = 0;
	molecule->Residues = new Residue[length];
	molecule->residueCount = length;
	while (atomsRead < length && trajectory.good())
	{
		trajectory.getline(line,512);
		//cout << "READ: '" << line << "'" <<  endl;
		if (line==NULL)
		{
			return false;
		}
		strncpy(tmpLine,line,strlen(line));
		char *token = strtok(tmpLine," ");

		if (token != NULL)
		{
			//discard MODEL, TITLE and TER elements, only read atoms
			if (strcmp(token,"ATOM")==0)
			{
				/*
				1-6 			"ATOM  "
				7 - 11         Integer         serial        Atom serial number.
				13 - 16        Atom            name          Atom name.
				17             Character       altLoc        Alternate location indicator.
				18 - 20        Residue name    resName       Residue name.
				22             Character       chainID       Chain identifier.
				23 - 26        Integer         resSeq        Residue sequence number.
				27             AChar           iCode         Code for insertion of residues.
				31 - 38        Real(8.3)       x             Orthogonal coordinates for X in
				39 - 46        Real(8.3)       y             Orthogonal coordinates for Y in
				47 - 54        Real(8.3)       z             Orthogonal coordinates for Z in
				55 - 60        Real(6.2)       occupancy     Occupancy.
				61 - 66        Real(6.2)       tempFactor    Temperature factor.
				73 - 76        LString(4)      segID         Segment identifier, left-justified.
				77 - 78        LString(2)      element       Element symbol, right-justified.
				79 - 80        LString(2)      charge        Charge on the atom.
				*/
				int serial(0);
				char name[5] = {' '};
				char altLoc = ' ';
				char resName[4] = {' '};
				char chainID = ' ';
				int resSeq(0);
				char iCode = ' ';
				float x(0);
				float y(0);
				float z(0);
				float occupancy(0);
				float tempFactor(0);
				char segID[5] = {' '};
				char element[3] = {' '};
				char charge[3] = {' '};

				char tmp[6] = {' '};

				// parse serial
				strncpy(tmp,line+6,4);
				sscanf (tmp,"%d",&serial);

				// parse name
				strncpy(tmp,line+12,4);
				sscanf (tmp,"%s",name);

				//parse altLoc, has to be this char if present
				altLoc = line[16];

				// parse resName
				strncpy(tmp,line+17,3);
				sscanf (tmp,"%s",resName);
				//cout << "resName" << resName << endl;

				strncpy(tmp,line+21,1);
				sscanf (tmp,"%c",&chainID);

				// parse resName
				strncpy(tmp,line+22,4);
				sscanf (tmp,"%d",&resSeq);

				strncpy(tmp,line+26,1);
				sscanf (tmp,"%c",&iCode);

				// parse x, y, z
				char tmpVals[36]  = {' '};

				//strncpy(tmpVals,line+30,35);
				//char *tmptokens = strtok(tmpVals," ");
				strncpy (tmpVals,line+30,8);
				//	sscanf (tmpVals,"%8.3f%8.3f%8.3f%6.2f%6.2f",&x,&y,&z,&occupancy,&tempFactor);
				sscanf (tmpVals,"%f",&x);
				//tmptokens = strtok(NULL," ");
				strncpy (tmpVals,line+38,8);
				sscanf (tmpVals,"%f",&y);
				strncpy (tmpVals,line+46,8);
				//tmptokens = strtok(NULL," ");
				sscanf (tmpVals,"%f",&z);
				strncpy (tmpVals,line+54,6);
				//tmptokens = strtok(NULL," ");
				sscanf (tmpVals,"%f",&occupancy);
				strncpy (tmpVals,line+60,6);
				//tmptokens = strtok(NULL," ");
				sscanf (tmpVals,"%f",&tempFactor);

				if (strlen(line)>=80)
				{
					// parse segID
					strncpy(tmp,line+72,4);
					sscanf (tmp,"%s",segID);
					// parse element
					strncpy(tmp,line+76,2);
					sscanf (tmp,"%s",element);
					// parse charge
					strncpy(tmp,line+78,2);
					sscanf (tmp,"%s",charge);
				}

				// push the atom onto the vector of our structure if it is a CA
				if (strcmp(name,"CA") == 0) // it is a CA atom. center of our bead
				{
					molecule->Residues[atomsRead].aminoAcidIndex = aminoAcidData.getAminoAcidIndex(resName);
					molecule->Residues[atomsRead].electrostaticCharge = float(aminoAcidData.get(molecule->Residues[atomsRead].aminoAcidIndex).electrostaticCharge);// * E_charge;
					molecule->Residues[atomsRead].vanderWaalRadius = aminoAcidData.data[molecule->Residues[atomsRead].aminoAcidIndex].vanderWaalRadius;
					if (molecule->Residues[atomsRead].vanderWaalRadius != aminoAcidData.get(molecule->Residues[atomsRead].aminoAcidIndex).vanderWaalRadius)
					{
						cout << molecule->Residues[atomsRead].aminoAcidIndex << " differing vdw radii:" << endl;
					}
					molecule->Residues[atomsRead].position = Vector3f(x,y,z);
					molecule->Residues[atomsRead].relativePosition = Vector3f(0,0,0);
					molecule->Residues[atomsRead].resSeq = atomsRead+1;

					molecule->Residues[atomsRead].moleculeID = index;
					atomsRead++;
				}


			}
			/*else if (strcmp(token,"ENDMDL")==0)
			{
				cout << "read " << atomsRead << " atoms (encountered ENDMDL)" << endl;
				return true;
			}*/
		} // null token
		else
		{
			cout << "Null string ends trajectory file." << endl;
			return false;
		}
	}// all atoms read in
	//cout << "read " << atomsRead << " atoms" << endl;
	if (atomsRead > 0)
		return true;
	else
		return false;

};

void writeToTrajectory(Replica *replica, ofstream &file, int frame)
{
	if (replica->moleculeCount == 0)
		return;
	char * output = new char[256];
	sprintf(output,"TITLE     frame t=%7.3f\n",float(frame));
	file << output;
	sprintf(output,"MODEL%9d\n",frame);
	file << output;
	sprintf(output,"REMARK potential,drms: %0.10f %0.10f \n",replica->potential, replica->drms);
	file << output;

	char chainId = 'A';
	int itemcount = 0;
	int lastSeqNo = 0;

	Molecule *molecules = replica->molecules;

	for (size_t m=0;m<replica->moleculeCount;m++)
	{
		size_t i=0;
		while (i<molecules[m].residueCount)
		{
			itemcount++;
			sprintf(output,"ATOM  %5d %4s%C%3s %C%4d%C  %8.3f%8.3f%8.3f%6.2f%6.2f\n",itemcount,"CA",' ',aminoAcidData.get(molecules[m].Residues[i].aminoAcidIndex).getSNAME(),chainId,molecules[m].Residues[i].resSeq,' ',molecules[m].Residues[i].position.x,molecules[m].Residues[i].position.y,molecules[m].Residues[i].position.z,1.0f,1.0f);
			file << output;

			lastSeqNo = molecules[m].Residues[i].resSeq;
			i++;
		}
		sprintf(output,"TER   %5d      %3s %C%4d \n",itemcount,aminoAcidData.get(molecules[m].Residues[i-1].aminoAcidIndex).getSNAME(),chainId,lastSeqNo);
		file << output;

		chainId++;
	}
	sprintf(output,"ENDMDL\n");
	file << output;

	file.flush();
};


int main(int argc, char **argv)
{
	cout <<	"Post process pdb trajectoies output by g_cluster." << endl;
	cout << "args: <clustertrajectory.pdb> processTraj.exp <molecules per frame> <list of molecule sizes, space delimited>" << endl;
	cout << "processTraj.exp contains a list of molecule files, 1 per line for drms" << endl;

	cout << "Loading amino acid data. " << AMINOACIDDATASOURCE;
	aminoAcidData.loadAminoAcidData(AMINOACIDDATASOURCE);
	cout << "Loading pair lookup table. " << LJPDSOURCE << endl;
	aminoAcidData.loadLJPotentialData(LJPDSOURCE);

	char *filename = new char[256];

	cout << "Opening: " << argv[1] << endl;
	ifstream trajectoryfile(argv[1]);
	if (!trajectoryfile.good())
	{
		cout << "Failed" << endl;
		exit(0);
	}

	ifstream referenceMolecules("processTraj.exp");
	if (!referenceMolecules.good())
	{
		cout << "Failed to open reference structure file: processTraj.exp" << endl;
		exit(0);
	}

	Replica reference;
	reference.aminoAcids = aminoAcidData;

	char line[512] = {0};
	while (referenceMolecules.good())
	{
		referenceMolecules.getline(line,512);
		cout << "Loading " << line << endl;
		if (strlen(line)>0)
			reference.loadMolecule(line);
	}
	referenceMolecules.close();


	int molecules;
	int *lengths;

	/*molecules = reference.moleculeCount;
	lengths = new int[molecules];
	for (int i=0;i<molecules;i++)
	{
		lengths[i] = reference.molecules[i].residueCount;
	}
	*/
	cout << "Number of Molecules per frame: ";

	molecules = atoi(argv[2]);
	cout << molecules << endl;
	if(molecules==0)
	{
		cout << "molecules empty" << endl;
		exit(-1);
	}	
	lengths = new int[molecules];

	for (int i=0;i<molecules;i++)
	{
		cout << "Length of molecule " << i+1 << ": ";
		lengths[i] = atoi(argv[3+i]);
		cout << lengths[i] << endl;
	}


	strcpy(filename,argv[1]);
	strcat(filename,".delim.pdb");
	ofstream completeTrajectoryFile(filename);
	strcpy(filename,argv[1]);
	strcat(filename,".drms.dat");
	ofstream drms_e_Plot(filename);
	drms_e_Plot.precision(6);


	// USE CUDA
	float * ljp_t;
	cudaMalloc((void**)&ljp_t,LJArraySize);
	cutilCheckMsg("Failed to cudaMalloc");
	copyLJPotentialDataToDevice(ljp_t,&aminoAcidData);
	bindLJTexture(ljp_t);

	int frame(1);
	bool success(true);

	while (success)
	{
		double drms = 0.0;
		Replica workingReplica;
		workingReplica.aminoAcids = aminoAcidData;
		workingReplica.molecules = new Molecule[molecules];
		workingReplica.moleculeCount = molecules;
		cout << "Processing frame " << frame << " " << endl;

		for (int m=0;m<molecules;m++)
		{
			// copy in each atom to a molecule
			success = extractMolecule(&workingReplica.molecules[m],trajectoryfile,lengths[m],m);
		}

		if (success)  // read a non 0 frame
		{
			// get DRMS
			drms = DRMS(&workingReplica,&reference,true);
			workingReplica.setDeviceLJPotentials(ljp_t);
			workingReplica.setBlockSize(64);
			workingReplica.ReplicaDataToDevice();
			workingReplica.potential = workingReplica.EonDevice();
			//workingReplica.potential = workingReplica.E();
			workingReplica.drms = drms;

			workingReplica.FreeDevice();

			drms_e_Plot << workingReplica.potential << "\t" << drms << endl;
			//cout << workingReplica.potential << " " << workingReplica.drms << endl;
			writeToTrajectory(&workingReplica,completeTrajectoryFile, frame);
			frame++;
		}
		else
		{
			cout << "Finished." << endl;
		}
	}
	drms_e_Plot.close();
	trajectoryfile.close();
	completeTrajectoryFile.close();


	unbindLJTexture();
	cudaFree(ljp_t);

}
