#ifndef AMINOACIDS_H_
#define AMINOACIDS_H_

#include "AminoAcid.h"

class AminoAcids
{
	public:  // we want public data, saves function call overhead.
		AminoAcids();
		AminoAcids(char* inputDataFileAminoAcids, char* inputDataFileLJPotentials);
		//~AminoAcids() {};
		AminoAcid data[20];
		float LJpotentials[20][20];
		bool loadAminoAcidData(const char* inputDataFile);
		bool loadLJPotentialData(const char* inputDataFile);
		bool clean;  // class contains data if true

		AminoAcid get(int aaIndex);
		AminoAcid get(char aminoAcidSymbol);
		AminoAcid get(const char* aminoAcidName);
		int getAminoAcidIndex (const char* name);
		int getAminoAcidIndex (char name);
		float getLJPairPotential(int indexA, int indexB);
		void printAcidData();
		void printPotentials();  // prints LJ lookup table
		float getLJp(int a1, int a2);
};

#endif /* AMINOACIDS_H_ */
