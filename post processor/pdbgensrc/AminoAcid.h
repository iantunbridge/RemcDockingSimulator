#ifndef AMINOACID_H_
#define AMINOACID_H_

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>

using namespace std;

class AminoAcid
{
public:
	AminoAcid();
	virtual ~AminoAcid();
	char name[20];
	char sname[4];
	char shortSymbol;
	float vanderWaalRadius;
	float electrostaticCharge;
	int index;
	char *getSNAME();
	char *getName();
	char getShortSymbol() { return shortSymbol;} ;
};
/*
AminoAcid *AminoAcidData;

bool loadAminoAcidData(AminoAcid data[], char* filename);
int getAcidIndex(AminoAcid data[], char* name);
int getAcidIndex(char* name);
int getAcidIndex(AminoAcid data[], char name);
void printAcidData(AminoAcid data[]);
*/
#endif /*AMINOACID_H_*/
