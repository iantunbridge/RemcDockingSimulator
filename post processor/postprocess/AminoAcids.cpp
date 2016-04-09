/*
 * AminoAcids.cpp
 *
 *  Created on: 2008/10/16
 *      Author: ian
 */

#include "AminoAcids.h"

AminoAcids::AminoAcids()
{
	clean = false;
}

AminoAcids::AminoAcids(char* inputDataFileAminoAcids, char* inputDataFileLJPotentials)
{
	clean = loadAminoAcidData(inputDataFileAminoAcids) && loadLJPotentialData(inputDataFileLJPotentials);
}

inline int min(size_t x, int y)
{
	return ((int(x)<y) ? int(x): y);
};

bool AminoAcids::loadAminoAcidData(const char* filename)
{
	ifstream input(filename);
	if (!input.good())
	{
		cout << "-!- Failed to open file: " << filename << "\n";
		exit(0);
		return false;
	}

	int i=0;
	char * line = new char[256];
	input.getline(line,255);
	cout << endl;
	while (!input.eof())
	{
		// ignore # lines
		if (strcmp(line,"")!=0 && strncmp(line,"#",1)!=0)
		{
			// read name
			memset(data[i].name,0,20);
			memset(data[i].sname,0,4);


			strncpy(data[i].name,line,min(strlen(line),19));
			// read short name char[4]
			input.getline(line,255);
			strncpy(data[i].sname,line,min(strlen(line),3));
			char buffer[64];
			input.getline(buffer,63);
			//sscanf(buffer,"%1c",&data[i].shortSymbol);
			data[i].shortSymbol = buffer[0];
			input.getline(buffer,63);
			//sscanf(buffer,"%f",&data[i].vanderWaalRadius);
			data[i].vanderWaalRadius = atof(buffer);
			input.getline(buffer,63);
			sscanf(buffer,"%f",&data[i].electrostaticCharge);
			data[i].index = i;
			//cout << data[i].name << " " << data[i].sname << " " << data[i].vanderWaalRadius << " " << data[i].electrostaticCharge << endl;
		    i++;
		}
		input.getline(line,255);
	}
	input.close();
	delete [] line;
	return true;
}

bool AminoAcids::loadLJPotentialData(const char* filename)
{
	//read file contents
	ifstream input(filename);
	if (!input.good())
	{
		cout << "-!- Failed to open file: " << filename << "\n";
		exit(0);
		return false;
	}

	/*
	 * map the label sequence to the stored data sequence, ie, read in " Cys  Phe  Leu...."
	 * and note that Cys is index 8 or something like that using getAcidIndex() so store labels[0] = 8
	 * then use this to fill in data[x][y] such that its in the same order as the aminoacid stored data
	 * so that lookups of values dont need an extra lookup step or hashmap which is slower than a direct
	 * query to an array possition
	 */

	int labels[20];
	bool labelline = true;
	int labelcount = 0;
	int startindex = 0;
	char line[256];
	input.getline(line,255);
	while (!input.eof())
	{
		// ignore # lines
		// if its not a blank line (=0) or a line beginiing with #
		if (strcmp(line,"")!=0 && strncmp(line,"#",1)!=0)
		{
			//cout << "line:" << line << endl;
			if (labelline)
			{
				char *label = strtok (line," ");
				while (label != NULL)
				{
					labels[labelcount] = getAminoAcidIndex(label);
					label = strtok (NULL, " ");
					labelcount++;
				}
				labelline = false;
			}
			else
			{
				// ignore the first token as it is the label, same order as first row
				char *value = strtok (line," ");

				for (int j=startindex;j<20;j++)
				{
					value = strtok (NULL, " ");
					float e;// = atof(value);
					sscanf(value,"%f",&e);
					//printf ("%0.2f\n",e);

					LJpotentials[labels[startindex]][labels[j]] = e;
					LJpotentials[labels[j]][labels[startindex]] = LJpotentials[labels[startindex]][labels[j]];
				}
				startindex++;
				value = strtok (NULL, " ");
			}
		}
		input.getline(line,255);
	}
	input.close();
	return true;
}

float AminoAcids::getLJp(int a, int b)
{
	return LJpotentials[a][b];
}

int AminoAcids::getAminoAcidIndex(const char* NAME)
{
	if (strlen(NAME)<1)
		return -1;
	int index=20; // set to max index allowed
	char name[4] = {0};
	name[0] = toupper(NAME[0]);
	// search short name
	if (strlen(NAME)==3)
	{
		name[0] = toupper(NAME[0]);
		name[1] = tolower(NAME[1]);
		name[2] = tolower(NAME[2]);
		while (--index > -1)
		{
			if (strcmp(name,data[index].sname)==0)
				return index;
		}
	}
	else if (strlen(NAME)==1)
	{
		while (--index > -1)
		{
			if (name[0]==data[index].shortSymbol)
				return index;
		}
	}
	else
	{
		while (--index > -1)
		{
			if (strcmp(name,data[index].name)==0)
				return index;
		}
	}
	return index;
}

int AminoAcids::getAminoAcidIndex(char name)
{
	int index=20; // set to max index allowed
	// search short name
	while (index > -1)
	{
		index--;
		if (name==data[index].shortSymbol)
			return index;
	}
	return index;
}

void AminoAcids::printAcidData()
{
	int i=0; // set to max index allowed
	// search short name
	while (i < 20)
	{
		cout << data[i].sname << " " << data[i].shortSymbol << " " << data[i].electrostaticCharge << " " << data[i].vanderWaalRadius << " " << endl;
		i++;
	}
}

AminoAcid AminoAcids::get(int index)
{
	return data[index];
}

AminoAcid AminoAcids::get(char aminoAcidSymbol)
{
	return get(getAminoAcidIndex(aminoAcidSymbol));
}

AminoAcid AminoAcids::get(const char* aminoAcidName)
{
	return get(getAminoAcidIndex(aminoAcidName));
}

void  AminoAcids::printPotentials()
{
	cout << "------------------------------------------- LJ look up table -----------------------------------------" << endl;
	cout << "LJx  ";
	for (int j = 0; j<20; j++)
	{
		cout << data[j].sname << "   ";
	}
	cout << endl;
	for (int i = 0; i<20; i++)
	{
		cout << data[i].sname << "  ";
		for (int j = 0; j<20; j++)
		{
			printf ("%1.2f ",LJpotentials[i][j]);
		}
		cout << endl;
	}
	cout << "---------------------------------------- End: LJ look up table ---------------------------------------" << endl;;

	return;
}
