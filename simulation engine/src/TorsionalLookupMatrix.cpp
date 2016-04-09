#include "TorsionalLookupMatrix.h"

TorsionalLookupMatrix::TorsionalLookupMatrix(AminoAcids *a)
{
	//setAminoAcidLookupData(a);
}

TorsionalLookupMatrix::TorsionalLookupMatrix()
{

}

TorsionalLookupMatrix::~TorsionalLookupMatrix()
{
	//aminoAcids = AminoAcids();
}

void TorsionalLookupMatrix::setAminoAcidLookupData(AminoAcids &a)
{
	aminoAcids = a;
}

double TorsionalLookupMatrix::getV(int ai, int aj,int n)
{
	return V[ai][aj][n];
}

double TorsionalLookupMatrix::getSigma(int ai, int aj,int n)
{
	return Sigma[ai][aj][n];
}

bool TorsionalLookupMatrix::loadData(char *filename)
{
	//read file contents
	ifstream input(filename);
	if (!input.good())
	{
		cout << "Failed to open file.\n";
		return false;
	}

	char line[256];
	char acid0[5] = {0};
	char acid1[5] = {0};
	//memset(acid0,0,5);
	//memset(acid1,0,5);

	double v;
	double sigma;
	int  n;
	input.getline(line,255);

	while (input.good())
	{

		// ignore # lines
		// if its not a blank line (=0) or a line beginiing with #
		if (strcmp(line,"")!=0 && strncmp(line,"#",1)!=0)
		{
			sscanf(line,"%3s %3s %lf %d %lf",acid0,acid1,&v,&n,&sigma);
			int a0 = aminoAcids.getAminoAcidIndex(acid0);
			int a1 = aminoAcids.getAminoAcidIndex(acid1);
			V[a0][a1][n] = v;
			Sigma[a0][a1][n] = sigma/57.29577951308232; // convert to radians, i assume they are in degrees as the numbers are big ?
		}
		input.getline(line,255);

	}
	input.close();
	return true;
}
