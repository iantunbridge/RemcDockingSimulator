#ifndef MOLECULE_H_
#define MOLECULE_H_

#include <fstream>
#include <cstring>
#include <vector>
#include "definitions.h"
#include "AminoAcids.h"
#include "vector3f.h"
#include "Quaternion.h"
//#include "TorsionalLookupMatrix.h"
#include "Residue.h"
//#include "Link.h"


using namespace std;

struct charStruct
{
	char * data;
	const char * c_str() { return data; };
};

class Molecule
{
	public:

	Molecule();
	~Molecule();
	Molecule(const Molecule& m);
	Molecule operator = (const Molecule& m) { return Molecule(m); };
	void deleteResidues();
	void init(AminoAcid residues[], char *description, char *data);
	bool initFromPDB(const char* pdbfilename);
	void copy(const Molecule& m);
	void saveBeforeStateChange(const Molecule* m);
	void undoStateChange(const Molecule* m);
	void reserveResidueSpace(int size);
    void shallowCopy(Molecule & r);
    bool update(AminoAcid AminoAcids[]);
    bool checkCollisions(const Molecule & m);
    void updatePositions();
    char *print();
    void saveAsPDB(const char *filename);
    Vector3f calculateCenter();
    void setPosition(Vector3f v);
    bool translate(Vector3f v);
    void setRotation(Quaternion q);
    bool rotate(const Vector3double Raxis, const double angle);
    bool rotateQ(const Vector3double Raxis, const double angle);
    bool mutate(const int type, const int position, const float amount);
    float Ebond();
    float Eangle();
    float Etorsion();
    float Emembrane();
    Residue *getResidues();
    //Link *getLinks();
    Vector3f getPosition();
    int getIndex();
    Quaternion getRotation();
    char *getDescription();
    char *getDatadescription();
    int getLength();
    float getTemperature();
    short getLabel();
    double getPotential();
    bool checkSelfCollisions();
    bool checkCollisions(Molecule *molecule2);
    bool update();
    float getMoleculeRoleIdentifier();
    void setMoleculeRoleIdentifier(float moleculeRoleIdentifier);

    float getVolume();
    float calculateVolume();
    float calculateSASA();

    AminoAcids AminoAcidsData;
    //    TorsionalLookupMatrix torsions;
    Residue *Residues;
    int residueCount;
    //    Link *Links;
    size_t linkCount;
    size_t chainCount;
    Vector3f position;
    Quaternion rotation;
    Vector3f center;
    int index;
    float volume;
	//    int length;
	//    float temperature;
	//    short label;
    float moleculeRoleIdentifier;
    float translationalStep;
    float rotationalStep;
    char *filename;
    bool hasFilename;
    Vector3f xAxis;
    Vector3f yAxis;
    Vector3f zAxis;
    Vector3f translation;
    /*charStruct HEADER;
    charStruct TITLE;
    vector<charStruct> COMPND;
    vector<charStruct> SOURCE;
    vector<charStruct> KEYWDS;
    vector<charStruct> EXPDTA;
    vector<charStruct> AUTHOR;
    vector<charStruct> REVDAT;
    vector<charStruct> REMARK;
    vector<charStruct> SEQRES;
    vector<charStruct> HELIX;
    vector<charStruct> SHEET;
    vector<charStruct> TURN;
    vector<charStruct> CRYST1;
	vector<charStruct> ORIGX1;
	vector<charStruct> ORIGX2;
	vector<charStruct> ORIGX3;
	vector<charStruct> SCALE1;
	vector<charStruct> SCALE2;
	vector<charStruct> SCALE3;
	vector<charStruct> MASTER;*/

	private:
		float EchargedMenbrane(float distanceToMembrane,float temperature);

};

#endif /*MOLECULE_H_*/



