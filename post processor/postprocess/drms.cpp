
#include "drms.h"

using namespace std;

// calculate the root mean square between replicas
double DRMS(Replica* sim, Replica* exp, bool reverseOrder)
{
	double drms = 0;
	double N = 0;
	double native_contacts = 0;
	double experimental_contacts = 0;
	double binding_contacts = 0;

	for (int moleculeI = 0; moleculeI < sim->moleculeCount; moleculeI++)
	{
		for (int moleculeJ = moleculeI+1; moleculeJ < sim->moleculeCount; moleculeJ++)
		{
			for (int residueI = 0; residueI	< sim->molecules[moleculeI].residueCount; residueI++)
			{
				for (int residueJ = 0; residueJ < sim->molecules[moleculeJ].residueCount; residueJ++)
				{
					N++;
					drms
					+= abs((sim->molecules[moleculeI].Residues[residueI].position - sim->molecules[moleculeJ].Residues[residueJ].position).magnitude() -
						   (exp->molecules[moleculeI].Residues[residueI].position - exp->molecules[moleculeJ].Residues[residueJ].position).magnitude());

				}
			}
		}
	}
	drms /= N;

	if (reverseOrder)
	{
		double drms2(0.0);
		for (int moleculeI = 0; moleculeI < sim->moleculeCount; moleculeI++)
		{
			for (int moleculeJ = moleculeI+1; moleculeJ < sim->moleculeCount; moleculeJ++)
			{
				for (int residueI = 0; residueI	< sim->molecules[moleculeI].residueCount; residueI++)
				{
					for (int residueJ = 0; residueJ< sim->molecules[moleculeJ].residueCount; residueJ++)
					{
						drms2 += abs( (sim->molecules[moleculeI].Residues[residueI].position - sim->molecules[moleculeJ].Residues[residueJ].position).magnitude() -
								      (exp->molecules[moleculeJ].Residues[residueI].position - exp->molecules[moleculeI].Residues[residueJ].position).magnitude());

					}
				}
			}
		}
		drms = min(drms,drms2/N);
	}
	//cout << "DRMS (angstrom): " << drms/Angstrom << endl;
	//cout << "Quality of interfacial residue data [0:1]: " << experimental_contacts/native_contacts << endl;
	//cout << "Binding contacts discovered: " << binding_contacts << endl;
	return drms;
}
