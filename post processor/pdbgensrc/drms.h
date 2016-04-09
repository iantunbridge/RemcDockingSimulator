/*
 * drms.h
 *
 *  Created on: 07 Jul 2010
 *      Author: ian
 */

#ifndef DRMS_H_
#define DRMS_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <cmath>
#include <sys/mman.h>
/*#include <gsl/gsl_qrng.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>*/
#include "definitions.h"
#include "AminoAcid.h"
#include "Replica.h"
#include "vector3f.h"
#include "Quaternion.h"

double DRMS(Replica* simulation, Replica* native, bool reverseMoleculeOrder);

#endif /* DRMS_H_ */
