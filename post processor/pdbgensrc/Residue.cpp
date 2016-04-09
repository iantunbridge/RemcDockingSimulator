#include "Residue.h"

Residue::Residue()
{
	aminoAcidIndex = int(PADDER_IDENTIFIER);
}

Residue::~Residue()
{
}

Residue::Residue(const Residue & r)
{
	aminoAcidIndex = r.aminoAcidIndex;
	electrostaticCharge = r.electrostaticCharge;
	vanderWaalRadius = r.vanderWaalRadius;
	position = r.position;
	relativePosition = r.relativePosition;
	//sasa = r.sasa;
	isCrowder = r.isCrowder;
	chainId = r.chainId;
	resSeq = r.resSeq;

}
