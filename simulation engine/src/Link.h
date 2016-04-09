#ifndef LINK_H_
#define LINK_H_

class Link
{
	public:
		Link();
		virtual ~Link();
		Link(const Link& l);
		Link operator = (Link l);

	//private:
		float Angle;  // angular displacement relative to previous bead
		float TorsionAngle; // twist relative to previous bead
		float BondLength; // distance between this and the previous bead
		float flexible; // 0 (inflexible) or 1 (flexible): determines fixed tertiary structures
		bool terminal; // true if a new chain starts from this link;
};

#endif /*LINK_H_*/
