#include "Link.h"

Link::Link()
{
	terminal = false;
}

Link::~Link()
{
}

Link::Link(const Link & l)
{
	Angle = l.Angle;
	BondLength = l.BondLength;
	TorsionAngle = l.TorsionAngle;
	flexible = l.flexible;
	terminal = l.terminal;
}

Link Link::operator =(Link l)
{
	return Link(l);
}
