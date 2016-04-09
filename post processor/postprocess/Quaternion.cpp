#include "Quaternion.h"


Quaternion::Quaternion(double w,double x,double y,double z)
{
	this->w = w;
	this->x = x;
	this->y = y;
	this->z = z;
}

double Quaternion::getMagnitude()
{
	magnitude = sqrt(w*w + x*x + y*y + z*z);
	return magnitude;
}

// normalize the quaternion: corrects numerical inaccuracy in rotations
Quaternion Quaternion::normalize()
{
	getMagnitude();
	w /= magnitude;
	x /= magnitude;
	y /= magnitude;
	z /= magnitude;
	return Quaternion(w,x,y,z);
}

Quaternion Quaternion::getNormalized()
{
	getMagnitude();
	return Quaternion(w/magnitude,x/magnitude,y/magnitude,z/magnitude);
}

Quaternion Quaternion::multiply (const Quaternion q2)
{
	//(Q1 * Q2).w = (w1w2 - x1x2 - y1y2 - z1z2)
	//(Q1 * Q2).x = (w1x2 + x1w2 + y1z2 - z1y2)
	//(Q1 * Q2).y = (w1y2 - x1z2 + y1w2 + z1x2)
	//(Q1 * Q2).z = (w1z2 + x1y2 - y1x2 + z1w2)
	Quaternion output;
	output.w = this->w*q2.w - this->x*q2.x - this->y*q2.y - this->z*q2.z;
	output.x = this->w*q2.x + this->x*q2.w + this->y*q2.z - this->z*q2.y;
	output.y = this->w*q2.y - this->x*q2.z + this->y*q2.w + this->z*q2.x;
	output.z = this->w*q2.z + this->x*q2.y - this->y*q2.x + this->z*q2.w;
	return output;
}

// quaternion multiplication
Quaternion Quaternion::operator * (const Quaternion q2)
{
	return multiply(q2);
}

inline Quaternion Quaternion::equals (const Quaternion q2)
{
	this->w = q2.w;
	this->x = q2.x;
	this->y = q2.y;
	this->z = q2.z;
	return *this;
}

Quaternion Quaternion::operator = (const Quaternion q2)
{
	return equals(q2);
}


// rotate a vector using a quaternion
Quaternion Quaternion::rotate (const double theta, const double x, const double y, const double z)
{
		double t = theta * 0.5;
		double temp = sin(t); // optimization (save 2 operations)
		return Quaternion (cos(t),x*temp, y*temp,z*temp).multiply(*this);
}

// rotate a vector using a quaternion
Quaternion Quaternion::rotate (const double theta, const Vector3f v)
{
	return rotate(theta,v.x,v.y,v.z);
}

Vector3f Quaternion::rotateVector(Vector3f v)
{
		//try to reduce the number of operations using temps (14 ops saved vs more tmp memory)
		double w2 = w*w; // save 2 a
		double x2 = x*x; //2 b
		double y2 = y*y; //2 c
		double z2 = z*z; //2 d
		double wx = w*x;
		double wy = w*y;
		double wz = w*z;
		double xy = x*y;
		double xz = x*z;
		double yz = y*z;

		return Vector3f(v.x*(w2+x2-y2-z2) + ( v.y*(xy-wz)   + v.z*(wy+xz) )*2.0,
						( v.x*(wz+xy)   + v.z*(yz-wx) )*2.0 + v.y*(w2-x2+y2-z2),
						( v.x*(xz-wy)   + v.y*(yz+wx) )*2.0 + v.z*(w2-x2-y2+z2));
}

RotationalMatrix creatRotationalMatrix(Quaternion q)
{
    RotationalMatrix r;

    double w2 = q.w*q.w; //2 ops
	double x2 = q.x*q.x; //2 ops
	double y2 = q.y*q.y; //2 ops
	double z2 = q.z*q.z; //2 ops
	double wx = q.w*q.x;
	double wy = q.w*q.y;
	double wz = q.w*q.z;
	double xy = q.x*q.y;
	double xz = q.x*q.z;
	double yz = q.y*q.z;

    r.r[0][0] = (w2+x2-y2-z2);
	r.r[0][1] = (xy-wz)*2.0;
	r.r[0][2] = (wy+xz)*2.0;
	r.r[1][0] = (wz+xy)*2.0;
	r.r[1][1] = (w2-x2+y2-z2);
	r.r[1][2] = (yz-wx)*2.0;
	r.r[2][0] = (xz-wy)*2.0;
	r.r[2][1] = (yz+wx)*2.0;
	r.r[2][2] = (w2-x2-y2+z2);
	return r;
}

Quaternion createRotationQuaternion(double theta, Vector3f v)
{
	///[w, x, y, z] = [cos(a/2), sin(a/2) * nx, sin(a/2)* ny, sin(a/2) * nz]
	theta *= 0.5;
	double temp = sin(theta); // optimisation (save 2 trig operations)
	return Quaternion(cos(theta),v.x*temp, v.y*temp, v.z*temp);
}





