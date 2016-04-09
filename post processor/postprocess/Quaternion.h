#ifndef QUATERNION_H_
#define QUATERNION_H_

using namespace std;

#include <cmath>

#define INCLUDE_OPTIMIZATIONS

#include "vector3f.h"



class Quaternion
{
	public:
		Quaternion()  : w(1.0), x(0.0), y(0.0), z(0.0), magnitude(1.0) {};  // default constuctor (1,0,0,0)
		Quaternion(double,double,double,double);  //explicit constructor
		double getMagnitude(); //get the magnitude of this
		Quaternion normalize(); //normalize the quarternion in place
		Quaternion getNormalized(); // get a normalised version without modifying the original
		Quaternion multiply (const Quaternion q2); //multiply current by another
		inline Quaternion equals (const Quaternion q2); //multiply current by another
		Quaternion operator * (const Quaternion);
		Quaternion operator = (const Quaternion); //overload previous operation
		Quaternion rotate(const Quaternion); //
		Quaternion rotate(const double theta, const Vector3f u);
		Quaternion rotate(const double theta, const double ux, const double uy, const double uz);
		void set(double wi,double xi,double yi,double zi) { w=wi;x=xi;y=yi;z=zi; };
		Vector3f rotateVector(Vector3f);
		inline double getW() { return w; };
		inline double getX() { return x; };
		inline double getY() { return y; };
		inline double getZ() { return z; };
		double w;
		double x;
		double y;
		double z;
		double magnitude; // store for speed reasons
};


class Quaternionf
{
	public:
	float w;
	float x;
	float y;
	float z;
	Quaternionf()  : w(1.0f), x(0.0f), y(0.0f), z(0.0f) {};  // default constuctor (1,0,0,0)
	Quaternion makeQuaternion() { return Quaternion(w,x,y,z); }
};

Quaternion createRotationQuaternion(double theta, Vector3f v);
RotationalMatrix creatRotationalMatrix(Quaternion q);

#endif /*QUATERNION_H_*/
