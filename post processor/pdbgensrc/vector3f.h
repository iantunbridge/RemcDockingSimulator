#ifndef vector3f_h
#define vector3f_h

using namespace std;

#include <cmath>
#include <iostream>

#define max_err 0.0000000000001

// a general purpose 3d vector of floats a counter variable and a pointer to the same
class Vector3f
{
	public:
		float x;
		float y;
		float z;
		Vector3f() : x(0.0f), y(0.0f), z(0.0f) { };
		Vector3f(float xi, float yi, float zi) { x=xi;y=yi;z=zi; };
		Vector3f operator = (Vector3f v) { x=v.x; y=v.y; z=v.z; return Vector3f(x,y,z); };
		Vector3f operator += (Vector3f v);
		Vector3f set(Vector3f v) { x = v.x; y = v.y; z = v.z; return Vector3f(x,y,z);  };
		Vector3f set(float xi, float yi, float zi) { x = xi; y = yi; z = zi; return Vector3f(x,y,z); };
		Vector3f operator * (const float s);
		Vector3f operator / (const float s);
		Vector3f operator + (const Vector3f v);
		Vector3f operator - (const Vector3f v);
		Vector3f project(const Vector3f destination);
		float dot(const Vector3f v);
		Vector3f dotV(Vector3f v);
		Vector3f dotV(float xi,float yi,float zi) { return dotV(Vector3f(xi,yi,zi));};
		Vector3f cross(const Vector3f v);
		float angle(const Vector3f v);
		Vector3f normalize();
		float magnitude() { return sqrt(x*x+y*y+z*z); };
		void zero() { x=y=z=0; };
		float sumSquares();
		//float distance(const Vector3f &v);
		inline float distance(const Vector3f &b)
		{
			return sqrtf((x-b.x)*(x-b.x)+(y-b.y)*(y-b.y)+(z-b.z)*(z-b.z));
		}
		const char * c_str();
};

// a 3 dimentional vector of integers
class Vector3d
{
	public:
		long x;
		long y;
		long z;
		Vector3d()	{ x=y=z=0; };
		Vector3d(long xi, long yi, long zi) { x=xi;y=yi;z=zi; };

};

class Vector3double
{
	public:
		double x;
		double y;
		double z;
		Vector3double()  : x(0.0), y(0.0), z(0.0) 	{ };
		Vector3double(double xi, double yi, double zi) { x=xi;y=yi;z=zi; };
		double magnitude() { return sqrt(x*x+y*y+z*z); };
		Vector3double operator += (const Vector3f v);
		void normalizeInPlace();
		Vector3double operator = (Vector3double v) { x=v.x; y=v.y; z=v.z; return Vector3double(x,y,z); };
		const char * c_str();

};


class RotationalMatrix
{
	public:
		float r[3][3];
		Vector3f rotate(Vector3f v);
		RotationalMatrix ();
		RotationalMatrix (Vector3f u, float alpha); // creates the rotational matrix about the arbitrary axis u for angle alpha
		RotationalMatrix set(Vector3f u, float alpha); // creates the rotational matrix about the arbitrary axis u for angle alpha
		RotationalMatrix operator = (RotationalMatrix rm);

};

Vector3f operator * (float s,Vector3f v);
Vector3f operator * (double s,Vector3f v);
Vector3f operator / (float s,Vector3f v);
Vector3f operator / (double s,Vector3f v);
Vector3f getCrossProduct(float ax, float ay, float az, float bx, float by, float bz);
Vector3f getCrossProduct(const Vector3f a, const Vector3f b);
float getDotProduct(float ax, float ay, float az, float bx, float by, float bz);
float getDotProduct(const Vector3f a, const Vector3f b);
Vector3f normalize(const Vector3f v);
Vector3f * pnormalize(Vector3f* v);
inline float scalar_distance(const Vector3f &a, const Vector3f &b)
{
	return sqrtf((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z));
}
#endif
