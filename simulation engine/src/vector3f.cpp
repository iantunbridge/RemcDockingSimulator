#include "vector3f.h"

Vector3f getCrossProduct(const float ax, const float ay, const float az, const float bx, const float by, const float bz)
{
	return Vector3f (ay*bz - az*by,az*bx - ax*bz,ax*by - ay*bx);
};

Vector3f getCrossProduct(const Vector3f a, const Vector3f b)
{
	return getCrossProduct(a.x,a.y,a.z,b.x,b.y,b.z);
};

float getDotProduct(const float ax, const float ay, const float az, const float bx, const float by, const float bz)
{
	return ax*bx + ay*by + az*bz;
};

float getDotProduct(const Vector3f a, const Vector3f b)
{
	return getDotProduct(a.x,a.y,a.z,b.x,b.y,b.z);
};

Vector3f normalize(const Vector3f v)
{
	Vector3f n(v.x,v.y,v.z);
	float magnitude = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
	if (magnitude != 0)
	{
		n.x /= magnitude;
		n.y /= magnitude;
		n.z /= magnitude;
	}
	return n;
};

void Vector3double::normalizeInPlace()
{
	double magnitude = sqrt(x*x + y*y + z*z);

	x /= magnitude;
	y /= magnitude;
	z /= magnitude;
};

Vector3f * pnormalize(Vector3f *v)
{
	Vector3f *n = new Vector3f();
	float magnitude = sqrt(v->x*v->x + v->y*v->y + v->z*v->z);
	if (magnitude != 0)
	{
		n->x = v->x/magnitude;
		n->y = v->y/magnitude;
		n->z = v->z/magnitude;
	}
	return n;
};

Vector3f Vector3f::normalize()
{
	float magnitude = sqrt(x*x + y*y + z*z);
	return magnitude != 0 ? Vector3f(x,y,z)/magnitude : Vector3f(x,y,z);
};

Vector3f operator * (const float s, const Vector3f v)
{
	return Vector3f(v.x*s, v.y*s, v.z*s);
};

Vector3f operator * (const double s, const Vector3f v)
{
	return Vector3f(v.x*s, v.y*s, v.z*s);
};

Vector3f operator / (const float s, const Vector3f v)
{
	return Vector3f(s/v.x, s/v.y, s/v.z);
};

Vector3f operator / (const double s, const Vector3f v)
{
	return Vector3f(s/v.x, s/v.y, s/v.z);
};

Vector3f Vector3f::operator * (const float s)
{
	return Vector3f(x*s, y*s, z*s);
};

Vector3f Vector3f::operator / (const float s)
{
	return Vector3f(x/s, y/s, z/s);
};

Vector3f Vector3f::operator + (const Vector3f v)
{
	return Vector3f(x+v.x, y+v.y, z+v.z);
};

Vector3f Vector3f::operator += (const Vector3f v)
{
	this->x = this->x + v.x;
	this->y = this->y + v.y;
	this->z = this->z + v.z;
};

Vector3f Vector3f::operator - (const Vector3f v)
{
	return Vector3f(x-v.x, y-v.y, z-v.z);
};

Vector3f Vector3f::project(Vector3f a)
{
	double amag = a.magnitude();
	return amag != 0 ? Vector3f(getDotProduct(a,*this)/(amag*amag)*a) : Vector3f(0,0,0);
	//return Vector3f(getDotProduct(a,*this)/(amag*amag)*a);
};

float Vector3f::angle(Vector3f v)
{
	double cmag = v.magnitude()*magnitude();
	if (cmag<max_err*max_err)
		return 0;//acos(0);
	return acos(dot(v)/cmag);
}

float Vector3f::dot(Vector3f v)
{
	return x*v.x + y*v.y + z*v.z;
}

Vector3f Vector3f::dotV(Vector3f v)
{
	return Vector3f(x*v.x,y*v.y,z*v.z);
}

Vector3f Vector3f::cross(Vector3f v)
{
	return Vector3f (y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x);
}

float Vector3f::sumSquares()
{
	return x*x+y*y+z*z;
}

const char * Vector3f::c_str()
{
//	printf("(%1.3e;%1.3e;%1.3e)",x,y,z);
	return "";
}

const char * Vector3double::c_str()
{
//	printf("(%1.3e;%1.3e;%1.3e)",x,y,z);
	return "";
}

RotationalMatrix::RotationalMatrix (Vector3f u, float alpha)
{
    set (u,alpha);
}

RotationalMatrix::RotationalMatrix ()
{

}

RotationalMatrix RotationalMatrix::set (Vector3f u, float alpha)
{
    float c = cos(alpha);
    float s = sin(alpha);
	r[0][0] = u.x * u.x + (1-u.x * u.x)*c;
	r[0][1] = u.x * u.y * (1-c) - u.z*s;
	r[0][2] = u.x * u.z * (1-c) + u.y*s;
	r[1][0] = u.x * u.y * (1-c) + u.z*s;
	r[1][1] = u.y * u.y + (1-u.y * u.y)*c;
	r[1][2] = u.y * u.z * (1-c) - u.x*s;
	r[2][0] = u.x * u.z * (1-c) - u.y*s;
	r[2][1] = u.y * u.z * (1-c) + u.x*s;
	r[2][2] = u.z * u.z + (1-u.z * u.z)*c;
    return *this;
};

RotationalMatrix RotationalMatrix::operator = (RotationalMatrix rm)
{
    r[0][0] = rm.r[0][0];
    r[0][1] = rm.r[0][1];
    r[0][2] = rm.r[0][2];
    r[1][0] = rm.r[1][0];
    r[1][1] = rm.r[1][1];
    r[1][2] = rm.r[1][2];
    r[2][0] = rm.r[2][0];
    r[2][1] = rm.r[2][1];
    r[2][2] = rm.r[2][2];
    return *this;
}

Vector3f RotationalMatrix::rotate (Vector3f v)
{
	return Vector3f (r[0][0]*v.x + r[0][1]*v.y + r[0][2]*v.z,
	                 r[1][0]*v.x + r[1][1]*v.y + r[1][2]*v.z,
	                 r[2][0]*v.x + r[2][1]*v.y + r[2][2]*v.z);

};

/*float Vector3f::distance(const Vector3f &b)
{
	return sqrtf((x-b.x)*(x-b.x)+(y-b.y)*(y-b.y)+(z-b.z)*(z-b.z));
}*/



