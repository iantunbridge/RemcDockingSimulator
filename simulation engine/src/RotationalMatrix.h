#ifndef ROTMATRIX
#define ROTMATRIX
class RotationalMatrixF
{
	public:
		float r[3][3];
		Vector3f rotate(Vector3f v);
		RotationalMatrixF (Vector3f u, float alpha); // creates the rotational matrix about the arbitrary axis u for angle alpha
		
};

#endif
