#include "Camera.h"

// NOTE: this code was originally written using code from nehe.gamedev.net tutorials, but a long time before this tut

Camera::Camera()
{
	position.x = 0.0f;
	position.y = 0.0f;
	position.z = -10.0f;
	view.x = 0.0f;
	view.y = 0.0f;
	view.z = 1.0f;
	up.x = 0.0f;
	up.y = 1.0f;
	up.z = 0.0f;
	screenwidth = 1024;
	screenheight = 768;
};

void Camera::setPosition(float xposition, float yposition, float zposition,  float xview, float yview, float zview, float upx, float upy, float upz)
{
	position.x = xposition; position.y = yposition; position.z = zposition;
	view.x = xview; view.y = yview; view.z = zview;
	up.x = upx; up.y = upy; up.z = upz;
};

void Camera::setPosition(Vector3f pos, Vector3f view, Vector3f up)
{
	setPosition (pos.x,pos.y,pos.z,view.x,view.y,view.z,up.x,up.y,up.z);
}

// rotates the camera about a point
void Camera::rotate(float x, float y, float z)
{
	Vector3f dir;

	dir.x = view.x - position.x;
	dir.y = view.y - position.y;
	dir.z = view.z - position.z;

	if(x) // if non zero
	{
		view.z = (float)(position.z + sin(x)*dir.y + cos(x)*dir.z);
		view.y = (float)(position.y + cos(x)*dir.y - sin(x)*dir.z);
	}
	if(y)
	{
		view.z = (float)(position.z + sin(y)*dir.x + cos(y)*dir.z);
		view.x = (float)(position.x + cos(y)*dir.x - sin(y)*dir.z);
	}
	if(z)
	{
		view.x = (float)(position.x + sin(y)*dir.y + cos(z)*dir.x);
		view.y = (float)(position.y + cos(y)*dir.y - sin(z)*dir.x);
	}
}

// enables camera rotation about a point

void Camera::pointRotate(Vector3f point, float x, float y, float z)
{
	Vector3f vector;
	vector.x = position.x - point.x;
	vector.y = position.y - point.y;
	vector.z = position.z - point.z;

	if(x)
	{
		position.z = (float)(point.z + sin(x)*vector.y + cos(x)*vector.z);
		position.y = (float)(point.y + cos(x)*vector.y - sin(x)*vector.z);
	}
	if(y)
	{
		position.z = (float)(point.z + sin(y)*vector.x + cos(y)*vector.z);
		position.x = (float)(point.x + cos(y)*vector.x - sin(y)*vector.z);
	}
	if(z)
	{
		position.x = (float)(point.x + sin(z)*vector.y + cos(z)*vector.x);
		position.y = (float)(point.y + cos(z)*vector.y - sin(z)*vector.x);
	}

};

//moves the camera
void Camera::move(float speed)
{
	Vector3f vector;

	vector.x = view.x - position.x;
	vector.y = view.y - position.y;
	vector.z = view.z - position.z;

	position.x += vector.x * speed;
	position.y += vector.y * speed;
	position.z += vector.z * speed;
	view.x += vector.x * speed;
	view.y += vector.y * speed;
	view.z += vector.z * speed;
};


// stafes the camera
void Camera::side(float speed)
{
	Vector3f vector = getCrossProduct(view.x-position.x,view.y-position.y,view.z-position.z,up.x,up.y,up.z);
	float norm = (float)sqrt(vector.x*vector.x+vector.y*vector.y+vector.z*vector.z);

	position.x += vector.x * speed / norm;
	position.y += vector.y * speed / norm;
	position.z += vector.z * speed / norm;
	view.x += vector.x * speed / norm;
	view.y += vector.y * speed / norm;
	view.z += vector.z * speed / norm;

};

// allows the camera to be pointed using the mouse
void Camera::mouseMove(int mpx, int mpy)
{
	// rotates wrt (0,1,0) as the up vector, will generalize at another time or if required
	// get the center of the window
	int middleX = screenwidth  >> 1;
	int middleY = screenheight >> 1;

	float dy;
	float ry;

	if( (mpx == middleX) && mpy == middleY )
		return;

	ry = (float)( (middleX - mpx) ) / (MOUSEROTATEFACTOR); // side
	dy  = (float)( (middleY - mpy) ) / (MOUSEROTATEFACTOR); // up/down

	view.y += dy * 2;


	//cap the maximum it can move, quick and dirty
	if((view.y - position.y) > 10)
		view.y = position.y + 10;
	if((view.y - position.y) < -10)
		view.y = position.y - 10;

	rotate(0, -ry, 0);

};
