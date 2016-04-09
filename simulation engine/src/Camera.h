#ifndef camera_h
#define camera_h

#include "vector3f.h"
#include <cmath>
#include <GL/glut.h>

#define MOUSEROTATEFACTOR 350*0.8f

class Camera
{
	public:
		Vector3f position;
		Vector3f view;
		Vector3f up;

		// used for mouse movement of the camera
		int screenwidth;
		int screenheight;

		Camera(void);
		void setPosition(float xposition, float yposition, float zposition,  float xview, float yview, float zview, float upx, float upy, float upz);
		void setPosition(Vector3f pos, Vector3f view, Vector3f up);
		void rotate(float x, float y, float z);
		void pointRotate(Vector3f point, float x, float y, float z);
		void move(float speed);
		void mouseMove(int x, int y);
		void side(float speed);
		void setScreenMetric(int w, int h) {screenwidth=w;screenheight=h;};
};

#endif
