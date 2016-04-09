// openglvis.h
// this file contains all global functions for displaying whats going on in opengl,
// its non essential for the contant of this code so its all been moved here to keep
// it separate from the other stuff

#ifndef OPENGL_H
#define OPENGL_H

extern Replica replica[REPLICA_COUNT];

bool displayminimum = false;

#define WIDTH      	 	800
#define HEIGHT       	800
#define MOVESPEED  0.05f  //movement
#define MSPEED     0.02f  //mouse
#define MACC       0.8f  //mouse acceleration
#include <GL/glut.h>
#include <GL/gl.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include "Camera.h"
#include <vector>
#include <ctime>
#include <linux/kernel.h>
//#include <linux/time.h>

#define M_ROTATE		1001
#define M_TRANSLATE 	1002
#define M_SCALE			1003
#define M_QUIT			1004
#define M_COLOR			1005
#define M_ADD			1006
#define M_REMOVE		1007
#define M_SAVE			1008
#define M_LOAD			1009
#define M_DESELECT  	1010
#define M_SELECT    	1011
#define M_WIREFRAME 	1012
#define M_FILLED		1013
#define M_NORMALS		1014
#define M_AXES			1015
#define M_INVSELECTION 	1016

void GlutInit(int windowheight,int windowwidth, char* title);
void GlutDisplay();
void GlutResize(int width, int height);
void GlutMouseFunction(int btn, int state, int x, int y);
void GlutIdleFunction();
void GLUTMouseActiveMotion(int x, int y);
void GLUTMousePassiveMotion(int x, int y);
void GLUTMouse(int btn, int state, int x, int y);
void GLUTKeyboardPress(unsigned char x , int y ,int z);
void GLUTKeyboardUp(unsigned char x , int y ,int z);
void GLUTKeyboardSpecial(int x , int y ,int z);
void GLUTcreateMenu();
void GLUTmenu(int);
void GLUTdrawaxes();

Camera camera;

//modes
bool b_mlook = true;
bool b_rotate = false;
bool b_translate = false;
bool b_scale  = false;

//draw options
bool drawAxes = false;
bool drawNormals = false;
bool drawFilled = true;
bool drawWire = false;
bool pulse = false;
bool glow = false;
bool fog = false;

int mousex;
int mousey;
int windowh;
int windoww;
bool showSelection = true;
int gl_replicaCount;
int gl_boundingValue;
int dr = 0; //replica to display
// input control variables 1,0,-1 indicate states
int moveFrontback = 0;
int moveLeftRight = 0;
float deltaRotate = 0;
unsigned char lastkey = 0;
int glutWindowID;

float LightPos0[] = { -100.0f, 100.0f, 100.0f, 1.0f};			// Light Position
float LightPos1[] = { -100.0f, -100.0f, 100.0f, 1.0f};			// Light Position
float LightPos2[] = { -100.0f, 0.0f, -100.0f, 1.0f};			// Light Position
float materialShininess = 128.0f;										// Material shininess

float COLOURS[5][3] =
{{ 0.8f, 0.2f, 0.0f },
{ 0.0f, 0.7f, 0.7f },
{ 0.0f, 0.0f, 0.7f },
{ 1.0f, 0.5f, 0.0f },
{ 1.0f, 1.0f, 0.0f }};

char *windowtitle;
// initialization function
void GlutInit(int windowheight,int windowwidth, char* title)
{
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(windowheight, windowwidth);
	windowh = windowheight;
	windoww = windowwidth;
	mousex = windowwidth >> 1;
	mousey = windowheight >> 1;
	camera.setScreenMetric(windowwidth,windowheight);
    glutWindowID = glutCreateWindow(title);
	windowtitle = title;
	GLUTcreateMenu();

	glutReshapeFunc(GlutResize);
	glutDisplayFunc(GlutDisplay);
	glutIdleFunc(GlutIdleFunction);

	glutMouseFunc(GLUTMouse);
	glutMotionFunc(GLUTMouseActiveMotion);
	glutPassiveMotionFunc(GLUTMousePassiveMotion);
	glutIgnoreKeyRepeat(false);
	glutKeyboardFunc(GLUTKeyboardPress);
	glutKeyboardUpFunc(GLUTKeyboardUp);
	glutSpecialFunc(GLUTKeyboardSpecial);

	float LightPos0[] = { 200.0f, 200.0f, 200.0f, 1.0f};			// Light Position
	float LightPos1[] = { -200.0f, 200.0f, 200.0f, 1.0f};			// Light Position
	float LightPos2[] = { 150.0f, 200.0f, -180.0f, 1.0f};			// Light Position
	float LightAmb[] = { 0.3f, 0.3f, 0.3f, 1.0f};					// Ambient Light Values
	float LightDif[] = { 0.7f, 0.7f, 0.7f, 1.0f};					// Diffuse Light Values
	float LightSpc[] = {0.5f, 0.5f, 0.5f, 1.0f};					// Specular Light Values
	float materialAmbient[] = {0.2f, 0.2f, 0.2f, 0.0f};						// Material ambient values
	float materialDiffuse[] = {0.27f, 0.27f, 0.27f, 0.0f};						// Material diffuse values
	float materialSpecular[] = {1.0f, 1.0f, 1.0f, 0.0f};						// Material specular values

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_NORMALIZE);

	glViewport(0,0,windowwidth,windowheight);
	glMatrixMode(GL_PROJECTION);
	gluPerspective(65.0f, 1.0, 0.001f, 100.0f);

	glEnable(GL_CULL_FACE);

	//glClearColor(0.25,0.25,0.25,0);
	glClearColor(1,1,1,0);
	glMatrixMode(GL_MODELVIEW);

	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
	glShadeModel(GL_SMOOTH);

	glEnable(GL_LIGHTING);
	glLightfv(GL_LIGHT0,GL_POSITION,LightPos0);
	glLightfv(GL_LIGHT0,GL_AMBIENT,LightAmb);
	glLightfv(GL_LIGHT0,GL_SPECULAR, LightSpc);
	glLightf(GL_LIGHT0,GL_CONSTANT_ATTENUATION,1.25f);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, LightDif);
	glEnable(GL_LIGHT0);

	glLightfv(GL_LIGHT1,GL_POSITION,LightPos1);
	glLightfv(GL_LIGHT1,GL_AMBIENT,LightAmb);
	glLightfv(GL_LIGHT1,GL_SPECULAR, LightSpc);
	glLightf(GL_LIGHT1,GL_CONSTANT_ATTENUATION,1.25f);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, LightDif);
	glEnable(GL_LIGHT1);

	glLightfv(GL_LIGHT2,GL_POSITION,LightPos2);
	glLightfv(GL_LIGHT2,GL_AMBIENT,LightAmb);
	glLightfv(GL_LIGHT2,GL_SPECULAR, LightSpc);
	glLightf(GL_LIGHT2,GL_CONSTANT_ATTENUATION,1.25f);
	glLightfv(GL_LIGHT2, GL_DIFFUSE, LightDif);
	glEnable(GL_LIGHT2);

	glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER,GL_TRUE);

	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT,GL_AMBIENT);
	glMaterialfv(GL_FRONT, GL_AMBIENT, materialAmbient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, materialDiffuse );
	glMaterialfv(GL_FRONT, GL_SPECULAR, materialSpecular );
	glMaterialf( GL_FRONT, GL_SHININESS, materialShininess );
}

#define AXISDIM 1500
void GLUTdrawBoundingSphere()
{
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glPushMatrix();
		glDisable(GL_LIGHTING);
		glDisable(GL_FOG);
		/*glBegin(GL_LINES);
			glLineWidth(2.0f);

			// x axis
			glColor3f(0.15,0.15,0.15);
			glVertex3f(-AXISDIM,0,0);

			glColor3f(1.0,0,0);
			glVertex3f(0,0,0);
			glVertex3f(0,0,0);

			glColor3f(0.15,0.15,0.15);
			glVertex3f(AXISDIM,0,0);

			// y axis
			glColor3f(0.15,0.15,0.15);
			glVertex3f(0,-AXISDIM,0);

			glColor3f(0,1.0,0);
			glVertex3f(0,0,0);
			glVertex3f(0,0,0);

			glColor3f(0.15,0.15,0.15);
			glVertex3f(0,AXISDIM,0);

			// z axis
			glColor3f(0.15,0.15,0.15);
			glVertex3f(0,0,-AXISDIM);

			glColor3f(0,0,1.0);
			glVertex3f(0,0,0);
			glVertex3f(0,0,0);

			glColor3f(0.15,0.15,0.15);
			glVertex3f(0,0,AXISDIM);

			glLineWidth(0.5f);
			glColor4f(0.75,0.75,0.75,0.2);
			for (int i=0; i<=gl_boundingValue;i+=2)
			{
				float l = gl_boundingValue * cos ( asin(i/gl_boundingValue));
				glVertex3f(-l,0,i);
				glVertex3f(l,0,i);
				glVertex3f(-l,0,-i);
				glVertex3f(l,0,-i);
				glVertex3f(i,0,-l);
				glVertex3f(i,0,l);
				glVertex3f(-i,0,-l);
				glVertex3f(-i,0,l);
			}

		glEnd();
*/
		// displays bounding sphere
		glBegin(GL_LINE_LOOP);
			glLineWidth(2.0f);
			glColor3f(1.0,1.0,0.0);
			for (int i=0;i<=628;i++)
			{
				float x = float(i)/100.0f;
				glVertex3f(sin(x)*gl_boundingValue,0,cos(x)*gl_boundingValue);

			}
		glEnd();
		GLUquadric * quadric = gluNewQuadric();
		glPushMatrix();
			glTranslatef(0,0,0);
			//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			glEnable (GL_BLEND); glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			glColor4f(0.75,0.75,0.75,0.2);

			glRotatef(90,1,0,0);
			gluSphere(quadric, gl_boundingValue, 30, 30);
		glPopMatrix();
		gluDeleteQuadric(quadric);
		glEnable(GL_LIGHTING);
	glPopMatrix();
	glPopAttrib();
}

void GLUTdrawBoundingBox()
{
	float lindim  = gl_boundingValue*0.5;
	glPushAttrib(GL_ALL_ATTRIB_BITS);
		glPushMatrix();
			glDisable(GL_LIGHTING);
			glDisable(GL_FOG);
			/*glBegin(GL_LINES);
				glLineWidth(1.0f);
				glColor3f(0.75,0.75,0.75);
				for (int i=0; i<=gl_boundingValue;i+=5)
					for (int j=0; j<=gl_boundingValue;j+=5)
					{
						glVertex3f(-lindim,0,j-lindim);
						glVertex3f(lindim,0,j-lindim);
						glVertex3f(i-lindim,0,-lindim);
						glVertex3f(i-lindim,0,lindim);
					}

			glEnd();
*/
			//glScalef(1e10,1e10,1e10);

			// displays bounding box
			glBegin(GL_LINES);
				glLineWidth(5.0f);
				//glColor3f(1.0,1.0,1.0);
				glColor3f(0,0,0);
					glVertex3f(-lindim,-lindim,-lindim);
					glVertex3f(lindim,-lindim,-lindim);

					glVertex3f(-lindim,-lindim,-lindim);
					glVertex3f(-lindim,-lindim,lindim);

					glVertex3f(-lindim,-lindim,-lindim);
					glVertex3f(-lindim,lindim,-lindim);

					glVertex3f(lindim,lindim,lindim);
					glVertex3f(-lindim,lindim,lindim);

					glVertex3f(lindim,lindim,lindim);
					glVertex3f(lindim,-lindim,lindim);

					glVertex3f(lindim,lindim,lindim);
					glVertex3f(lindim,lindim,-lindim);


					glVertex3f(lindim,-lindim,lindim);
					glVertex3f(-lindim,-lindim,lindim);

					glVertex3f(lindim,-lindim,lindim);
					glVertex3f(lindim,-lindim,-lindim);

					glVertex3f(lindim,-lindim,-lindim);
					glVertex3f(lindim,lindim,-lindim);


					glVertex3f(lindim,lindim,-lindim);
					glVertex3f(-lindim,lindim,-lindim);

					glVertex3f(-lindim,lindim,-lindim);
					glVertex3f(-lindim,lindim,lindim);

					glVertex3f(-lindim,lindim,lindim);
					glVertex3f(-lindim,-lindim,lindim);

			glEnd();

			GLUquadric * quadric = gluNewQuadric();
			glPushMatrix();
				glTranslatef(0,0,0);
				//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
				glEnable (GL_BLEND); glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
				//glColor4f(0.75,0.75,0.75,0.2);
				glColor4f(0.25,0.25,0.25,0.15);
					glRotatef(90,1,0,0);
					glutSolidCube(gl_boundingValue);
					glPopMatrix();
					gluDeleteQuadric(quadric);
			glEnable(GL_LIGHTING);
		glPopMatrix();
		glPopAttrib();
}


// glut display call
void GlutDisplay()
{
	GLUquadric * quadric;
	float lindim  = gl_boundingValue*0.5;

	glLoadIdentity();
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	gluLookAt(camera.position.x,camera.position.y,camera.position.z, camera.view.x,camera.view.y,camera.view.z, camera.up.x,camera.up.y,camera.up.z);

	glLightfv(GL_LIGHT1,GL_POSITION,LightPos0);
	glLightfv(GL_LIGHT1,GL_POSITION,LightPos1);
	glLightfv(GL_LIGHT1,GL_POSITION,LightPos2);

	glScaled(0.1,0.1,0.1);

	// drawing protein

	quadric = gluNewQuadric();
	Vector3f positionS;

	glPushAttrib(GL_ALL_ATTRIB_BITS);
		glPushMatrix();
			glEnable(GL_DEPTH_TEST);
			glDisable(GL_BLEND);


			/*glPushMatrix();
				glColor3f(1.0,0.0,0.0);
				glTranslated (1,0,0);
				gluSphere(quadric, 0.15, 20, 20);
			glPopMatrix();
			glPushMatrix();
				glColor3f(0.0,1.0,0.0);
				glTranslated (0,1,0);
				gluSphere(quadric, 0.15, 20, 20);
			glPopMatrix();
			glPushMatrix();
				glColor3f(0.0,0.0,1.0);
				glTranslated (0,0,1);
				gluSphere(quadric, 0.15, 20, 20);
			glPopMatrix();	*/


			for(size_t m=0;m<GLreplica->moleculeCount;m++)
			{
				positionS = GLreplica->molecules[m].center;
				#if BOUNDING_METHOD == PERIODIC_BOUNDRY
					positionS = positionS + Vector3f(-lindim,-lindim,-lindim);
				#endif

				glColor3f(1.0,1.0,0.0);
				glBegin(GL_LINES);
					glLineWidth(3.0f);
					glVertex3f(positionS.x,0,positionS.z);
					glVertex3f(positionS.x,positionS.y,positionS.z);
				glEnd();

				glPushMatrix();
					glTranslatef(positionS.x,positionS.y,positionS.z);
					glCullFace(GL_BACK);
					//glColor4d ( 1.0f, 1.0f, 0.0f, 1.0f);
					//gluSphere(quadric, 20, 12, 10);
					if (m==0)
					{
						glColor4d ( 0.0f, 0.7f, 0.1f, 1.0f);
					}
					else if (m==1)
					{
						glColor4d ( 0.0f, 0.1f, 0.7f, 1.0f);
					}
					else
					{
						//glDisable(GL_DEPTH_TEST);
						glEnable (GL_BLEND);
						glBlendFunc (GL_SRC_ALPHA,GL_SRC_COLOR);
						glColor4d ( 0.8f, 0.1f, 0.1f, 0.95f);
					}

					for (uint r=0;r<GLreplica->molecules[m].residueCount;r++)
					{
						glPushMatrix();
							Vector3f p = GLreplica->molecules[m].Residues[r].relativePosition;
							glTranslatef(p.x,p.y,p.z);
							// gluSphere(quadric, radius, slices, rings)
							gluSphere(quadric, GLreplica->molecules[m].Residues[r].vanderWaalRadius, 12, 10);
						glPopMatrix();
					}
				glPopMatrix();
			}


			/*  glPushMatrix();
				glColor3d(1.0,0,0);
				glBegin(GL_LINES);
					glVertex3d(replica[dr].pResidues[0].position.x,replica[dr].pResidues[0].position.y,replica[dr].pResidues[0].position.z);
					for (int i=1;i<replica[dr].length;i++)
					{
							glColor3d(1.0-i*1.0/(double(replica[dr].length)),0.0,i*1.0/(double(replica[dr].length)));
							glVertex3d(replica[dr].pResidues[i].position.x,replica[dr].pResidues[i].position.y,replica[dr].pResidues[i].position.z);
							//ated(replica[0].pResidues[i].position.x-replica[0].pResidues[i-1].position.x,
									//	 replica[0].pResidues[i].position.y-replica[0].pResidues[i-1].position.y,
									//	 replica[0].pResidues[i].position.z-replica[0].pResidues[i-1].position.z);
							// gluCylinder(quadric, bottomRadius, topRadius, height, slices, rings)
							//gluCylinder(quadric, 1e-10,  1e-10, (replica[0].pResidues[i].position-replica[0].pResidues[i-1].position).magnitude(), 7, 7);
					}
				glEnd();
			glPopMatrix();*/


	//	glPopMatrix();
	//glPopAttrib();
    //glPushAttrib(GL_ALL_ATTRIB_BITS);
	//	glPushMatrix();
	//		glEnable(GL_DEPTH_TEST);

		// show minimum
/*		if (displayminimum)
		{
			//glScaled(0.1,0.1,0.1);

			for(size_t m=0;m<replica[dr].moleculeCount;m++)
			{
				positionS = replica[dr].molecules[m].center;

				glColor3f(1.0,1.0,0.0);
				glBegin(GL_LINES);
					glLineWidth(3.0f);
					glVertex3f(positionS.x,0,positionS.z);
					glVertex3f(positionS.x,positionS.y,positionS.z);
				glEnd();

				glPushMatrix();
					glTranslatef(positionS.x,positionS.y,positionS.z);
					glColor3f(0.7,0.7,0.7);
				    for (uint r=0;r<replica[lowestEnergy].molecules[m].residueCount;r++)
					{
						glPushMatrix();
							Vector3f p = replica[lowestEnergy].molecules[m].Residues[r].relativePosition;
							glTranslatef(p.x,p.y,p.z);
							// gluSphere(quadric, radius, slices, rings)
							gluSphere(quadric, 2, 20, 20);
						glPopMatrix();
					}
				glPopMatrix();
			}
		}*/
		glPopMatrix();
	glPopAttrib();

	gluDeleteQuadric(quadric);
#if GL_AXES
	#if BOUNDING_METHOD == BOUNDING_SPHERE
		GLUTdrawBoundingSphere();
	#else
		GLUTdrawBoundingBox();
	#endif
#endif
	glFlush();
	glutSwapBuffers();
}

// idle glut function
void GlutIdleFunction()
{
	// introduce time dependent functions
	static long last_time =  clock()-100;
	long current_time = clock();

	float timescale = ((float)(current_time-last_time))/1000.0f;

	if (moveFrontback!=0)
		camera.move(MOVESPEED*timescale*moveFrontback);

	if (moveLeftRight!=0)
		camera.side(MOVESPEED*timescale*moveLeftRight);

	last_time = current_time;
	glutPostRedisplay();
}

void GlutResize(int w, int h)
{
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0f,(GLfloat)w/(GLfloat)h,0.1f,1000.0f);
	glMatrixMode(GL_MODELVIEW);
	gluLookAt(-2.0f,0.0f,-5.0f, -2.0f,0.0f,-4.0f , 0.0f,1.0f,0.0f);
}

void GLUTMouseActiveMotion(int x, int y)
{
	if (b_mlook)
	{
		int middleX = windoww  >> 1;
		int middleY = windowh >> 1;

		if (x<middleX)
			camera.pointRotate(Vector3f(0,0,0),0,1*MSPEED,0);
		if (x>middleX)
			camera.pointRotate(Vector3f(0,0,0),0,-1*MSPEED,0);
		if (y<middleY)
			camera.pointRotate(Vector3f(0,0,0),1*MSPEED,0,0);
		if (y>middleY)
			camera.pointRotate(Vector3f(0,0,0),-1*MSPEED,0,0);


		//camera.mouseMove(x,y);
		glutWarpPointer(windoww >> 1, windowh >> 1);
	}
}

void GLUTMousePassiveMotion(int x, int y)
{

}

void GLUTMouse(int btn, int state, int x, int y)
{
	if(btn==GLUT_LEFT_BUTTON && state == GLUT_DOWN)
	{
		mousex = x;
		mousey = y;
		glutWarpPointer(windoww >> 1, windowh >> 1);
		b_mlook = true;
		return;
	}
	if(btn==GLUT_LEFT_BUTTON && state == GLUT_UP)
	{
		if (b_mlook)
		{
			glutWarpPointer(mousex,mousey);
			moveFrontback = 0;
			moveLeftRight = 0;
			//ShowCursor(true);
			b_mlook = false;
		}
		return;
	}
}

void GLUTKeyboardPress(unsigned char key, int x, int y)
{
	switch(key)
		{
			case 'w':case 'W':
				camera.pointRotate(Vector3f(0,0,0),1*MOVESPEED,0,0);
				break;
			case 's':case 'S':
				camera.pointRotate(Vector3f(0,0,0),-1*MOVESPEED,0,0);
				break;
			case 'a':case 'A':
				camera.pointRotate(Vector3f(0,0,0),0,1*MOVESPEED,0);
				break;
			case 'd':case 'D':
				camera.pointRotate(Vector3f(0,0,0),0,-1*MOVESPEED,0);
				break;
			case 'n':case 'N':
				dr++;
				dr = dr%gl_replicaCount;
				GLreplica = & replica[dr];
				cout << "displaying replica: " << dr << " E=" << replica[dr].potential << endl;
				break;
			case 'm':case 'M':
				displayminimum = !displayminimum;
				break;
			case 'q':case 'Q':
				camera.position = camera.position * 1.01f;
				break;
			case 'e':case 'E':
				camera.position = camera.position * 0.99f;
				break;
		}

}

void GLUTKeyboardUp(unsigned char key, int x, int y)
{
	switch(key)
	{
		case 'w': case 'W':
			moveFrontback = 0;
			break;
		case 's': case 'S':
			moveFrontback = 0;
			break;
		case 'a': case 'A':
			moveLeftRight = 0;
			break;
		case 'd': case 'D':
			moveLeftRight = 0;
			break;
	}
}

void GLUTKeyboardSpecial(int key, int x, int y)
{
	switch(key)
	{
		case GLUT_KEY_F5 :
			//saveScene();
			break;

		case GLUT_KEY_F9 :
			//loadScene();
			break;
	}
}

void GLUTmenu(int call)
{
	b_mlook = true;
	switch(call)
	{
		case M_QUIT:
			glutDestroyWindow(glutWindowID);
			exit(0);
			break;
	}
}

void GLUTcreateMenu()
{
	glutCreateMenu(GLUTmenu);
	glutAddMenuEntry("Quit", M_QUIT);
	glutAttachMenu(GLUT_RIGHT_BUTTON);
}
// gl stuff ends

#endif
