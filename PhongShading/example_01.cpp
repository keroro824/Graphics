#include <stdlib.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

#ifdef OSX
#include <GLUT/glut.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#include <GL/glu.h>
#endif

#include <time.h>
#include <math.h>


#define PI 3.14159265  // Should be used from mathlib
inline float sqr(float x) { return x*x; }

using namespace std;

//****************************************************
// Some Classes
//****************************************************

class Viewport;
class Color;

class Viewport {
  public:
    int w, h; // width and height
};

class Vector3D {
public:

    float x, y, z;

    Vector3D (float X=0.0, float Y=0.0, float Z=0.0) {
        float ary[] = {X, Y, Z};
        std::vector<float> color (ary, ary+sizeof(ary)/sizeof(float));
        x = X;
        y = Y;
        z = Z;
    }
    
    Vector3D operator* (Vector3D const &v) {
        return Vector3D(x * v.x, y * v.y, z * v.z);
    }
    
    Vector3D operator* (float scaler) {
        return Vector3D(x * scaler, y * scaler, z * scaler);
    }
    
    Vector3D operator+ (Vector3D const &v) {
        return Vector3D(x + v.x, y + v.y,  z + v.z);
    }
    
    Vector3D operator- (Vector3D const &v) {
        return Vector3D(x - v.x, y - v.y, z - v.z);
    }
    
    float dotVector3D (const Vector3D &v) {
        return (float) x * v.x + y * v.y + z * v.z;
    }
    
    void rgb_normalize() {
        if (x!=0.0f || y!=0.0f || z!=0.0f) {
            x = x/(x+y+z);
            y = y/(x+y+z);
            z = z/(x+y+z);
        }
    }
    
    void xyz_normalize() {
        if (x!=0.0f || y!=0.0f || z!=0.0f) {
            float sc = sqrt(x*x+y*y+z*z);
            x = (x)/sc;
            y = (y)/sc;
            z = (z)/sc;
        }
    }
};


//****************************************************
// Global Variables
//****************************************************
Viewport	viewport;

int pt_counter = 0;
int dl_counter = 0;
int multis_counter = 0;

Vector3D pt_rgb[5];
Vector3D pt_xyz[5];

Vector3D dl_rgb[5];
Vector3D dl_xyz[5];
Vector3D multis_xyz[10];
Vector3D multis_r[10];


//Initialization
Vector3D ka = Vector3D(0.0f, 0.0f, 0.0f);
Vector3D kd = Vector3D(0.0f, 0.0f, 0.0f);
Vector3D ks = Vector3D(0.0f, 0.0f, 0.0f);
float p;
bool multis = false;


//****************************************************
// Simple init function
//****************************************************
void initScene(){

  // Nothing to do here for this simple example.

}


//****************************************************
// reshape viewport if the window is resized
//****************************************************
void myReshape(int w, int h) {
  viewport.w = w;
  viewport.h = h;

  glViewport (0,0,viewport.w,viewport.h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0, viewport.w, 0, viewport.h);

}


//****************************************************
// A routine to set a pixel by drawing a GL point.  This is not a
// general purpose routine as it assumes a lot of stuff specific to
// this example.
//****************************************************

void setPixel(int x, int y, GLfloat r, GLfloat g, GLfloat b) {
  glColor3f(r, g, b);
  glVertex2f(x + 0.5, y + 0.5);   // The 0.5 is to target pixel
  // centers 
  // Note: Need to check for gap
  // bug on inst machines.
}

//****************************************************
// Draw a filled circle.  
//****************************************************



Vector3D diffuseTerm (Vector3D kd, Vector3D intens, Vector3D normal, Vector3D light){
    float dotProduct = max(light.dotVector3D(normal), 0.0f);
    if (dotProduct == 0.0f) {
        return Vector3D(0.0f,0.0f,0.0f);
    } else {
        return (kd*intens)*dotProduct;
    }
};


Vector3D specularTerm (Vector3D ks, Vector3D intens, Vector3D normal, Vector3D light, Vector3D viewer, float p){
    Vector3D refl = light*(-1.0f) + normal*(2.0f*(light.dotVector3D(normal)));
    float dotProduct = max(refl.dotVector3D(viewer), 0.0f);
    if (dotProduct == 0.0f){
        return Vector3D(0.0f,0.0f,0.0f);
    } else {
        return (ks * intens) * pow (dotProduct, p);
    }
};



void circle(float centerX, float centerY, float radius) {
  // Draw inner circle
  glBegin(GL_POINTS);

  // We could eliminate wasted work by only looping over the pixels
  // inside the sphere's radius.  But the example is more clear this
  // way.  In general drawing an object by loopig over the whole
  // screen is wasteful.

    int i,j;  // Pixel indices

    int minI = max(0,(int)floor(centerX-radius));
    int maxI = min(viewport.w-1,(int)ceil(centerX+radius));

    int minJ = max(0,(int)floor(centerY-radius));
    int maxJ = min(viewport.h-1,(int)ceil(centerY+radius));
    

    for (i=0;i<viewport.w;i++) {
        for (j=0;j<viewport.h;j++) {

            // Location of the center of pixel relative to center of sphere
            float x = (i+0.5-centerX);
            float y = (j+0.5-centerY);

            float dist = sqrt(sqr(x) + sqr(y));

            if (dist<=radius) {
                
                if (pt_counter != 0.0f || dl_counter != 0.0f){
                    
                    Vector3D viewer = Vector3D(0.0f,0.0f,1.0f);

                    // This is the front-facing Z coordinate
                    float z = sqrt(radius*radius-dist*dist);

                    Vector3D normal = Vector3D(x, y, z);
                    normal.xyz_normalize();
//                    cout<<normal.z<<endl;
                    
                    
                    
                    // Default color of the sphere is black
                    Vector3D R = Vector3D(0.0f, 0.0f, 0.0f);
                    

                    if (pt_counter != 0) {
                        for (int l = 0; l < pt_counter; l++) {
                            Vector3D pt_light_xyz = pt_xyz[l];
                            Vector3D pt_light_rgb = pt_rgb[l];
                            
                            Vector3D light = pt_light_xyz - normal;
                            light.xyz_normalize();
                            
                            Vector3D diffuse = diffuseTerm(kd, pt_light_rgb, normal, light);
                            
                            Vector3D specular = specularTerm(ks, pt_light_rgb, normal, light, viewer, p);

                            Vector3D ambient = ka*pt_light_rgb;

                            R = R + (diffuse + specular + ambient);
//                            if (normal.x>0.99370 ){
//                                cout<<R.z<<endl;
//                                break;
//                            }
                            
                            
                            
                        }
                    }
                    
                    if (dl_counter != 0) {
                        for (int l = 0; l < dl_counter; l++) {
                            Vector3D dl_light_rgb = dl_rgb[l];
                            
                            Vector3D light = dl_xyz[l];
                            light.xyz_normalize();
                            
                            Vector3D diffuse = diffuseTerm(kd, dl_light_rgb, normal, light);
                            
                            Vector3D specular = specularTerm(ks, dl_light_rgb, normal, light, viewer, p);

                            Vector3D ambient = ka*dl_light_rgb;

                            R = R + (diffuse + specular + ambient);
                        }
                    }
                    
                    setPixel(i,j, R.x, R.y, R.z);
                
                }


            // This is amusing, but it assumes negative color values are treated reasonably.
            // setPixel(i,j, x/radius, y/radius, z/radius );
            }
        }
  }


  glEnd();
}
//****************************************************
// function that does the actual drawing of stuff
//***************************************************
void myDisplay() {

  glClear(GL_COLOR_BUFFER_BIT);				// clear the color buffer

  glMatrixMode(GL_MODELVIEW);			        // indicate we are specifying camera transformations
  glLoadIdentity();				        // make sure transformation is "zero'd"


  // Start drawing
  
 
  if (multis){
      float x,y,r;
      for(int i=0;i<multis_counter;i++){
          x = multis_xyz[i].x;
          y = multis_xyz[i].y;
         
          r = multis_r[i].x;
          circle(viewport.w / 2.0 + x, viewport.h / 2.0 +y, r);
      }
    }
  else{
       circle(viewport.w / 2.0 , viewport.h / 2.0 , min(viewport.w, viewport.h) / 3.0);
  }

  glFlush();
  glutSwapBuffers();					// swap buffers (we earlier set double buffer)
}


float convertor(char *argv){
    std::string strin= string(argv);
    return atof( strin.c_str());
}


void keyboard (unsigned char key, int x, int y)
{
    switch (key) {
        case 32:
            exit(0);
            break;
            
        default:
            break;
    }
    glutPostRedisplay();
}



//****************************************************
// the usual stuff, nothing exciting here
//****************************************************
int main(int argc, char *argv[]) {
//****************************************************************
// Things to do in parse: construct vectors for global variables *
//****************************************************************

    //This is parsing the command line
    for (int i = 1; i < argc; ++i) {
        float r,g,b,x,y,z;
        
        
        if (string(argv[i]) == "-ka"){
            r=convertor(argv[i+1]);
            g=convertor(argv[i+2]);
            b=convertor(argv[i+3]);
            ka = Vector3D(r,g,b);
            i+=3;
        }
        else if (string(argv[i]) == "-kd"){
            r=convertor(argv[i+1]);
            g=convertor(argv[i+2]);
            b=convertor(argv[i+3]);
            kd = Vector3D(r,g,b);
            i+=3;
        }
        else if (string(argv[i]) == "-ks"){
            r=convertor(argv[i+1]);
            g=convertor(argv[i+2]);
            b=convertor(argv[i+3]);
            ks = Vector3D(r,g,b);
            i+=3;
        }
        else if (string(argv[i]) == "-sp"){
            p=convertor(argv[i+1]);
            i+=1;
        }
        else if (string(argv[i]) == "-pl"){
            x=convertor(argv[i+1]);
            y=convertor(argv[i+2]);
            z=convertor(argv[i+3]);
            r=convertor(argv[i+4]);
            g=convertor(argv[i+5]);
            b=convertor(argv[i+6]);
            pt_xyz[pt_counter] = Vector3D(x,y,z);
            pt_rgb[pt_counter] = Vector3D(r,g,b);
            pt_counter+=1;
            i+=6;
            
        }
        else if (string(argv[i]) == "-dl"){
            x=convertor(argv[i+1]);
            y=convertor(argv[i+2]);
            z=convertor(argv[i+3]);
            r=convertor(argv[i+4]);
            g=convertor(argv[i+5]);
            b=convertor(argv[i+6]);
            dl_xyz[dl_counter] = Vector3D(x,y,z)*(-1);
            dl_rgb[dl_counter] = Vector3D(r,g,b);
            dl_counter+=1;
            i+=6;
            
        }
        else if (string(argv[i]) == "-multis"){
            multis = true;
            x=convertor(argv[i+1]);
            y=convertor(argv[i+2]);
            r=convertor(argv[i+3]);
            multis_xyz[multis_counter] = Vector3D(x,y,0);
            multis_r[multis_counter] = Vector3D(r,0,0);
            multis_counter+=1;
            i+=3;
            
        }
    }

    //This initializes glut
    glutInit(&argc, argv);
    
    //This tells glut to use a double-buffered window with red, green, and blue channels
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);

    // Initalize theviewport size
    viewport.w = 400;
    viewport.h = 400;

    //The size and position of the window
    glutInitWindowSize(viewport.w, viewport.h);
    glutInitWindowPosition(0,0);
    glutCreateWindow(argv[0]);
    glutKeyboardFunc( keyboard );
    
    initScene();							// quick function to set up scene

    glutDisplayFunc(myDisplay);				// function to run when its time to draw something
    glutReshapeFunc(myReshape);				// function to run when the window gets resized


    glutMainLoop();							// infinite loop that will keep drawing and resizing
    // and whatever else


  return 0;
}









