//
//  RayTracer.cpp
//  
//
//  Developed by Ian Chen and Betty Chen at Mar 2013
//
//
//  version Mar10 21:02
//
//  added features:
//      -fix Up-direction, ray generator now looks at perspective view from anywhere
//      -add emission and reflection features to all objects: spheres, tris, and trinormals
//      -fix reflection, replace reflection weight as the specularity ks.
//
//  Problems:
//      attenuation not supported
//
//  Future Development:
//      -add accelerating features to the RayTracer, use bounding box for intersection test.
//      -add refraction feature to the RayTracer.
//
//


#include "FreeImage.h"
#include <iostream>
#include "Eigen/Core"
#include "Eigen/Dense"
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>

#define BPP 24
#define PI 3.1415926535897932384626

using namespace std;
using namespace Eigen;

//***********************************************
//  Utility Function Claim
//***********************************************

Vector4f normalize(Vector4f);
Vector4f error_avoid(Vector4f, Vector4f);
Vector4f cross_product(Vector4f, Vector4f);


//***********************************************
//  Classes
//***********************************************
class Color{
    public:
        float r;
        float g;
        float b;
        Vector3f intensity;
        Color (){
        }
        Color (float r, float g, float b){
            this->r = r;
            this->g = g;
            this->b = b;
            this->intensity = Vector3f(r, g, b);
        }
    ~ Color(){}

        Color operator+ (Color const &c){
        	return Color(r+c.r, g+c.g, b+c.b);
        } 
        Color operator* (Color const &c){
        	return Color(r*c.r, g*c.g, b*c.b);
        }
        Color operator* (float scaler){
        	return Color(r*scaler, g*scaler, b*scaler);
        }
};



class Ray{
    public:
        Vector4f origin;       
        Vector4f direction;     
        Ray(){
        }

        Ray(Vector4f camara, Vector4f direc){
            this->origin = camara;
            this->direction = direc;
        }

        Vector4f shootpoint(float t){
            return this->origin + this->direction * t;
        }

        float get_t(Vector4f point){
        	Vector4f dir = point -origin;
        	for (int i = 0; i < 3; i ++){
        		if (dir[i] !=0 && direction[i] !=0){
					float t = dir[i]/direction[i];
					if (t > 0){
						return t;
					} else {
						// printf(" Error: you have a negative t value. fuck off. \n");
					}
        		}
        	}
        }
};


class Material {
    public:
        Color ka;
        Color kr;
        Color kd;
        Color ks;
        float sp;
        Color emission;
        Material(){
            // do nothing
        }
        Material(Color ambient, Color diffuse, Color specular, float shiness, Color reflection = Color(0,0,0), Color emission = Color(0.0, 0.0, 0.0)){
            this->ka = ambient;
            this->kd = diffuse;
            this->ks = specular;
            this->sp = shiness;
            this->kr = reflection;
            this->emission = emission;
        }
};


class Vertex{
    public:
        Vector4f coordinate;
        Vertex(){
        }
        Vertex(Vector4f position) {
            this->coordinate = position;

        }
        Vector4f operator+ (Vertex vplus){
            Vector4f sum = coordinate + vplus.coordinate;
            sum[3] = 0;
            return sum;
        }

        Vector4f operator- (Vertex vminus){
            Vector4f dif = coordinate - vminus.coordinate;
            dif[3] = 0;
            return dif;
        }
};

class VertexNormal : public Vertex{
    public:
        Vector4f normal;
        VertexNormal(){
        }
        VertexNormal(Vector4f coordinate, Vector4f normal) : Vertex(coordinate){
            if (normal[3] != 0){
                cout << "note that the normal has to be a direction to maintain consistency" << endl;
            }
            this->normal = normal;
        }
};


class Surface{
    public:
        string name;
        Matrix4f transformation;
        Material material;

        // for spheres only;
        Vector4f center;
        float radius;

        // for regular triangles only
        Vertex v1;
        Vertex v2;
        Vertex v3;
        Vector4f normal;

        // for triangles with normals only
        VertexNormal vn1;
        VertexNormal vn2;
        VertexNormal vn3;

        Surface(){
            // do nothing
        }

        Surface(string name, Matrix4f transformation, Material m){
            // this transformation matrix is from object space to world space;
            this->name = name;
            this->transformation = transformation;
            this->material = m;
        }

        Vector4f getNormal(Vector4f point){

        }
};

class Sphere : public Surface{
    public:
        Sphere() : Surface(){
        }
        Sphere (Matrix4f transformation_value, Vector4f center, float radius, Material m) : Surface ("Sphere", transformation_value, m){
            this->radius = radius;
            this->center = center;
        }
    ~Sphere(){}
};

class triangle : public Surface{
    public:
        triangle(){
        }
        triangle(Vertex v1, Vertex v2, Vertex v3, Matrix4f transformation, Material m) : Surface("Triangle", transformation, m){
            this->v1 = v1;
            this->v2 = v2;
            this->v3 = v3;

            Vector4f a1 = v2 - v1;
            Vector4f a2 = v3 - v1;
            Vector4f normal_direction = cross_product(a1, a2);
            this->normal = normalize(normal_direction);
        }
};


float find_area(Vector4f v1, Vector4f v2, Vector4f v3){
    Vector4f a = v2 - v1;
    Vector4f b = v3 - v1;
    Vector4f cx = cross_product(a, b);
    float area = sqrt(pow(cx[0], 2)+pow(cx[1], 2)+pow(cx[2], 2));
    return area/2;
};

class triangleNormal : public Surface{
    public:
        triangleNormal(){
        }
        triangleNormal(VertexNormal vn1, VertexNormal vn2, VertexNormal vn3, Matrix4f transformation, Material m) : Surface("TriangleNormal", transformation, m){
            this->vn1 = vn1;
            this->vn2 = vn2;
            this->vn3 = vn3;
        }



        Vector4f getNormal(Vector4f point){
            float AreaT = find_area(vn1.coordinate, vn2.coordinate, vn3.coordinate);
            float AreaB = find_area(vn1.coordinate, point, vn3.coordinate);
            float AreaC = find_area(vn1.coordinate, point, vn2.coordinate);
            float AreaA = find_area(vn2.coordinate, point, vn3.coordinate);

            float c1 = AreaA/AreaT;
            float c2 = AreaB/AreaT;
            float c3 = AreaC/AreaT;

            float nx = vn1.normal[0] * c1 + vn2.normal[0] * c2 + vn3.normal[0] * c3;
            float ny = vn1.normal[1] * c1 + vn2.normal[1] * c2 + vn3.normal[1] * c3;
            float nz = vn1.normal[2] * c1 + vn2.normal[2] * c2 + vn3.normal[2] * c3;

            Vector4f result = Vector4f(nx, ny, nz, 0);
            return result;
            
        }
};


class transformation_stack{
    public:
        std::vector<Matrix4f> storage;

        transformation_stack(){
            Matrix4f identical; identical << 1.0f, 0.0f, 0.0f, 0.0f,
                                        0.0f, 1.0f, 0.0f, 0.0f,
                                        0.0f, 0.0f, 1.0f, 0.0f,
                                        0.0f, 0.0f, 0.0f, 1.0f;
            storage.push_back(identical);
        }

        void push_transformation(Matrix4f trans){
            Matrix4f cur_transformation = storage.back();
            Matrix4f new_transformation = cur_transformation*trans;
            storage.pop_back();

            storage.push_back(new_transformation);

        }

        void push_transformation(){
            Matrix4f new_top = storage.back();
            storage.push_back(new_top);
        }

        Matrix4f get_transformation(){
            Matrix4f identical; identical << 1.0f, 0.0f, 0.0f, 0.0f,
                                        0.0f, 1.0f, 0.0f, 0.0f,
                                        0.0f, 0.0f, 1.0f, 0.0f,
                                        0.0f, 0.0f, 0.0f, 1.0f;
            if (storage.size() == 0){
                return identical;
            }
            Matrix4f result = storage.back();
            return result;
        }

        Matrix4f pop_transformation(){
            storage.pop_back();

        }
};




//     ~triangle();

//     /* data */
// };

//***********************************************
//  Gloable variables
//***********************************************
// Global Constant
const Color skycolor = Color(0.0, 0.0, 0.0);

// Global Marker
int hit_index = -1;

// Scene Specification
int TraceDepth = 5;
int WIDTH;
int HEIGHT;
float aspect_ratio;
float fov;
transformation_stack tfmStack;
string fname;

// Geometry Specification
int obj_counter;
int maxvertex;
int maxvertexnormal;
std::vector<Vertex> vertices;
std::vector<VertexNormal> verticesNormal;
std::vector<Surface> objects;

// Light Source Specification
int pt_light_counter;
int dl_light_counter;
std::vector<Color> pt_rgb;
std::vector<Vector4f> pt_xyz;
std::vector<Color> dl_rgb;
std::vector<Vector4f> dl_xyz;

// Camara Specification
Vector4f camara_position; 
Vector4f camara_looking_direction;
Vector4f camara_up_direction;
Vector4f camara_right_direction;

// Input Material
Color Ambient_input;
Color Diffuse_input;
Color Specular_input;
float Shiness_input;
Color reflection_input;
Color emission_input;


//***********************************************
//  Utility functions
//***********************************************

Vector4f error_avoid(Vector4f point, Vector4f normal){
    return point + normal*0.001;
};

Vector4f normalize(Vector4f v) {
	// if (v[3] != 0){
	// 	cout << "the vector you want to normalize is not a direction, v is given by: " << v << endl;
	// }

    if (v[0] !=0|| v[1] !=0 || v[2] !=0) {
        float sc = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
        v[0] = (v[0])/sc;
        v[1] = (v[1])/sc;
        v[2] = (v[2])/sc;
        v[3] = 0;
    }
    return v;
};

Vector4f cross_product(Vector4f a, Vector4f b){
	Vector3f aa = Vector3f(a[0], a[1], a[2]);
	Vector3f bb = Vector3f(b[0], b[1], b[2]);
	Vector3f product = aa.cross(bb);
	Vector4f result = Vector4f(product[0], product[1], product[2], 0);
	return result;
};

Matrix4f translation(Vector3f a){
    Matrix4f translate; translate<<1,0,0,a[0],0,1,0,a[1],0,0,1,a[2],0,0,0,1;
    return translate;
};

Matrix4f scale(Vector3f scaler){
    Matrix4f scales; scales<<scaler[0],0,0,0,0,scaler[1],0,0,0,0,scaler[2],0,0,0,0,1;
    return scales;
};

Matrix4f rotate(Vector4f rotate){
    float x = rotate[0];
    float y = rotate[1];
    float z = rotate[2];
    float angle = (float)rotate[3]/180*PI;
    Matrix4f I; I<<1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1;
    if (x==1.0){
        Matrix4f mx; mx<<1,0,0,0,0,cos(angle),(-sin(angle)),0,0, sin(angle), cos(angle),0,0,0,0,1;
        I = mx*I;
    }
    if (y==1.0){
        Matrix4f my; my<<cos(angle),0,sin(angle),0,0,1,0,0,(-sin(angle)),0,cos(angle),0,0,0,0,1;
        I = my*I;
    }
    if (z==1.0){    
        Matrix4f mz; mz<<cos(angle),(-sin(angle)),0,0,sin(angle),cos(angle),0,0,0,0,1,0,0,0,0,1;
        I = mz*I;
    }
    return I;
};



// Shading
Color diffuseTerm(Color kd, Color intens, Vector4f normal, Vector4f light){
    float dotProduct = light.adjoint()*normal;
    dotProduct = max(dotProduct,(float)0);
    if (dotProduct == (float)0) {
        return Color(0.0f, 0.0f, 0.0f);
    } 
    else{
        return (kd*intens)*dotProduct;
    }
};

Color specularTerm(Color ks, Color intens, Vector4f normal, Vector4f light, Vector4f viewer, float p){
    
    Vector4f refl = light*(-1) + normal*((float)2*(light.adjoint()*normal));
    float dotProduct = refl.adjoint()*viewer;
    // cout<<dotProduct<<endl;
    dotProduct = max(dotProduct,0.0f);
    if (dotProduct==0.0f){
        return Color(0.0f,0.0f,0.0f);
    }else{
        Color result = (ks*intens)* (pow(dotProduct, p));
        if (dotProduct > 1){
            cout << dotProduct<<endl;
        }
        // if (result.intensity[0] > 1 || result.intensity[1] > 1 || result.intensity[2] > 1){
        //     printf ("why you got larger than 1?");
        // }
        return result;
    }
};

MatrixXf Find_nearest(Ray, std::vector<Surface>);


Color get_color(Vector4f viewer, Vector4f normal, Vector4f intersect){
    Color R = Color(0,0,0);
    MatrixXf None(4,2); None<<0,0,0,0,0,0,0,0;
    Ray testshadow;
    MatrixXf is_shadow;
    Color ka = (objects[hit_index].material).ka;
    Color kd = (objects[hit_index].material).kd;
    Color ks = (objects[hit_index].material).ks;
    Color em = (objects[hit_index].material).emission;
    float p = (objects[hit_index].material).sp;

    // cout << "em is " << em.intensity << endl;

    
    if (pt_light_counter != 0) {
        for (int l = 0; l < pt_light_counter; l++) {
            Vector4f pt_light_xyz = pt_xyz[l];
            Color pt_light_rgb = pt_rgb[l];
            Vector4f light = pt_light_xyz - intersect;
            light = normalize(light);
            normal = normalize(normal);
            
            
            testshadow = Ray(error_avoid(intersect,normal), light);
            is_shadow = Find_nearest(testshadow, objects);

            Vector4f is_shadow_point = Vector4f(is_shadow(4), is_shadow(5), is_shadow(6), is_shadow(7));

            
            if(is_shadow==None || (testshadow.get_t(pt_light_xyz)<testshadow.get_t(is_shadow_point))){
                Color diffuse1 = diffuseTerm(kd, pt_light_rgb, normal, light);
                
                Color specular1 = specularTerm(ks, pt_light_rgb, normal, light, viewer, p);
                // cout<<specular.intensity<<endl;

                R = R + (diffuse1 + specular1);

            }
        }
    }

    if (dl_light_counter != 0) {
    
        for (int l = 0; l < dl_light_counter; l++) {
            Color dl_light_rgb = dl_rgb[l];
//            cout<<l<<endl;
            
//            cout<<dl_light_rgb.b<<endl;
//            cout << dl_light_rgb.intensity << endl;
//            cout << "\n" << endl;
            
            Vector4f light1 = dl_xyz[l];
            
            
            light1 = normalize(light1);
//
            
            normal = normalize(normal);
            
//            cout<<light1<<endl;
             testshadow = Ray(error_avoid(intersect,normal), light1);
             is_shadow = Find_nearest(testshadow, objects);
            
             if(is_shadow==None){
//            cout<<dl_light_rgb.intensity<<endl;
                Color diffuse = diffuseTerm(kd, dl_light_rgb, normal, light1);
                
                Color specular = specularTerm(ks, dl_light_rgb, normal, light1, viewer, p);
            
//                cout<<diffuse.intensity<<endl;
                R = R + (diffuse + specular);
                 
             }
        }
    }
    return R + ka + em;
};



// PointIntersection
// in obj space
MatrixXf PointIntersection(Ray ray, Surface surface){
    Vector4f e=ray.origin;
    Vector4f d=ray.direction;

    Vector4f n,intersection,c;
    float t_1,t_2,t_3,t1,t2,t,discriminant,discriminant1,discriminant2,discriminant3,R;
    bool Flag=false;
    if(surface.name=="Sphere"){

        c = ((Sphere*) &surface) -> center;
        R = ((Sphere*) &surface) -> radius;


        discriminant1 = (d.adjoint()*(e-c))*(d.adjoint()*(e-c));
        discriminant2 = (d.adjoint()*d);
        discriminant3 = ((e-c).adjoint()*(e-c)-R*R);
        discriminant = discriminant1-discriminant2*discriminant3;

        if (discriminant<(float)0){
        }

        else if(discriminant==(float)0){
            t_1 = (d.adjoint()*(e-c)); // B/2
            t_2 = (d.adjoint()*d);     // A
            t = -t_1/t_2;
            if (t < 0){
                Flag = false;
            } else {
                Flag = true;
            }
        }
        else{

            t_1 = d.adjoint()*(e-c);
            t_2 = sqrt(discriminant);
            t_3 = (d.adjoint()*d);
            t1 = (-t_1+t_2)/t_3;
            t2 = (-t_1-t_2)/t_3;

            if (t1 < 0){
                Flag = false;
            } else if (t2 < 0){
                Flag = true;
                t = t1;
            } else {
                Flag = true;
                t = t2;
            }    
        }
    }
    
    if (Flag){
        n = ((e+d*t)-c)/R;
        intersection = e+d*t;
        MatrixXf NAndP(4,2); NAndP<<n[0],intersection[0],n[1],intersection[1],n[2],intersection[2],n[3],intersection[3];
        return NAndP;
        }
    
    else if (surface.name=="Triangle"){
        Vector4f v0 = surface.v1.coordinate;
        Vector4f v1 = surface.v2.coordinate;
        Vector4f v2 = surface.v3.coordinate;
            
        Vector4f N = cross_product((v1-v0),(v2-v0));
        float is_parallel = N.adjoint()*d;
        if (is_parallel!=0.0){
            float D = -N.adjoint()*v0;
            float tt = -(N.adjoint()*e + D) / (N.adjoint()*d);
            if (tt>0){
                Vector4f P = e+tt*d;
                Vector4f C;
                Vector4f edge0 = v1-v0;
                Vector4f pv0 = P-v0;
                C = cross_product(edge0,pv0);
                if (N.adjoint()*C>0 or N.adjoint()*C==0){
                    Vector4f edge1 = v2-v1;
                    Vector4f pv1 = P-v1;
                    C = cross_product(edge1,pv1);
                    if (N.adjoint()*C>0 or N.adjoint()*C==0){
                        Vector4f edge2 = v0-v2;
                        Vector4f pv2 = P-v2;
                        C = cross_product(edge2,pv2);
                        if (N.adjoint()*C>0 or N.adjoint()*C==0){

                            MatrixXf returnV(4,2); returnV<<N[0],P[0],N[1],P[1],N[2],P[2],0,1;
                            return returnV;
                        }
                    }
                }
            }   
        }
    }

    else if (surface.name=="TriangleNormal"){
        Vector4f v0 = surface.v1.coordinate;
        Vector4f v1 = surface.v2.coordinate;
        Vector4f v2 = surface.v3.coordinate;
            
        Vector4f N = cross_product((v1-v0),(v2-v0));
        float is_parallel = N.adjoint()*d;
        if (is_parallel!=0.0){
            float D = -N.adjoint()*v0;
            float tt = -(N.adjoint()*e + D) / (N.adjoint()*d);
            if (tt>0){
                Vector4f P = e+tt*d;
                Vector4f C;
                Vector4f edge0 = v1-v0;
                Vector4f pv0 = P-v0;
                C = cross_product(edge0,pv0);
                if (N.adjoint()*C>0 or N.adjoint()*C==0){
                    Vector4f edge1 = v2-v1;
                    Vector4f pv1 = P-v1;
                    C = cross_product(edge1,pv1);
                    if (N.adjoint()*C>0 or N.adjoint()*C==0){
                        Vector4f edge2 = v0-v2;
                        Vector4f pv2 = P-v2;
                        C = cross_product(edge2,pv2);
                        if (N.adjoint()*C>0 or N.adjoint()*C==0){
                            Vector4f point = Vector4f(P[0], P[1], P[2], 1);
                            Vector4f norm = surface.getNormal(point);
                            MatrixXf returnV(4,2); returnV<<norm[0],P[0],norm[1],P[1],norm[2],P[2],0,1;
                            return returnV;
                        }
                    }
                }
            }   
        }
    }
        
        
        
    
    MatrixXf None(4,2); None<<0,0,0,0,0,0,0,0;
    return None;
};

MatrixXf Find_nearest(Ray ray, std::vector<Surface> surface){
    float t;
    float compare=100000000;
    bool Flag = false;
    Vector4f finalpoint, finalnormal;
    MatrixXf None(4,2); None<<0,0,0,0,0,0,0,0;
    
    Vector4f returnP, returnN;

    for (int i=0; i < obj_counter; i++){

        // if (surface[i].transformation.inverse() != surface[i].transformation){
        //     cout << surface[i].transformation.inverse() << endl;
        // }

        Vector4f origin = (surface[i].transformation).inverse()* ray.origin;
        
        Vector4f direction = (surface[i].transformation).inverse()* ray.direction;
        
        direction = normalize(direction);
        
        
        
        
        Ray newRay = Ray(origin,direction);

        MatrixXf intersection = PointIntersection(newRay,surface[i]);

        if (intersection!=None){
            // Flag = true;
            Vector4f point(intersection(4),intersection(5),intersection(6),intersection(7));
            Vector4f normal(intersection(0),intersection(1),intersection(2),intersection(3));
            finalpoint = surface[i].transformation* point;
            finalnormal = surface[i].transformation.inverse().transpose()* normal;
            finalnormal = normalize(finalnormal);
            t = ray.get_t(finalpoint);
            if (t>0 && t<compare){
                Flag = true;
                hit_index = i;
                returnP = finalpoint;
                returnN = finalnormal;
                compare = t;
            }
        }
    }
    if (Flag){
        MatrixXf returnValue(4,2); returnValue<<returnN[0],returnP[0],returnN[1], returnP[1],returnN[2],returnP[2],returnN[3], returnP[3];
        return returnValue;
    }
    return None;
};

Vector4f find_reflection(Ray ray, Vector4f normal){
    
    // if (ray.direction[3] != 0 || ray.origin[3] != 1){
    //     cout << ray.direction << endl;
    //     cout << ray.origin << endl;
    // }
    
    Vector4f direction = ray.direction;
    float c = -normal.adjoint()*direction;
    float divider = normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2];
    c = c/divider;
    Vector4f Reflect = direction + (2*normal*c);
    
    
    return Reflect;
};

Color trace(Ray ray,int TraceDepth){
	Color R = Color(0.0, 0.0, 0.0);
    MatrixXf None(4,2); None<<0,0,0,0,0,0,0,0;
    if (TraceDepth < 0){
        return R;
    }

    // find the neares hit
    MatrixXf result = Find_nearest(ray, objects);

    Color ks = objects[hit_index].material.ks;
    if (result == None){
        return skycolor;
    }

    // printf("hit something\n\n\n\n\n");

    Vector4f normal = Vector4f(result(0),result(1),result(2),result(3));
    Vector4f intersect = Vector4f(result(4),result(5),result(6),result(7));

   	// cout << "normal is \n" << normal << endl;
   	// cout << "intersect is \n" << intersect << endl;
   	// cout << "\n\n";
    
    // generate another ray
    normal = normalize(normal);
    Vector4f rflct = find_reflection(ray, normal);
    rflct = normalize(rflct);
    
    Ray reflection_ray = Ray(error_avoid(intersect,normal), rflct);
    

    Vector4f viewer = ray.origin-intersect;
    viewer = normalize(viewer);

    Color pp = get_color(viewer, normal, intersect);
    // cout << pp.intensity << endl;
    R = R + pp;
    // cout << R.intensity << endl;
    
    hit_index = -1;
    R = R + ks*trace(reflection_ray,TraceDepth-1);
    return R;
};

//***********************************************
//  Input Parser
//***********************************************
void loadScene(std::string file) {
  //store variables and set stuff at the end
    fname = "output.png";
    std::ifstream inpfile(file.c_str());
    if(!inpfile.is_open()) {
        std::cout << "Unable to open file" << std::endl;
    } else {
        std::string line;
        //MatrixStack mst;
        while(inpfile.good()) {
            std::vector<std::string> splitline;
            std::string buf;
            std::getline(inpfile,line);
            std::stringstream ss(line);
            while (ss >> buf) {
                splitline.push_back(buf);
            }
            //Ignore blank lines
            if(splitline.size() == 0) {
                continue;
            }
            //Ignore comments
            if(splitline[0][0] == '#') {
                continue;
            }
            //Valid commands:
            //size width height
            //  must be first command of file, controls image size
            else if(!splitline[0].compare("size")) {
                WIDTH = atoi(splitline[1].c_str());
                HEIGHT = atoi(splitline[2].c_str());
            }
            //maxdepth depth
            //  max # of bounces for ray (default 5)
            else if(!splitline[0].compare("maxdepth")) {
                TraceDepth = atoi(splitline[1].c_str());
            }
            //output filename
            //  output file to write image to 
            else if(!splitline[0].compare("output")) {
                fname = splitline[1];
            }
            //camera lookfromx lookfromy lookfromz lookatx lookaty lookatz upx upy upz fov
            //  speciﬁes the camera in the standard way, as in homework 2.
            else if(!splitline[0].compare("camera")) {
                // lookfrom:
                    float camara_x = atof(splitline[1].c_str());
                    float camara_y = atof(splitline[2].c_str());
                    float camara_z = atof(splitline[3].c_str());

                    camara_position = Vector4f(camara_x, camara_y, camara_z, 1.0f);
                // lookat:
                    float look_x = atof(splitline[4].c_str());
                    float look_y = atof(splitline[5].c_str());
                    float look_z = atof(splitline[6].c_str());

                    camara_looking_direction = normalize(Vector4f(look_x-camara_x, look_y-camara_y , look_z-camara_z, 0.0f));
                // up:
                    float up_x = atof(splitline[7].c_str());
                    float up_y = atof(splitline[8].c_str());
                    float up_z = atof(splitline[9].c_str());

                    camara_up_direction = Vector4f(up_x, up_y, up_z, 0.0f);
                // fov:
                    fov = atof(splitline[10].c_str());
            }
            //sphere x y z radius
            //  Deﬁnes a sphere with a given position and radius.
            else if(!splitline[0].compare("sphere")) {
                float sphere_x = atof(splitline[1].c_str());
                float sphere_y = atof(splitline[2].c_str());
                float sphere_z = atof(splitline[3].c_str());
                float radius = atof(splitline[4].c_str());

                // Create new sphere:
                //   Store 4 numbers
                //   Store current property values
                //   Store current top of matrix stack

                Vector4f center = Vector4f(sphere_x, sphere_y, sphere_z, 1.0f);
                Material material = Material(Ambient_input, Diffuse_input, Specular_input, Shiness_input, reflection_input, emission_input);
                Matrix4f transformation_value = tfmStack.get_transformation();

                Sphere new_sphere = Sphere(transformation_value, center, radius, material);
                objects.push_back(new_sphere);
            }
            //maxverts number
            //  Deﬁnes a maximum number of vertices for later triangle speciﬁcations. 
            //  It must be set before vertices are deﬁned.
            else if(!splitline[0].compare("maxverts")) {
                // Care if you want
                // Here, either declare array size
                // Or you can just use a STL vector, in which case you can ignore this
                continue;
            }
            //maxvertnorms number
            //  Deﬁnes a maximum number of vertices with normals for later speciﬁcations.
            //  It must be set before vertices with normals are deﬁned.
            else if(!splitline[0].compare("maxvertnorms")) {
                // Care if you want
                continue;
            }
            //vertex x y z
            //  Deﬁnes a vertex at the given location.
            //  The vertex is put into a pile, starting to be numbered at 0.
            else if(!splitline[0].compare("vertex")) {
                float v_x = atof(splitline[1].c_str());
                float v_y = atof(splitline[2].c_str());
                float v_z = atof(splitline[3].c_str());
                // Create a new vertex with these 3 values, store in some array
                Vector4f point = Vector4f(v_x, v_y, v_z, 1.0f);
                Vertex v = Vertex(point);
                vertices.push_back(v);
            }
            //vertexnormal x y z nx ny nz
            //  Similar to the above, but deﬁne a surface normal with each vertex.
            //  The vertex and vertexnormal set of vertices are completely independent
            //  (as are maxverts and maxvertnorms).
            else if(!splitline[0].compare("vertexnormal")) {
                float x = atof(splitline[1].c_str());
                float y = atof(splitline[2].c_str());
                float z = atof(splitline[3].c_str());
                float nx = atof(splitline[4].c_str());
                float ny = atof(splitline[5].c_str());
                float nz = atof(splitline[6].c_str());
                
                Vector4f point = Vector4f(x, y, z, 1.0f);
                Vector4f normal = Vector4f(nx, ny, nz, 0.0f);
                VertexNormal vn = VertexNormal(point, normal);
                verticesNormal.push_back(vn);
                // Create a new vertex+normal with these 6 values, store in some array
            }
            // tri v1 v2 v3
            //  Create a triangle out of the vertices involved (which have previously been speciﬁed with
            //  the vertex command). The vertices are assumed to be speciﬁed in counter-clockwise order. Your code
            //  should internally compute a face normal for this triangle.
            else if(!splitline[0].compare("tri")) {
                int v1 = atof(splitline[1].c_str());
                int v2 = atof(splitline[2].c_str());
                int v3 = atof(splitline[3].c_str());

                // Create new triangle:
                //   Store pointer to array of vertices
                //   Store 3 integers to index into array
                //   Store current property values
                //   Store current top of matrix stack
                Vertex vtx1 = vertices[v1];
                Vertex vtx2 = vertices[v2];
                Vertex vtx3 = vertices[v3];

                Material material = Material(Ambient_input, Diffuse_input, Specular_input, Shiness_input, reflection_input, emission_input);
                Matrix4f transformation_value = tfmStack.get_transformation();

                triangle tri = triangle(vtx1, vtx2, vtx3, transformation_value, material);
                objects.push_back(tri);
            }
            //trinormal v1 v2 v3
            //  Same as above but for vertices speciﬁed with normals.
            //  In this case, each vertex has an associated normal, 
            //  and when doing shading, you should interpolate the normals 
            //  for intermediate points on the triangle.
            else if(!splitline[0].compare("trinormal")) {
                int v1 = atof(splitline[1].c_str());
                int v2 = atof(splitline[2].c_str());
                int v3 = atof(splitline[3].c_str());
                // Create new triangle:
                //   Store pointer to array of vertices (Different array than above)
                //   Store 3 integers to index into array
                //   Store current property values
                //   Store current top of matrix stack
                VertexNormal vtxn1 = verticesNormal[v1];
                VertexNormal vtxn2 = verticesNormal[v2];
                VertexNormal vtxn3 = verticesNormal[v3];

                Material material = Material(Ambient_input, Diffuse_input, Specular_input, Shiness_input, reflection_input, emission_input);
                Matrix4f transformation_value = tfmStack.get_transformation();

                triangleNormal trin = triangleNormal(vtxn1, vtxn2, vtxn3, transformation_value, material);
                objects.push_back(trin);

            }
            //translate x y z
            //  A translation 3-vector
            else if(!splitline[0].compare("translate")) {
                float x = atof(splitline[1].c_str());
                float y = atof(splitline[2].c_str());
                float z = atof(splitline[3].c_str());
                // Update top of matrix stack
                Matrix4f transformation = translation(Vector3f(x, y, z));

                tfmStack.push_transformation(transformation);

                cout << "now translate" << endl;
                cout << tfmStack.get_transformation() << endl;
                cout << "\n";

            }
            //rotate x y z angle
            //  Rotate by angle (in degrees) about the given axis as in OpenGL.
            else if(!splitline[0].compare("rotate")) {
                float x = atof(splitline[1].c_str());
                float y = atof(splitline[2].c_str());
                float z = atof(splitline[3].c_str());
                float angle = atof(splitline[4].c_str());
                // Update top of matrix stack
                Matrix4f transformation = rotate(Vector4f(x, y, z, angle));
 
                tfmStack.push_transformation(transformation);  
                cout << "now rotate" << endl;
                cout << tfmStack.get_transformation() << endl;
                cout << "\n";

            }
            //scale x y z
            //  Scale by the corresponding amount in each axis (a non-uniform scaling).
            else if(!splitline[0].compare("scale")) {
                float x = atof(splitline[1].c_str());
                float y = atof(splitline[2].c_str());
                float z = atof(splitline[3].c_str());
                //  date top of matrix stack
                Matrix4f transformation = scale(Vector3f(x, y, z));
                
                tfmStack.push_transformation(transformation);
                cout << "now scale" << endl;
                cout << tfmStack.get_transformation() << endl;
                cout << "\n";

            }
            //pushTransform
            //  Push the current modeling transform on the stack as in OpenGL. 
            //  You might want to do pushTransform immediately after setting 
            //   the camera to preserve the “identity” transformation.
            else if(!splitline[0].compare("pushTransform")) {
                tfmStack.push_transformation();
                cout << "now pushTransform" << endl;
                cout << tfmStack.get_transformation() << endl;
                cout << "\n";

            }
            //popTransform
            //  Pop the current transform from the stack as in OpenGL. 
            //  The sequence of popTransform and pushTransform can be used if 
            //  desired before every primitive to reset the transformation 
            //  (assuming the initial camera transformation is on the stack as 
            //  discussed above).
            else if(!splitline[0].compare("popTransform")) {
                tfmStack.pop_transformation();
                cout << "now popTransform" << endl;
                cout << tfmStack.get_transformation() << endl;
                cout << "\n";

            }
            //directional x y z r g b
            //  The direction to the light source, and the color, as in OpenGL.
            else if(!splitline[0].compare("directional")) {
                float x = atof(splitline[1].c_str());
                float y = atof(splitline[2].c_str());
                float z = atof(splitline[3].c_str());
                float r = atof(splitline[4].c_str());
                float g = atof(splitline[5].c_str());
                float b = atof(splitline[6].c_str());
                // add light to scene...
                Vector4f directional_xyz = Vector4f(x, y, z, 0.0f);
                Color directional_rgb = Color(r, g, b);
                dl_xyz.push_back(directional_xyz);
                dl_rgb.push_back(directional_rgb);

            }
            //point x y z r g b
            //  The location of a point source and the color, as in OpenGL.
            else if(!splitline[0].compare("point")) {
                float x = atof(splitline[1].c_str());
                float y = atof(splitline[2].c_str());
                float z = atof(splitline[3].c_str());
                float r = atof(splitline[4].c_str());
                float g = atof(splitline[5].c_str());
                float b = atof(splitline[6].c_str());
                // add light to scene...
                Vector4f point_xyz = Vector4f(x, y, z, 1.0f);
                Color point_rgb = Color(r, g, b);
                pt_xyz.push_back(point_xyz);
                pt_rgb.push_back(point_rgb);

            }
            //attenuation const linear quadratic
            //  Sets the constant, linear and quadratic attenuations 
            //  (default 1,0,0) as in OpenGL.
            else if(!splitline[0].compare("attenuation")) {
                // const: atof(splitline[1].c_str())
                // linear: atof(splitline[2].c_str())
                // quadratic: atof(splitline[3].c_str())
            }
            //ambient r g b
            //  The global ambient color to be added for each object 
            //  (default is .2,.2,.2)
            else if(!splitline[0].compare("ambient")) {
                float r = atof(splitline[1].c_str());
                float g = atof(splitline[2].c_str());
                float b = atof(splitline[3].c_str());

                Ambient_input = Color(r, g, b);
            }
            //diﬀuse r g b
            //  speciﬁes the diﬀuse color of the surface.
            else if(!splitline[0].compare("diffuse")) {
                float r = atof(splitline[1].c_str());
                float g = atof(splitline[2].c_str());
                float b = atof(splitline[3].c_str());

                // Update current properties
                Diffuse_input = Color(r, g, b);
            }
            //specular r g b 
            //  speciﬁes the specular color of the surface.
            else if(!splitline[0].compare("specular")) {
                float r = atof(splitline[1].c_str());
                float g = atof(splitline[2].c_str());
                float b = atof(splitline[3].c_str());
                
                // Update current properties
                Specular_input = Color(r, g, b);
            }
            //shininess s
            //  speciﬁes the shininess of the surface.
            else if(!splitline[0].compare("shininess")) {
                // Update current properties
                Shiness_input = atof(splitline[1].c_str());
            }
            //emission r g b
            //  gives the emissive color of the surface.
            else if(!splitline[0].compare("emission")) {
                float r = atof(splitline[1].c_str());
                float g = atof(splitline[2].c_str());
                float b = atof(splitline[3].c_str());

                // Update current properties
                emission_input = Color(r, g, b);
            } else {
                std::cerr << "Unknown command: " << splitline[0] << std::endl;
            }
        }
        inpfile.close();
    }
};

//***********************************************
//  Main Funciton
//***********************************************
int main(int args, char* argv[]){

    char* inputname;

    for (int i = 1; i < args; ++i){
        if(string(argv[i]) == "-input"){
            inputname = argv[i+1];
        }
    }







    loadScene(inputname);

    // specify camara:
    // WIDTH = 480;
    // HEIGHT = 480;




    aspect_ratio = (float)WIDTH/(float)HEIGHT;

    camara_right_direction = cross_product(camara_looking_direction, camara_up_direction);
    camara_right_direction = normalize(camara_right_direction);
    camara_up_direction = cross_product(camara_right_direction, camara_looking_direction);
    camara_up_direction = normalize(camara_up_direction);
    // TraceDepth = 5;
    
    float fovV = fov/180.0f*PI;
    float vertical_offset = tan(fovV/2);
    
    float horizontal_offset = vertical_offset*aspect_ratio;
    
    float rr = horizontal_offset;
    float ll = -horizontal_offset;
    float tt = vertical_offset;
    float bb = -vertical_offset;

    obj_counter = objects.size();
    maxvertex = vertices.size();
    maxvertexnormal = verticesNormal.size();
    pt_light_counter = pt_xyz.size(); 
    dl_light_counter = dl_xyz.size();


    cout << "WIDTH is " << WIDTH << endl;
    // HEIGHT = 200;
    cout << "HEIGHT is " << HEIGHT << endl;
    // aspect_ratio = WIDTH/HEIGHT;
    cout << "aspect_ratio is " << aspect_ratio << endl;
    // fov = 36.8f;
    cout << "fov is " << fov << endl;
    // camara_position = Vector4f(0.0f, 0.0f, 0.0f, 1.0f);
    cout << "camara_position is " << camara_position << endl;
    // camara_looking_direction = Vector4f (0.0f, 0.0f, -1.0f, 0.0f);
    cout << "camara_looking_direction is " << camara_looking_direction << endl;
    // camara_up_direction = Vector4f (0.0f, 1.0f, 0.0f, 0.0f);
    cout << "camara_up_direction is " << camara_up_direction << endl;
    // camara_right_direction = cross_product(camara_looking_direction, camara_up_direction);
    cout << "camara_right_direction is " << camara_right_direction << endl;
    // TraceDepth = 5;

    // obj_counter = objects.size();
    cout << "obj_counter is " << obj_counter << endl;
    // pt_light_counter = pt_xyz.size();
    cout << "dl_light_counter is " << dl_light_counter << endl;
    // dl_light_counter = dl_xyz.size();
    cout << "pt_light_counter is " << pt_light_counter << endl;
    // maxvertex = vertices.size();
    cout << "vertices number is  " << maxvertex << endl;
    // maxvertexnormal = verticesNormal.size();
    cout << "normal_vertecies number is " << maxvertexnormal << endl;
    // float fovV = fov/180.0f*PI;
    cout << "fovV is " << fovV << endl;
    // float vertical_offset = tan(fovV/2);
    cout << "vertical_offset is " << vertical_offset << endl;
    
    // float horizontal_offset = vertical_offset*aspect_ratio;
    cout << "horizontal_offset is " << horizontal_offset << endl;

    cout << "Ambient_input is " << Ambient_input.intensity << endl;
    // Diffuse_input = Color(1, 0, 0);
    cout << "Diffuse_input is " << Diffuse_input.intensity << endl;
    // Specular_input = Color(0, 0, 0);
    cout << "Specular_input is " << Specular_input.intensity << endl;
    // Shiness_input = 0;
    cout << "Shiness_input is " << Shiness_input << endl;


    cout << "Material is " << objects[0].material.ka.intensity << endl;


    FreeImage_Initialise();

    FIBITMAP *bitmap = FreeImage_Allocate(WIDTH, HEIGHT, BPP);
    RGBQUAD color;

    if (!bitmap)
        exit(1);

    for (int i=0; i<WIDTH; i++){
        for (int j=0; j<HEIGHT; j++){

        	// Ray Generation according to camara geometry
            if (i == 0 &&  j == 1){
                cout << "complete 1 ray"<<endl;
            }

            if (i == WIDTH/4 && j == HEIGHT/4){
                cout << "12.5 percent complete"<<endl;
            }

            if (i == WIDTH/4 && j == HEIGHT/2){
                cout << "25 percent complete"<<endl;
            }


            if (i == WIDTH/2 && j == HEIGHT/2){
                cout << "50 percent complete"<<endl;
            }

            if (i == WIDTH/2 && j == HEIGHT){
                cout << "62.5 percent complete"<<endl;
            }

            if (i == WIDTH && j == HEIGHT/4){
                cout << "75 percent complete"<<endl;
            }

            if (i == WIDTH && j == HEIGHT/2){
                cout << "87.5 percent complete"<<endl;
            }


        	float u = ll + (rr-ll)*(i+0.5)/WIDTH;
        	float v = bb + (tt-bb)*(j+0.5)/HEIGHT;

        	Vector4f direction = camara_looking_direction + u * camara_right_direction + v * camara_up_direction;
        	direction = normalize(direction);

        	Vector4f origin = camara_position;

        	Ray initial_ray = Ray (origin, direction);

        	// cout << "initial ray origin\n" << initial_ray.origin << endl;
       		// cout << "initial ray direction\n" << initial_ray.direction << endl;
       		// cout << "\n\n";

        	Color result = trace(initial_ray, TraceDepth);

        	color.rgbRed = (result.intensity[2]*255 > 255 ? 255 : result.intensity[2]*255);
        	color.rgbGreen = result.intensity[1]*255 > 255 ? 255 : result.intensity[1]*255;
       		color.rgbBlue = result.intensity[0]*255 > 255 ? 255 : result.intensity[0]*255;

			FreeImage_SetPixelColor (bitmap, i, j, &color);
        }
    }

    const char * filename = fname.c_str();
    if (FreeImage_Save(FIF_PNG, bitmap, filename, 0))
        cout << "Image successfully saved!" << endl;
    
    FreeImage_DeInitialise();
}


       




