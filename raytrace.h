///////////////////////////////////////////////////////////////////////
// A framework for a raytracer.
////////////////////////////////////////////////////////////////////////
#pragma once
#include "acceleration.h"

class Shape;
class AccelerationBvh;

const float PI = 3.14159f;
const float Radians = PI/180.0f;    // Convert degrees to radians

////////////////////////////////////////////////////////////////////////
// Texture: encapsulates a BRDF and communication with a shader.
////////////////////////////////////////////////////////////////////////
class Texture
{
public:
    unsigned int id;
    int width, height, depth;
    unsigned char* image;
    Texture(const std::string &path);
};

////////////////////////////////////////////////////////////////////////
// Material: encapsulates a BRDF and communication with a shader.
////////////////////////////////////////////////////////////////////////
class Material
{
 public:
    vec3 Kd, Ks, Kt;
    float alpha;
    float ior;
    Texture* tex;

    virtual bool isLight() { return false; }

    Material()  : Kd(vec3(1.0, 0.5, 0.0)), Ks(vec3(1,1,1)), alpha(1.0), tex(NULL), Kt(0, 0, 0), ior(1) {}
    Material(const vec3 d, const vec3 s, const float a, const vec3 t, const float ior) 
        : Kd(d), Ks(s), alpha(a), tex(NULL), Kt(t), ior(ior) {}
    Material(const Material& o) { Kd = o.Kd;  Ks = o.Ks;  alpha = o.alpha;  tex = o.tex; Kt = o.Kt; ior = o.ior; }

    void setTexture(const std::string path) { tex = new Texture(path); };
};

////////////////////////////////////////////////////////////////////////
// Data structures for storing meshes -- mostly used for model files
// read in via ASSIMP.
//
// A MeshData holds two lists (stl::vector) one for vertices
// (VertexData: consisting of point, normal, texture, and tangent
// vectors), and one for triangles (ivec3: consisting of three
// indices into the vertex array).
    
class VertexData
{
 public:
    vec3 pnt;
    vec3 nrm;
    vec2 tex;
    vec3 tan;
    VertexData(const vec3& p, const vec3& n, const vec2& t, const vec3& a) 
        : pnt(p), nrm(n), tex(t), tan(a) 
    {}
};

struct MeshData
{
    std::vector<VertexData> vertices;
    std::vector<ivec3> triangles;
    Material *mat;
};

////////////////////////////////////////////////////////////////////////
// Light: encapsulates a light and communiction with a shader.
////////////////////////////////////////////////////////////////////////
class Light: public Material
{
public:

    Light(const vec3 e) : Material() { Kd = e; }
    virtual bool isLight() { return true; }
    //virtual void apply(const unsigned int program);
};

////////////////////////////////////////////////////////////////////////////////
// Scene


struct Camera {
    vec3 eye;
    quat orient;
    vec3 ambient;
    float ry;
    float front, back;
    float spin, tilt;
    float cDist;
};

class Scene {
public:
    int width, height;
    //Realtime* realtime;         // Remove this (realtime stuff)
    std::vector<Shape*> shapes;
    std::vector<Shape*> lights;
    AccelerationBvh* bvh;
    Camera camera;
    Material* currentMat;

    Scene();
    void Finit();

    // The scene reader-parser will call the Command method with the
    // contents of each line in the scene file.
    void Command(const std::vector<std::string>& strings,
                 const std::vector<float>& f);

    // To read a model file into the scene via ASSIMP, call ReadAssimpFile.  
    void ReadAssimpFile(const std::string& path, const mat4& M);

    // Once ReadAssimpFile parses the information from the model file,
    // it will call:
    void triangleMesh(MeshData* mesh);

    // The main program will call the TraceImage method to generate
    // and return the image.  This is the Ray Tracer!
    void TraceImage(Color* image, const int pass);

    Color TracePath(const Ray& _ray);

    Intersection SampleLight();
    float PdfLight(Shape* L);
    Color EvalRadiance(Shape* L);

    vec3 SampleBrdf(Shape* shape, vec3 N, vec3 w);


    float GeometryFactory(const Intersection& A, const Intersection& B);
    vec3 SampleLobe(vec3 A, float c, float theta);
    Intersection SampleSphere(Shape* _shape);

    void sphere(const vec3 _center, const float _r, Material* _mat);
    void box(const vec3 _base, const vec3 _diag, Material* _mat);
    void cylinder(const vec3 _base, const vec3 _axis, const float _radius, Material* _mat);
};
