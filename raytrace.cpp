//////////////////////////////////////////////////////////////////////
// Provides the framework for a raytracer.
////////////////////////////////////////////////////////////////////////

#include <vector>

#ifdef _WIN32
// Includes for Windows
#include <crtdbg.h>
#include <cstdlib>
#include <limits>
#include <windows.h>
#else
// Includes for Linux
#endif

#include "acceleration.h"
#include "geom.h"
#include "raytrace.h"


#define STB_IMAGE_IMPLEMENTATION
#define STBI_FAILURE_USERMSG
#include "stb_image.h"

// A good quality *thread-safe* Mersenne Twister random number generator.
#include <random>
#include <fstream>
std::random_device device;
std::mt19937_64 RNGen(device());
std::uniform_real_distribution<float> myrandom(0.0, 1.0);
// Call myrandom(RNGen) to get a uniformly distributed random number in [0,1].

Scene::Scene()
{
    // realtime = new Realtime();
}

void Scene::Finit()
{
    bvh = new AccelerationBvh(shapes);
}

void Scene::triangleMesh(MeshData* mesh)
{
    for (auto& tri : mesh->triangles) {
        Shape* triangle = new Triangle(
            mesh->vertices[tri[0]].pnt,
            mesh->vertices[tri[1]].pnt,
            mesh->vertices[tri[2]].pnt,
            mesh->vertices[tri[0]].nrm,
            mesh->vertices[tri[1]].nrm,
            mesh->vertices[tri[2]].nrm);
        triangle->set_material(currentMat);
        shapes.push_back(triangle);
        if (currentMat->isLight())
            lights.push_back(shapes[shapes.size() - 1]);
    }
}

Texture::Texture(const std::string& bpath)
    : id(0)
{
    // Replace backslashes with forward slashes -- Good for Linux, and maybe Windows?
    std::string path = bpath;
    std::string bs = "\\";
    std::string fs = "/";
    while (path.find(bs) != std::string::npos) {
        path.replace(path.find(bs), 1, fs);
    }

    // Does the file exist?
    std::ifstream find_it(path.c_str());
    if (find_it.fail()) {
        std::cerr << "Texture file not found: " << path << std::endl;
        exit(-1);
    } else {
        // Read image, and check for success
        stbi_set_flip_vertically_on_load(true);
        image = stbi_load(path.c_str(), &width, &height, &depth, 4);
        printf("%d %d %d %s\n", depth, width, height, path.c_str());
        if (!image) {
            printf("\nRead error on file %s:\n  %s\n\n", path.c_str(), stbi_failure_reason());
            exit(-1);
        }
    }
}

quat Orientation(int i,
    const std::vector<std::string>& strings,
    const std::vector<float>& f)
{
    quat q(1, 0, 0, 0); // Unit quaternion
    while (i < strings.size()) {
        std::string c = strings[i++];
        if (c == "x")
            q *= angleAxis(f[i++] * Radians, Xaxis());
        else if (c == "y")
            q *= angleAxis(f[i++] * Radians, Yaxis());
        else if (c == "z")
            q *= angleAxis(f[i++] * Radians, Zaxis());
        else if (c == "q") {
            q *= quat(f[i + 0], f[i + 1], f[i + 2], f[i + 3]);
            i += 4;
        } else if (c == "a") {
            q *= angleAxis(f[i + 0] * Radians, normalize(vec3(f[i + 1], f[i + 2], f[i + 3])));
            i += 4;
        }
    }
    return q;
}

void Scene::Command(const std::vector<std::string>& strings,
    const std::vector<float>& f)
{
    if (strings.size() == 0)
        return;
    std::string c = strings[0];

    if (c == "screen") {
        // syntax: screen width height
        // realtime->setScreen(int(f[1]), int(f[2]));

        width = int(f[1]);
        height = int(f[2]);
    }

    else if (c == "camera") {
        // syntax: camera x y z   ry   <orientation spec>
        // Eye position (x,y,z),  view orientation (qw qx qy qz),  frustum height ratio ry
        camera.eye = vec3(f[1], f[2], f[3]);
        camera.orient = Orientation(5, strings, f);
        camera.ry = f[4];
        // realtime->setCamera(vec3(f[1], f[2], f[3]), Orientation(5, strings, f), f[4]);
    }

    else if (c == "ambient") {
        // syntax: ambient r g b
        // Sets the ambient color.  Note: This parameter is temporary.
        // It will be ignored once your raytracer becomes capable of
        // accurately *calculating* the true ambient light.
        camera.ambient = vec3(f[1], f[2], f[3]);
        // realtime->setAmbient(vec3(f[1], f[2], f[3]));
    }

    else if (c == "brdf") {
        // syntax: brdf  r g b   r g b  alpha
        // later:  brdf  r g b   r g b  alpha  r g b ior
        // First rgb is Diffuse reflection, second is specular reflection.
        // third is beer's law transmission followed by index of refraction.
        // Creates a Material instance to be picked up by successive shapes
        currentMat = new Material(vec3(f[1], f[2], f[3]), vec3(f[4], f[5], f[6]), f[7], vec3(f[8], f[9], f[10]), f[11]);
    }

    else if (c == "light") {
        // syntax: light  r g b
        // The rgb is the emission of the light
        // Creates a Material instance to be picked up by successive shapes
        currentMat = new Light(vec3(f[1], f[2], f[3]));
    }

    else if (c == "sphere") {
        // syntax: sphere x y z   r
        // Creates a Shape instance for a sphere defined by a center and radius
        sphere(vec3(f[1], f[2], f[3]), f[4], currentMat);
        if (currentMat->isLight())
            lights.push_back(shapes[shapes.size() - 1]);
    }

    else if (c == "box") {
        // syntax: box bx by bz   dx dy dz
        // Creates a Shape instance for a box defined by a corner point and diagonal vector
        box(vec3(f[1], f[2], f[3]), vec3(f[4], f[5], f[6]), currentMat);
        if (currentMat->isLight())
            lights.push_back(shapes[shapes.size() - 1]);
    }

    else if (c == "cylinder") {
        // syntax: cylinder bx by bz   ax ay az  r
        // Creates a Shape instance for a cylinder defined by a base point, axis vector, and radius
        cylinder(vec3(f[1], f[2], f[3]), vec3(f[4], f[5], f[6]), f[7], currentMat);
        if (currentMat->isLight())
            lights.push_back(shapes[shapes.size() - 1]);
    }

    else if (c == "mesh") {
        // syntax: mesh   filename   tx ty tz   s   <orientation>
        // Creates many Shape instances (one per triangle) by reading
        // model(s) from filename. All triangles are rotated by a
        // quaternion (qw qx qy qz), uniformly scaled by s, and
        // translated by (tx ty tz) .
        mat4 modelTr = translate(vec3(f[2], f[3], f[4]))
            * scale(vec3(f[5], f[5], f[5]))
            * toMat4(Orientation(6, strings, f));
        ReadAssimpFile(strings[1], modelTr);
    }

    else {
        fprintf(stderr, "\n*********************************************\n");
        fprintf(stderr, "* Unknown command: %s\n", c.c_str());
        fprintf(stderr, "*********************************************\n\n");
    }
}

void Scene::TraceImage(Color* image, const int pass)
{
    if (image == nullptr)
        return;
    const float rx = camera.ry * (float)width / (float)height;
    const vec3 right = rx * transformVector(camera.orient, { 1, 0, 0 });
    const vec3 up = camera.ry * transformVector(camera.orient, { 0, 1, 0 });
    const vec3 forward = transformVector(camera.orient, { 0, 0, 1 });

#pragma omp parallel for schedule(dynamic, 1) // Magic: Multi-thread y loop
    for (int y = 0; y < height; y++) {

        // fprintf(stderr, "Rendering %4d\r", y);
        for (int x = 0; x < width; x++) {
            const float dx = 2.0f * (x + myrandom(RNGen)) / width - 1;
            const float dy = 2.0f * (y + myrandom(RNGen)) / height - 1;

            const Ray ray(camera.eye, glm::normalize(dx * right + dy * up - forward));

           
            Color c = TracePath(ray);
            
            constexpr float stuff = 30;
            if (c.x > -stuff && c.x <= stuff && c.y > -stuff && c.y <= stuff && c.z > -stuff && c.z <= stuff)
                if (!isnan(c.x) && !isinf(c.x) && !isnan(c.y) && !isinf(c.y) && !isnan(c.z) && !isinf(c.z))
                    image[y * width + x] += c;// glm::clamp(c, -15.f, 15.f);
            //image[y * width + x] /= (pass == 0) ? 1 : 2;
        }
    }
    //fprintf(stderr, "\n");
}

Color Scene::TracePath(const Ray& _ray)
{
    Color C = { 0, 0, 0 };
    Color W = { 1, 1, 1 };

    Intersection P{};

    P = bvh->intersect(_ray);

    if (P.object) {
        if (P.object->get_material()->isLight()) {
            return EvalRadiance(P.object);
        }

        vec3 wo = -normalize(_ray.d);
        while (myrandom(RNGen) < 0.8f) {

            P.normal = normalize(P.normal);
            const Intersection L = SampleLight();
            const vec3 N = P.normal;
            float p = PdfLight(L.object) / GeometryFactory(P, L);
            vec3 wi = normalize(L.point - P.point);
            float q = P.object->PdfBrdf(wo, P.normal, wi) * 0.8f;
            float wmis = (p * p) / (p * p + q * q);
            Ray r(P.point, wi);
            if (p > 10e-6f) {
                const Intersection I = bvh->intersect(r);
                const vec3 d = glm::abs(I.point - L.point);
                if (I.object && d.x < 10e-3f && d.y < 10e-3f && d.z < 10e-3f) {
                    const vec3 f = P.object->EvalScattering(wi, N, wo, 0);
                    const vec3 radiance = EvalRadiance(L.object);
                    C += W * wmis * (f / p) * radiance;
                }
            }

            wi = normalize(SampleBrdf(P.object, N, wo));
            r = Ray(P.point, wi);

            Intersection Q = bvh->intersect(r);
            if (!Q.object)
                break;

            const vec3 f = P.object->EvalScattering(wi, P.normal, wo, Q.t);


            p = P.object->PdfBrdf(wo, P.normal, wi) * 0.8f;
            if (p < 10e-6)
                break;

            W *= f / p;

            if (Q.object->get_material()->isLight()) {
                q = PdfLight(Q.object) / GeometryFactory(P,Q);
                wmis = (p * p) / (p * p + q * q);
                const vec3 radiance = EvalRadiance(Q.object);
                C += W * wmis * radiance;
                break;
            }

            P = Q;
            wo = -wi;
        }
    }

    return C;
}

Intersection Scene::SampleLight()
{
    const int r = (int)roundf((float)myrandom(RNGen) * (lights.size() - 1));

    return SampleSphere(lights[r]);
}

float Scene::PdfLight(Shape* L)
{
    const float r = static_cast<Sphere*>(L)->Raidus();
    const float LightArea = 4.f * PI * r * r;
    return 1 / (LightArea * lights.size());
}

Color Scene::EvalRadiance(Shape* L)
{
    return L->get_material()->Kd;
}

vec3 Scene::SampleBrdf(Shape* shape, vec3 N, vec3 w)
{
    const float r = myrandom(RNGen);
    const float Kdl = length(shape->get_material()->Kd);
    const float Ksl = length(shape->get_material()->Ks);
    const float Ktl = length(shape->get_material()->Kt);
    const float s = Kdl + Ksl + Ktl;
    const float pd = Kdl / s;
    const float pr = Ksl / s;
    const float pt = Ktl / s;

    const float r1 = myrandom(RNGen);
    const float r2 = myrandom(RNGen);
    if (r < pd)
    {
        return SampleLobe(N, sqrtf(r1), 2.f * PI * r2);
    }
    else if (r < (pd + pr))
    {
        //constexpr float roughness = 0.5f;
        const float alpha =  shape->get_material()->alpha;

        //const float theta = cosf(atanf((alpha * sqrtf(r1)) / sqrtf(1 - r1)));
        const float theta = powf(r1, 1 / (alpha + 1));
        const vec3 m = SampleLobe(N, theta, 2.f * PI * r2);
        return 2.f * abs(dot(w, m)) * m - w;
    }
    else {
        const float alpha = shape->get_material()->alpha;
        const float theta = powf(r1, 1 / (alpha + 1));
        const vec3 m = SampleLobe(N, theta, 2.f * PI * r2);
        float n = 0;
        const float WN = dot(w, N);
        if (WN > 0) {
            n = 1.f / shape->get_material()->ior;
        }
        else {
            n = shape->get_material()->ior / 1.f;
        }

        const float WM = dot(w, m);
        const float re = 1.f - (n * n) * (1.f - (WM * WM));
        if (re < 0) {
            return 2.f * abs(dot(w, m)) * m - w;
        }
        else {
            const float swn = (WN >= 0) ? 1 : -1;
            return (n * WM - swn * sqrtf(re)) * m - n * w;
        }
    }
}

float Scene::GeometryFactory(const Intersection& A, const Intersection& B)
{
    const vec3 D = (A.point - B.point);

    const float DDotD = dot(D, D);
    return std::abs(dot(A.normal, D) * dot(B.normal, D) / (DDotD * DDotD));
}

vec3 Scene::SampleLobe(vec3 A, float c, float theta)
{
    const float s = sqrtf(1.f - (c * c));
    const vec3 K = { s * cosf(theta), s * sinf(theta), c };
    if (abs(A.z - 1) < 10e-3f)
        return K;
    if (abs(A.z + 1) < 10e-3f)
        return { K.x, -K.y, -K.z };

    A = normalize(A);
    const vec3 B = normalize(vec3 { -A.y, A.x, 0 });
    const vec3 C = cross(A, B);
    return K.x * B + K.y * C + K.z * A;
}

Intersection Scene::SampleSphere(Shape* _shape)
{
    const Sphere* sphere = dynamic_cast<Sphere*>(_shape);
    if (!sphere)
        return Intersection();
    const vec3 C = sphere->Center();
    const float R = sphere->Raidus();

    const float r1 = myrandom(RNGen);
    const float r2 = myrandom(RNGen);
    const float z = 2.f * r1 - 1.f;
    const float r = sqrtf(1.f - (z * z));
    const float a = 2.f * PI * r2;

    Intersection intersection;
    intersection.normal = normalize(vec3 { r * cosf(a), r * sinf(a), z });
    intersection.point = C + R * intersection.normal;
    intersection.object = _shape;

    return intersection;
}

void Scene::sphere(const vec3 _center, const float _r, Material* _mat)
{
    Shape* sphere = new Sphere(_center, _r);
    sphere->set_material(_mat);
    shapes.push_back(sphere);
}

void Scene::box(const vec3 _base, const vec3 _diag, Material* _mat)
{
    Shape* shape = new Box(_base, _diag);
    shape->set_material(_mat);
    shapes.push_back(shape);
}

void Scene::cylinder(const vec3 _base, const vec3 _axis, const float _radius, Material* _mat)
{
    Shape* shape = new Cylinder(_base, _axis, _radius);
    shape->set_material(_mat);
    shapes.push_back(shape);
}
