#ifndef _ACC_H
#define _ACC_H

// This uses a library called bvh for a bounding volume hierarchy acceleration
// structure. See https://github.com/madmann91/bvh

// The bvh library uses its own version of vectors, bounding boxes,
// rays, and list of shapes, so there will be some conversions to/from
// the raytracer's versions of these structures.  Alternatively, the
// raytracer could be written to directly use bvh's version of these
// structures.  This seems feasible, but has not been tested.

#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif

#include <bvh/bvh.hpp>
#include <bvh/ray.hpp>
#include <bvh/vector.hpp>
#include <optional>
#include "geom.h"

class Shape;
// FIX THIS:  Dummy ray to allow for compilation
class Ray {
public:
    vec3 o, d;
    Ray(const vec3 _o, const vec3 _d)
        : o(_o)
        , d(_d)
    {
    }
    vec3 eval(float _t) const;
};

// Vectors:
// Expectation: The raytracer uses glm::vec3 throughout
// These convert between glm::vec3 and bvh::Vector3<float>
bvh::Vector3<float> vec3ToBvh(const vec3& v);
vec3 vec3FromBvh(const bvh::Vector3<float>& v);

// Bounding boxes
// Expectation: The raytracer uses SimpleBox throughout.
// This class derives from the bvh bounding box, overloading two methods to take
// glm::vec3. Use: Construct a BB with zero or one points (vec3), extend it with
// further points.
//      Access its bounds as two bvh vectors box.min and box.max.
//bvh::Vector3<float> vec3ToBvh(const vec3& v);
//vec3 vec3FromBvh(const bvh::Vector3<float>& v);

class SimpleBox : public bvh::BoundingBox<float> {
public:
    SimpleBox();
    SimpleBox(const vec3 v);
    SimpleBox& extend(const vec3 v);
};

// FIX THIS:
// Rays:  Rays have an origin and a direction.  bvh::Ray also has a tmin and a
// tmax. Supply your own ray class and convert.
bvh::Ray<float> RayToBvh(const Ray& r);
Ray RayFromBvh(const bvh::Ray<float>& r);

// FIX THIS:  This dummy Intersection record is defficient -- just barely enough
// to compile.
class Intersection {
public:
    float t = INFINITY;
    Shape* object = nullptr;
    vec3 point{};
    vec3 normal{};

    void operator=(const Intersection& _rhs)
    {
        t = _rhs.t;
        object = _rhs.object;
        point = _rhs.point;
        normal = _rhs.normal;
    }

    float distance() const noexcept
    {
        return t;
    } // A function the BVH traversal needs to be supplied.
};

class Interval {
public:
    float t0, t1;
    vec3 N0{};
    vec3 N1{};

    Interval() noexcept
        : t0(0)
        , t1(INFINITY)
    {
    }

    Interval(float _t0, float _t1, vec3 _N0, vec3 _N1)
        : t0(_t0)
        , t1(_t1)
        , N0(_N0)
        , N1(_N1)
    {
        if (t0 > t1) {
            std::swap(t0, t1);
            std::swap(N0, N1);
        }
    }

    void empty() noexcept
    {
        t0 = 0;
        t1 = -1;
    }

    void intersect(const Interval& _interval)
    {
        if (t0 < _interval.t0) {
            t0 = _interval.t0;
            N0 = _interval.N0;
        }

        if (t1 > _interval.t1) {
            t1 = _interval.t1;
            N1 = _interval.N1;
        }
    }

    void intersect(const vec3& N, const float& d0, const float& d1, const Ray& _ray)
    {
        const float nDotD = glm::dot(N, _ray.d);
        const float nDotQ = glm::dot(N, _ray.o);
        if (nDotD != 0) {
            float T0 = -(d0 + nDotQ) / nDotD;
            float T1 = -(d1 + nDotQ) / nDotD;
            vec3 TN0 = -N;
            vec3 TN1 = N;
            if (T0 > T1) {
                std::swap(T0, T1);
                std::swap(TN0, TN1);
            }
            if (t0 < T0) {
                t0 = T0;
                N0 = TN0;
            }
            if (t1 > T1) {
                t1 = T1;
                N1 = TN1;
            }
        } else {
            float s0 = nDotQ + d0;
            float s1 = nDotQ + d1;
            if ((s0 < 0 && s1 > 0) || (s0 > 0 && s1 < 0)) {
                t0 = std::max(0.0f, t0);
                t1 = std::min(t1, INFINITY);
            } else {
                t0 = std::max(1.0f, t0);
                t1 = std::min(0.0f, t1);
            }
        }
    }
};

class Material;
// FIX THIS; A dummy class -- just barely enough to compile.
class Shape {
protected:
    SimpleBox boundingBox;
    Material* material;

public:
    virtual ~Shape() {}
    _Success_(return == true) virtual bool intersect(_In_ const Ray& _ray,
        _Out_ Intersection& _intersection)
        = 0;
    SimpleBox& bounding_box()
    {
        return boundingBox;
    }

    Material* get_material()
    {
        return material;
    }

    void set_material(Material* _material)
    {
        material = _material;
    }

    vec3 FresnelFactor(const float d, const float WoN);
    
    float Distribution(const vec3& m, const vec3& N);
    float PdfBrdf(const vec3& wo, const vec3& N, const vec3& wi);

    float GeometryGGX(const vec3& wi, const vec3& w0, const vec3& m, const vec3& N)
    {
        return G1GGX(wi, m, N) * G1GGX(w0, m, N);
    }
    float G1GGX(const vec3& v, const vec3& m, const vec3& N);
    vec3 EvalScattering(const vec3& wi, const vec3& N, const vec3& wo, const float t);
};

class Sphere : public Shape {
    vec3 center;
    float radius;

public:
    const vec3 Center() const
    {
        return center;
    }

    const float Raidus() const
    {
        return radius;
    }
    Sphere(const vec3& _center, const float& _radius)
        : center(_center)
        , radius(_radius)
    {
        boundingBox = SimpleBox(center - radius).extend(center + radius);
    }

    _Success_(return == true) bool intersect(_In_ const Ray& _ray,
        _Out_ Intersection& _intersection) final
    {
        const vec3 Q = _ray.o - center;
        const float QD = dot(Q, normalize(_ray.d));
        const float v = sqrtf((QD * QD) - dot(Q, Q) + (radius * radius));
        const float t0 = -QD - v;
        const float t1 = -QD + v;

        if (t0 > 0) {
            _intersection.point = _ray.o + (t0 * _ray.d);
            _intersection.normal = normalize(_intersection.point - center);
            _intersection.t = t0;
            _intersection.object = this;
            return true;
        } else if (t1 > 0) {
            _intersection.point = _ray.o + (t1 * _ray.d);
            _intersection.normal = normalize(_intersection.point - center);
            _intersection.t = t1;
            _intersection.object = this;
            return true;
        }

        return false;
    }
};

class Box : public Shape {
    vec3 center;
    vec3 diagonal;

public:
    Box(const vec3& _center, const vec3& _diagonal)
        : center(_center)
        , diagonal(_diagonal)
    {
        boundingBox = SimpleBox(center).extend(center + diagonal);
    }
    _Success_(return == true) bool intersect(_In_ const Ray& _ray, _Out_ Intersection& _intersection) final
    {
        Interval interval;
        {
            static const vec3 N{ 1, 0, 0 };
            const float d0 = -center.x;
            const float d1 = -center.x - diagonal.x;
            interval.intersect(N, d0, d1, _ray);
        }
        {
            static const vec3 N2 = { 0, 1, 0 };
            const float d0 = -center.y;
            const float d1 = -center.y - diagonal.y;
            interval.intersect(N2, d0, d1, _ray);
        }
        {
            static const vec3 N3 = { 0, 0, 1 };
            const float d0 = -center.z;
            const float d1 = -center.z - diagonal.z;
            interval.intersect(N3, d0, d1, _ray);
        }
        if (interval.t0 > interval.t1)
            return false;
        else if (interval.t0 > 10e-4f) {
            _intersection.t = interval.t0;
            _intersection.object = this;
            _intersection.normal = normalize(interval.N0);
            _intersection.point = _ray.o + _intersection.t * _ray.d;
            return true;
        } else if (interval.t1 > 10e-4f) {
            _intersection.t = interval.t1;
            _intersection.object = this;
            _intersection.normal = normalize(interval.N1);
            _intersection.point = _ray.o + _intersection.t * _ray.d;
            return true;
        }
        return false;
    }
};

class Cylinder : public Shape {
    vec3 base;
    vec3 axis;
    float radius;

public:
    Cylinder(const vec3& _base, const vec3& _axis, const float& _radius)
        : base(_base)
        , axis(_axis)
        , radius(_radius)
    {
        boundingBox = SimpleBox(base + radius)
                          .extend(base - radius)
                          .extend(axis + base + radius)
                          .extend(axis + base - radius);
    }

    _Success_(return == true) bool intersect(_In_ const Ray& _ray, _Out_ Intersection& _intersection) final
    {
        // vec3 N = { 0, 0, 1 };
        // return false;
        const float d0 = 0;
        Interval interval;

        const vec3 A = normalize(axis);

        const float d1 = -length(axis);
        const bool zaxis = (A == vec3(0, 0, 1));
        const vec3 B = normalize(cross(
            dot(A, { 0, 0, 1 }) == 1 ? Xaxis() : Zaxis(), A));
        const vec3 C = cross(A, B);
        const glm::mat3 R = zaxis ? glm::mat3(1.0f) : glm::transpose(glm::mat3(B, C, A));

        const vec3 Q = R * (_ray.o - base);
        const vec3 D = R * _ray.d;

        const Ray tempRay(Q, D);
        interval.intersect({ 0, 0, 1 }, d0, d1, tempRay);
        if (interval.t0 > interval.t1)
            return false;

        const float a = (D.x * D.x + D.y * D.y);
        const float b = 2.0f * (D.x * Q.x + D.y * Q.y);
        const float c = (Q.x * Q.x + Q.y * Q.y - radius * radius);
        const float discriminant = (b * b) - (4 * a * c);
        if (discriminant < 0)
            return false;

        const float disc = sqrt(discriminant);

        Interval interval2;

        const float b0 = (-b - disc) / (2 * a);
        const float b1 = (-b + disc) / (2 * a);

        if (b0 > b1)
            return false;
       
        interval2 = Interval(b0, b1,
            normalize(vec3(Q.x + b0 * D.x, Q.y + b0 * D.y, 0)),
            normalize(vec3(Q.x + b1 * D.x, Q.y + b1 * D.y, 0)));

        interval2.intersect(interval);
        const float t0 = interval2.t0;
        const float t1 = interval2.t1;
        if (t0 > t1)
            return false;

        if (t0 >= 10e-4f) {
            // Interval imin = interval.t0 > interval2.t0 ? interval : interval2;
            _intersection.t = t0;
            _intersection.normal = normalize(glm::transpose(R)* interval2.N0);
        } else if (t1 >= 10e-4f) {
            _intersection.t = t1;
            _intersection.normal = normalize(glm::transpose(R) * interval2.N1);
        } else {
            return false;
        }

        _intersection.point = _ray.eval(_intersection.t);
        //_intersection.intersectionNormal = glm::transpose(R) * _intersection.intersectionNormal;
        _intersection.object = this;
        return true;
    }
};

class Triangle : public Shape {
    vec3 v0, v1, v2;
    vec3 n0, n1, n2;

public:
    Triangle(const vec3& _v0, const vec3& _v1, const vec3& _v2)
        : v0(_v0)
        , v1(_v1)
        , v2(_v2)
    {
        boundingBox = SimpleBox(v0).extend(v1).extend(v2);
    }

    Triangle(const vec3& _v0, const vec3& _v1, const vec3& _v2,
        const vec3& _n0, const vec3& _n1, const vec3& _n2)
        : Triangle(_v0, _v1, _v2)
    {
        n0 = _n0;
        n1 = _n1;
        n2 = _n2;
    }
    _Success_(return == true) bool intersect(_In_ const Ray& _ray, _Out_ Intersection& _intersection) final
    {
        const vec3 E1 = v1 - v0;
        const vec3 E2 = v2 - v0;
        const vec3 p = glm::cross(_ray.d, E2);
        const float d = glm::dot(p, E1);
        if (d == 0)
            return false;

        const vec3 S = _ray.o - v0;
        const float u = glm::dot(p, S) / d;
        if (u < 0 || u > 1)
            return false;

        const vec3 q = glm::cross(S, E1);
        const float v = glm::dot(_ray.d, q) / d;

        if (v < 0 || (u + v) > 1)
            return false;

        const float t = glm::dot(E2, q) / d;
        if (t < 10e-4f)
            return false;

        _intersection.point = _ray.o + t * _ray.d;
        _intersection.normal = normalize((1 - u - v) * n0 + u * n1 + v * n2);
        _intersection.object = this;
        _intersection.t = t;
        return true;
    }
};

// Wrapper of a single Shape*:

// The ray tracer stores the scene objects as a list of Shape*.  BVH
// expects a list of NON-POINTERS with various methods and type
// declarations.

class BvhShape {
    Shape* shape;

public:
    using ScalarType = float; // Float or double for rays and vectors.
    using IntersectionType = Intersection; // Specify the intersection record type.

    BvhShape(Shape* s)
        : shape(s) {}; // Constructor given the Shape to wrap

    SimpleBox bounding_box() const; // Returns the bounding box of the shape
    bvh::Vector3<float> center() const; // Returns the bounding_box().center()

    // The intersection routine.
    // Given a bvh::Ray, intersect it with the shape and return either of:
    //     an Intersection
    //         ( it exists, and is between the ray's tmin and tmax.
    //     std::nullopt ((otherwise)
    std::optional<Intersection> intersect(const bvh::Ray<float>& bvhray) const;
};

// Encapsulates the BVH structure, the list of shapes it's built from,
// and method to intersect a ray with the full scene and return the
// front most intersection point.
class AccelerationBvh {
    bvh::Bvh<float> bvh;
    std::vector<BvhShape> shapeVector;

public:
    AccelerationBvh(std::vector<Shape*>& objs);
    Intersection intersect(const Ray& ray);
};

#endif
