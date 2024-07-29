
#include <vector>
#include "geom.h"
#include "raytrace.h"
#include "acceleration.h"

#include <bvh/sweep_sah_builder.hpp>
#include <bvh/single_ray_traverser.hpp>
#include <bvh/primitive_intersectors.hpp>

/////////////////////////////
// Vector and ray conversions
Ray RayFromBvh(const bvh::Ray<float>& r) {
    return Ray(vec3FromBvh(r.origin), vec3FromBvh(r.direction));
}
bvh::Ray<float> RayToBvh(const Ray& r) {
    return bvh::Ray<float>(vec3ToBvh(r.o), vec3ToBvh(r.d));
}


/////////////////////////////
// SimpleBox
bvh::Vector3<float> vec3ToBvh(const vec3& v) {
    return bvh::Vector3<float>(v[0], v[1], v[2]);
}

vec3 vec3FromBvh(const bvh::Vector3<float>& v) {
    return vec3(v[0], v[1], v[2]);
}

SimpleBox::SimpleBox() : bvh::BoundingBox<float>() {}
SimpleBox::SimpleBox(const vec3 v) : bvh::BoundingBox<float>(vec3ToBvh(v)) {}

SimpleBox& SimpleBox::extend(const vec3 v) {
    bvh::BoundingBox<float>::extend(vec3ToBvh(v));
    return *this;
}


/////////////////////////////
// BvhShape

SimpleBox BvhShape::bounding_box() const {
    //  Return the shape's bounding box.
    return shape->bounding_box(); // FIX THIS
}

bvh::Vector3<float> BvhShape::center() const {
    return bounding_box().center();
}

std::optional<Intersection> BvhShape::intersect(const bvh::Ray<float>& bvhray) const {
    // Intersect RayFromBvh(bvhray) with shape;  store result in an Intersection
    // If no intersection,
    //    return std::nullopt;
    // If intersection's t value < bvhray.tmin  or > bvhray.tmax
    //    return std::nullopt;
    // else return
    //    return the Intersection

    const Ray ray = RayFromBvh(bvhray);
    Intersection intersection;
    if (shape->intersect(ray, intersection)) {
        if (intersection.t < bvhray.tmin || intersection.t > bvhray.tmax)
            return std::nullopt;
        else
            return intersection;
    }

    return std::nullopt;
}

AccelerationBvh::AccelerationBvh(std::vector<Shape*>& objs) {
    // Wrap all Shape*'s with a bvh specific instance
    for (Shape* shape : objs) {
        shapeVector.emplace_back(shape);
    }

    // Magic found in the bvh examples:
    auto [bboxes, centers] = bvh::compute_bounding_boxes_and_centers(shapeVector.data(),
        shapeVector.size());
    auto global_bbox = bvh::compute_bounding_boxes_union(bboxes.get(), shapeVector.size());

    bvh::SweepSahBuilder<bvh::Bvh<float>> builder(bvh);
    builder.build(global_bbox, bboxes.get(), centers.get(), shapeVector.size());
}

Intersection AccelerationBvh::intersect(const Ray& ray) {
    bvh::Ray<float> bvhRay = RayToBvh(ray);

    // Magic found in the bvh examples:
    bvh::ClosestPrimitiveIntersector<bvh::Bvh<float>, BvhShape> intersector(bvh, shapeVector.data());
    bvh::SingleRayTraverser<bvh::Bvh<float>> traverser(bvh);

    auto hit = traverser.traverse(bvhRay, intersector);
    if (hit) {
        return hit->intersection;
    }
    else
        return  Intersection();  // Return an IntersectionRecord which indicates NO-INTERSECTION


}

vec3 Ray::eval(float _t) const {
    return o + (_t * d);
}


vec3 Shape::FresnelFactor(const float d, const float WoN) {
    if (0) {//(material->Kt != vec3(0)) {

        const float c = abs(d);
        float no;
        float ni;

        if (WoN < 0) {
            ni = 1.f;
            no = material->ior;
        }
        else {
            no = 1.f;
            ni = material->ior;
        }
        const float n = ni / no;
        const float g = sqrtf(((no * no) / (ni * ni)) - 1 + (c * c));
        const float gmc = g - c;
        const float gpc = g + c;
        const float one = (gmc * gmc) / (gpc * gpc);
        const float num = c * gpc - 1;
        const float den = c * gmc + 1;
        const float two = 1.f + ((num * num) / (den * den));
        return vec3(0.5f * one * two);
    }
    else {
        const vec3 Ks = material->Ks;
        return Ks + (1.f - Ks) * powf(1.f - abs(d), 5.f);
    }
}

float Shape::Distribution(const vec3& m, const vec3& N) {
    // const float d = dot(m, N) > 0 ? 1 : 0;
    // //constexpr float alpha = 0.5f * 0.5f;
    // const float alpha = length(material->Ks);
    // const float alpha2 = alpha * alpha;
    //
    // const float mDotN = dot(m, N);
    // const float theta = sqrtf(1.0f - powf(mDotN, 2)) / mDotN;
    // const float denom = PI * powf(mDotN, 4) * powf(alpha2 + (theta * theta), 2);
    //
    // return d * alpha2 / denom;
    const float mDotN = dot(m, N);
    const float d = mDotN > 0 ? 1 : 0;
    return d * ((material->alpha + 2) / (2 * PI)) * powf(mDotN, material->alpha);
}

float Shape::G1GGX(const vec3& v, const vec3& m, const vec3& N) {
    const float vDotN = dot(v, N);
    if (vDotN > 1)
        return 1;
    const float d = (dot(v, m) / vDotN) > 0 ? 1 : 0;
    const float theta = sqrtf(1.0f - (vDotN * vDotN)) / vDotN;// std::max(vDotN, 0.001f);
    if (abs(theta) < 10e-6)
        return 1;
    //const float alpha = length(material->Ks);
    //const float denom = 1.f + sqrtf(1 + (alpha * alpha) * (theta * theta));
    //return d * 2.f / denom;

    const float a = sqrtf(material->alpha / 2 + 1) / theta;
    float stuff = 1;
    if (a < 1.6f) {
        const float num = 3.535 * a + 2.181 * a * a;
        const float den = 1.f + 2.276 * a + 2.577 * a * a;
        stuff = num / den;
    }

    return d * stuff;
}

float Shape::PdfBrdf(const vec3& wo, const vec3& N, const vec3& wi) {
    const float Pd = abs(dot(wi, N)) / PI;
    vec3 m = normalize(wo + wi);
    const float Pr = Distribution(m, N) * abs(dot(m, N)) * (1.f / (4 * abs(dot(wi, m))));
    const float s = length(material->Kd) + length(material->Ks) + length(material->Kt);
    const float pd = length(material->Kd) / s;
    const float pr = length(material->Ks) / s;
    const float pt = length(material->Kt) / s;

    const float WoN = dot(wo, N);
    float no;
    float ni;

    if (WoN > 0) {
        ni = 1.f;
        no = material->ior;
    }
    else {
        no = 1.f;
        ni = material->ior;
    }
    const float n = ni / no;

    m = -(no * wi + ni * wo);
    m = normalize(m);
    const float wom = dot(wo, m);
    const float r = 1.f - (n * n) * (1.f - (wom * wom));
    float Pt = 0;
    if (r < 0) {
        Pt = Pr;
    }
    else {
        const float den = (no * dot(wi, m) + ni * dot(wo, m));
        Pt = Distribution(m, N) * abs(dot(m, N)) * (no * no) * abs(dot(wi, m)) / (den * den);
    }

    return pd * Pd + pr * Pr + pt * Pt;
}

vec3 Shape::EvalScattering(const vec3& wi, const vec3& N, const vec3& wo, const float t) {
    const float s = length(material->Kd) + length(material->Ks) + length(material->Kt);
    const float pd = length(material->Kd) / s;
    const float pr = length(material->Ks) / s;
    const float pt = length(material->Kt) / s;
    const vec3 Ed = (pd > 0) ? material->Kd / PI : vec3(0, 0, 0);
    const vec3 m = normalize(wo + wi);
    const vec3 Er = (pr > 0) ? (
        Distribution(m, N)
        * GeometryGGX(wi, wo, m, N)
        * FresnelFactor(dot(wi, m), dot(wo, N)))
        / (4.f * abs(dot(wi, N)) * abs(dot(wo, N))) : vec3(0, 0, 0);
    vec3 Et = { 0,0,0 };
    if (pt > 0) {
        const float WoN = dot(wo, N);
        float no;
        float ni;
        vec3 A;

        if (WoN >= 0) {
            ni = 1.f;
            no = material->ior;
            A = vec3(1);
        }
        else {
            no = 1.f;
            ni = material->ior;
            A = glm::pow(vec3(2.71828f), t * glm::log(material->Kt));
            //float temp = powf(2.71828f,  t * std::log(material->Kt.x));
            //float temp2 = powf(2.71828f, t * std::log(material->Kt.y));
            //float temp3 = powf(2.71828f, t * std::log(material->Kt.z));
            //A = vec3{ temp, temp2, temp3 };
        }
        const float n = ni / no;
        vec3 m2 = -(no * wi + ni * wo);
        m2 = normalize(m2);
        const float wom = dot(wo, m2);
        const float r = 1 - (n * n) * (1.f - (wom * wom));
        if (r < 0) {
            Et = Distribution(m, N)
                * GeometryGGX(wi, wo, m, N)
                * FresnelFactor(dot(wi, m), WoN)
                / (4.f * abs(dot(wi, N)) * abs(dot(wo, N)));
        }
        else {
            Et = Distribution(m2, N) * GeometryGGX(wi, wo, m2, N) * (vec3(1, 1, 1) - FresnelFactor(dot(wi, m2), WoN)) /
                (abs(dot(wi, N)) * abs(dot(wo, N)));
            const float den = (no * dot(wi, m2) + ni * dot(wo, m2));
            Et *= abs(dot(wi, m2)) * abs(dot(wo, m2)) * (no * no) / (den * den);
        }
        Et = { Et.x * A.x, Et.y * A.y, Et.z * A.z };
    }
    return abs(dot(N, wi))* (Ed + Er + Et);
}

