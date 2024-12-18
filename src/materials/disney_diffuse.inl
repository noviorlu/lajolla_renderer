#include "lajolla.h"

Spectrum eval_op::operator()(const DisneyDiffuse &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return make_zero_spectrum();
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    // Homework 1: implement this!
    Vector3 h = normalize(dir_in + dir_out);
    Real hdout = max(0.0, dot(h, dir_out));
    Real ndin = max(0.0, dot(frame.n, dir_in));
    Real ndout = max(0.0, dot(frame.n, dir_out));
    
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real ss = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);

    auto FresnelSchlick = [](Real F0, Real cosTheta) {
        return 1.0, F0 + (1.0 - F0) * pow(1.0 - cosTheta, 5);
    };

    // 2012 Disney BRDF
    // Real FSS90 = roughness * hdout * hdout;
    // Real FD90 = 0.5 + 2 * FSS90;
    // Spectrum baseDiffuse = base_color / c_PI * FresnelSchlick(FD90, ndin) * FresnelSchlick(FD90, ndout) * ndout;
    // Spectrum subsurface = 1.25 * base_color / c_PI 
    //         * (FresnelSchlick(FSS90, ndin) * FresnelSchlick(FSS90, ndout) * (1/(ndin+ndout) - 0.5) + 0.5 ) * ndout;
    // return (1 - ss) * baseDiffuse + ss * subsurface;
    
    // 2015 Disney BSDF
    Real FL = pow((1 - ndin), Real(5));
    Real FV = pow((1 - ndout), Real(5));
    Real FSS90 = roughness * hdout * hdout;
    Real Rr = 2 * FSS90;

    Spectrum lambert = base_color / c_PI;
    Spectrum subsurface = 1.25 * base_color / c_PI 
            * (FresnelSchlick(FSS90, ndin) * FresnelSchlick(FSS90, ndout) * (1/(ndin+ndout) - 0.5) + 0.5 ) * ndout;
    Spectrum retro_reflection = lambert * Rr * (FL + FV + FL * FV * (Rr - 1));
    return ((1 - ss) * lambert + ss * subsurface) * (1 - 0.5 * FL) * (1 - 0.5 * FV) + retro_reflection;
}

Real pdf_sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return 0;
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    
    // Homework 1: implement this!
    return fmax(dot(frame.n, dir_out), Real(0)) / c_PI;
}

std::optional<BSDFSampleRecord> sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    
    // Homework 1: implement this!
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    return BSDFSampleRecord{
        to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
        Real(0) /* eta */, Real(roughness) /* roughness */};
}

TextureSpectrum get_texture_op::operator()(const DisneyDiffuse &bsdf) const {
    return bsdf.base_color;
}
