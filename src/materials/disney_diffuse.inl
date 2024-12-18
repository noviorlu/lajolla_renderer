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
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real ss = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);

    Vector3 h = dir_in + dir_out;
    if(length_squared(h) == 0) h = frame.n;
    h = normalize(h);
    Real hdout = max(0.0, dot(h, dir_out));
    Real ndin = max(0.0, dot(frame.n, dir_in));
    Real ndout = max(0.0, dot(frame.n, dir_out));

    auto FresnelSchlick = [](Real F0, Real cosTheta) {
        return 1.0 + (F0 - 1.0) * pow(1.0 - cosTheta, 5);
    };

    Real FSS90 = roughness * hdout * hdout;
    Real ssCoeff = 1.25 * (FresnelSchlick(FSS90, ndin) * FresnelSchlick(FSS90, ndout) * (1.0 / (ndin + ndout) - 0.5) + 0.5);

    /* 2012 Disney */
    // if(ndout == 0) return make_zero_spectrum();
    // Real FD90 = 0.5 + 2 * FSS90;
    // Real baseCoeff = FresnelSchlick(FD90, ndin) * FresnelSchlick(FD90, ndout);
    // return base_color / c_PI * ( (1.0 - ss) * baseCoeff + ss * ssCoeff ) * ndout;
    
    /* 2015 Disney */
    if(ndout == 0) return make_zero_spectrum();
    Real FL = pow((1 - ndin), Real(5));
    Real FV = pow((1 - ndout), Real(5));
    Real Rr = 2 * FSS90;
    Real retroCoeff = Rr * (FL + FV + FL * FV * (Rr - 1));
    Real lambertCoeff = (1.0 - ss) + ss * ssCoeff;
    return base_color / c_PI * (lambertCoeff * (1 - 0.5 * FL) * (1 - 0.5 * FV) + retroCoeff) * ndout;

    /* Autodesk Diffuse */
    // Real sigma = roughness * roughness;
    // Real A = 1.0 - (sigma / (2.0 * (sigma + 0.33)));
    // Real B = 0.45 * sigma / (sigma + 0.09);
    // Real thetaI = std::acos(clamp(ndin, -1.0, 1.0));
    // Real thetaO = std::acos(clamp(ndout, -1.0, 1.0));
    // Real alpha = std::max(thetaI, thetaO);
    // Real beta = std::min(thetaI, thetaO);
    // Vector3f wo_proj = normalize(dir_out - frame.n * ndout);
    // Vector3f wi_proj = normalize(dir_in - frame.n * ndin);
    // Real cosPhiDiff = dot(wo_proj, wi_proj);
    // cosPhiDiff = max(0.0, cosPhiDiff);
    // Real sinAlpha = std::sin(alpha);
    // Real tanBeta = std::tan(beta);

    // // Oren-Nayar BRDF
    // Real orenNayarCoeff = A + B * cosPhiDiff * sinAlpha * tanBeta;
    // return base_color / c_PI * ( (1.0 - ss) * orenNayarCoeff + ss * ssCoeff ) * ndout;
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
