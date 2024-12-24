#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneySheen &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return make_zero_spectrum();
    }
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) <= 0) {
        return make_zero_spectrum();
    }

    // Homework 1: implement this!
    Spectrum base_clr = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen_tint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real l = luminance(base_clr);
    Spectrum tint = Spectrum(1);
    if (l > 0) tint = base_clr / l;

    Spectrum sheen = (1 - sheen_tint) + sheen_tint * tint;
    Vector3 half_vector = normalize(dir_in + dir_out);
    
    return sheen * pow(1 - abs(dot(half_vector, dir_out)), 5) * max(Real(0), dot(frame.n, dir_out));
}

Real pdf_sample_bsdf_op::operator()(const DisneySheen &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return 0;
    }
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) <= 0) {
        return 0;
    }

    // Homework 1: implement this!
    return fmax(dot(frame.n, dir_out), Real(0)) / c_PI;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneySheen &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) <= 0) {
        return {};
    }

    // Homework 1: implement this!
    return BSDFSampleRecord{
        to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
        Real(0) /* eta */, Real(1) /* roughness */ };
}

TextureSpectrum get_texture_op::operator()(const DisneySheen &bsdf) const {
    return bsdf.base_color;
}
