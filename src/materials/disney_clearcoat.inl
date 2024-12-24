#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyClearcoat &bsdf) const {
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
    Vector3 half_vector = normalize(dir_in + dir_out);
    Real n_dot_h = dot(frame.n, half_vector);
    Real n_dot_in = dot(frame.n, dir_in);
    Real n_dot_out = dot(frame.n, dir_out);
    Real h_dot_out = dot(half_vector, dir_out);
    if (n_dot_out <= 0 || n_dot_h <= 0) {
        return make_zero_spectrum();
    }

    Real clearcoatGloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha_g = (1 - clearcoatGloss) * 0.1 + clearcoatGloss * 0.001;

    Real R_0 = 0.04; //(1.5 - 1)^2 / (1.5 + 1)^2

    Real F_c = schlick_fresnel(R_0, h_dot_out);
    Real D_c = GTR1(n_dot_h, alpha_g);
    Real G_c = smith_masking_gtr2_aniso(to_local(frame, dir_in), 0.25, 0.25) *
               smith_masking_gtr2_aniso(to_local(frame, dir_out), 0.25, 0.25);

    return F_c * D_c * G_c / (4 * n_dot_in);
}

Real pdf_sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return 0;
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) <= 0) {
        return 0;
    }
    // Homework 1: implement this!
    Vector3 half_vector = normalize(dir_in + dir_out);
    Real n_dot_out = dot(frame.n, dir_out);
    Real n_dot_h = dot(frame.n, half_vector);
    if (n_dot_out <= 0 || n_dot_h <= 0) {
        return 0;
    }

    Real clearcoatGloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha_g = (1 - clearcoatGloss) * 0.1 + clearcoatGloss * 0.001;

    Real D_c = GTR1(n_dot_h, alpha_g);

    return D_c * n_dot_h / (4 * n_dot_out);
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) <= 0) {
        return {};
    }
    // Homework 1: implement this!
    // rnd_param_uv
    Real clearcoatGloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha = (1 - clearcoatGloss) * 0.1 + clearcoatGloss * 0.001;
    Real a2 = alpha * alpha;

    Real h_elevation = acos(sqrt((1-pow(a2, 1-rnd_param_uv[0])) / (1-a2)));
    Real h_azimuth = 2 * c_PI * rnd_param_uv[1];

    Vector3 h_local;
    h_local.x = sin(h_elevation) * cos(h_azimuth);
    h_local.y = sin(h_elevation) * sin(h_azimuth);
    h_local.z = cos(h_elevation);

    Vector3 half_vector = to_world(frame, h_local);
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
    return BSDFSampleRecord{
        reflected,
        Real(0) /* eta */, clearcoatGloss /* roughness */
    };
}

TextureSpectrum get_texture_op::operator()(const DisneyClearcoat &bsdf) const {
    return make_constant_spectrum_texture(make_zero_spectrum());
}
