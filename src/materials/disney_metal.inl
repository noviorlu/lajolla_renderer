#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyMetal &bsdf) const {
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
    if(length_squared(half_vector) == 0) half_vector = frame.n;

    Real n_dot_h = dot(frame.n, half_vector);
    Real n_dot_in = dot(frame.n, dir_in);
    Real n_dot_out = dot(frame.n, dir_out);
    Real h_dot_out = dot(half_vector, dir_out);
    if (n_dot_out <= 0 || n_dot_h <= 0) {
        return make_zero_spectrum();
    }

    Spectrum base_clr = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real aniso = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha_x, alpha_y;
    AnisoTransform(roughness, aniso, alpha_x, alpha_y);


    Spectrum tint = Spectrum(0.0, 1.0, 0.0);
    Real tint_strength = 1.0;
    
    Spectrum F_m = schlick_fresnel(base_clr, h_dot_out);
    // Spectrum F_m = schlick_generalized_fresnel(base_clr, h_dot_out, 5.0, tint_strength, tint);

    Real D_m = GTR2Aniso(to_local(frame, half_vector), alpha_x, alpha_y);
    Real G_m = smith_masking_gtr2_aniso(to_local(frame, dir_in), alpha_x, alpha_y) *
               smith_masking_gtr2_aniso(to_local(frame, dir_out), alpha_x, alpha_y);

    return F_m * D_m * G_m / (4 * n_dot_in);
}

Real pdf_sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
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
    Vector3 half_vector = normalize(dir_in + dir_out);
    Real n_dot_out = dot(frame.n, dir_out);
    Real n_dot_h = dot(frame.n, half_vector);
    Real h_dot_in = dot(half_vector, dir_in);
    if (n_dot_out <= 0 || n_dot_h <= 0) {
        return 0;
    }

    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real aniso = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha_x, alpha_y;
    AnisoTransform(roughness, aniso, alpha_x, alpha_y);

    Real D = GTR2Aniso(to_local(frame, half_vector), alpha_x, alpha_y);

    return D / (4 * h_dot_in);
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) <= 0) {
        return {};
    }
    
    // Homework 1: implement this!
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = clamp(roughness, Real(0.01), Real(1));
    Real aniso = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);    
    Real alpha_x, alpha_y;
    AnisoTransform(roughness, aniso, alpha_x, alpha_y);

    Vector3 local_micro_normal = sample_visible_normals(to_local(frame, dir_in), alpha_x, alpha_y, rnd_param_uv);
    
    Vector3 half_vector = to_world(frame, local_micro_normal);
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
    return BSDFSampleRecord{
        reflected,
        Real(0) /* eta */, roughness /* roughness */
    };
}

TextureSpectrum get_texture_op::operator()(const DisneyMetal &bsdf) const {
    return bsdf.base_color;
}
