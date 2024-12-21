#include "../microfacet.h"

inline Real square(Real x) { return x * x; }

Spectrum eval_op::operator()(const DisneyGlass &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Spectrum base_clr = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real aniso = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha_x, alpha_y;
    AnisoTransform(roughness, aniso, alpha_x, alpha_y);

    bool inner = dot(vertex.geometric_normal, dir_in) < 0;
    Real eta = inner ? 1 / bsdf.eta : bsdf.eta;

    Real n_dot_in = fabs(dot(frame.n, dir_in));
    Real n_dot_out = fabs(dot(frame.n, dir_out));

    Vector3 half_vector;
    if(reflect)half_vector = normalize(dir_in + dir_out); 
    else {
        half_vector = normalize(-dir_in - eta * dir_out);
        if(length_squared(-dir_in - eta * dir_out) < 1e-3) {
            half_vector = normalize(cross(cross(dir_out, frame.n), dir_out));
        }
    }
    if (dot(half_vector, frame.n) < 0) half_vector = -half_vector;

    Real h_dot_in = dot(half_vector, dir_in); // cos(theta_i)
    Real h_dot_out = dot(half_vector, dir_out); // cos(theta_t)

    
    Real F_g =  fresnel_dielectric(std::abs(h_dot_in), eta);
    Real G_g = smith_masking_gtr2_aniso(to_local(frame, dir_in), alpha_x, alpha_y) *
               smith_masking_gtr2_aniso(to_local(frame, dir_out), alpha_x, alpha_y);
    Real D_g = GTR2Aniso(to_local(frame, half_vector), alpha_x, alpha_y);
    if (reflect)
    {
        return base_clr * F_g * G_g * D_g / (4 * n_dot_in);
    }
    else
    {
        Real denominator = h_dot_in + h_dot_out * eta;
        // return base_clr * (1 - F_g) * G_g * D_g * (h_dot_in * h_dot_out) / n_dot_in * eta * eta / ( denominator * denominator );
        return base_clr * (1 - F_g) * G_g * D_g * (h_dot_in * h_dot_out) / n_dot_in / ( denominator * denominator );
    }
}

Real pdf_sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Spectrum base_clr = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real aniso = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha_x, alpha_y;
    AnisoTransform(roughness, aniso, alpha_x, alpha_y);


    bool inner = dot(vertex.geometric_normal, dir_in) < 0;
    Real eta = inner ? 1 / bsdf.eta : bsdf.eta;

    Vector3 half_vector;
    if(reflect)half_vector = normalize(dir_in + dir_out); 
    else half_vector = normalize(dir_in + eta * dir_out);
    if (dot(half_vector, frame.n) < 0) half_vector = -half_vector;

    Real h_dot_in = dot(half_vector, dir_in);
    Real h_dot_out = dot(half_vector, dir_out);


    Real F_g = fresnel_dielectric(h_dot_in, eta);
    Real D_g = GTR2Aniso(to_local(frame, half_vector), alpha_x, alpha_y);
    if (reflect) return (F_g * D_g) / (4 * fabs(h_dot_in));
    else
    {
        Real numerator = h_dot_in + h_dot_out * eta;
        Real jacobian = eta * eta / (numerator * numerator);
        return (1 - F_g) * D_g * jacobian * fabs(h_dot_in);
    }
}


Vector3 reflect(const Vector3 &i, const Vector3 &h) {
    return normalize(-i + 2 * dot(i, h) * h);
}

Vector3 refract(const Vector3 &i, const Vector3 &h, Real eta) {
    Real cos_theta_i = dot(i, h);
    Real sin2_theta_i = max(Real(0), 1 - cos_theta_i * cos_theta_i);
    Real sin2_theta_t = sin2_theta_i / (eta * eta);

    // Total Internal Reflection
    if (sin2_theta_t >= 1) return Vector3(0);
    auto half = h;
    if(cos_theta_i < 0) half = -h;

    Real cos_theta_t = sqrt(1 - sin2_theta_t);
    return normalize(-i / eta + (fabs(cos_theta_i) / eta - cos_theta_t) * half);
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real aniso = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha_x, alpha_y;
    AnisoTransform(roughness, aniso, alpha_x, alpha_y);

    // Sample a micro normal and transform it to world space -- this is our half-vector.
    bool inner = dot(vertex.geometric_normal, dir_in) < 0;
    Real eta = inner ? 1 / bsdf.eta : bsdf.eta;
    Vector3 local_in = to_local(frame, dir_in);
    Vector3 half_vector = to_world(frame, sample_visible_normals(local_in, alpha_x, alpha_y, rnd_param_uv));
    if (dot(half_vector, frame.n) < 0) half_vector = -half_vector;
    Real F_g = fresnel_dielectric(dot(half_vector, dir_in), eta);
    if (F_g > rnd_param_w)
    {
        Vector3 reflected = reflect(dir_in, half_vector);
        return BSDFSampleRecord{ reflected, Real(0), roughness };
    }
    else
    {
        Vector3 refracted = refract(dir_in, half_vector, eta);
        if(refracted == Vector3(0)) return {};
        return BSDFSampleRecord{ refracted, eta, roughness };
    }
}
TextureSpectrum get_texture_op::operator()(const DisneyGlass &bsdf) const {
    return bsdf.base_color;
}
