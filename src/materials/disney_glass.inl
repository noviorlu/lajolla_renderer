#include "../microfacet.h"

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


    bool inner = dot(frame.n, dir_in) < 0;
    Real eta = inner ? bsdf.eta : Real(1) / bsdf.eta; // internal IOR t / externalIOR i (dir_in is external to material)
    
    Real n_dot_in = fabs(dot(frame.n, dir_in));
    Real n_dot_out = fabs(dot(frame.n, dir_out));
    
    Vector3 half_vector;
    if(reflect)half_vector = normalize(dir_in + dir_out); 
    else half_vector = normalize(-dir_in - eta * dir_out);
    half_vector = dot(frame.n, half_vector) > 0 ? half_vector : -half_vector; // sync half_vec in the same side with frame normal

    Real h_dot_in = fabs(dot(half_vector, dir_in)); // cos(theta_i)
    Real h_dot_out = fabs(dot(half_vector, dir_out)); // cos(theta_t)


    bool criticCheck = (1 - h_dot_out * h_dot_out) * eta * eta < 1;
    Real F_g = (criticCheck) ? fresnel_dielectric(h_dot_in, h_dot_out, eta) : Real(1);
    // doesnt matter if reflection or refraction since always squared
    Real G_g = smith_masking_gtr2_aniso(to_local(frame, dir_in), alpha_x, alpha_y) *
               smith_masking_gtr2_aniso(to_local(frame, dir_out), alpha_x, alpha_y);
    Real D_g = GTR2Aniso(to_local(frame, half_vector), alpha_x, alpha_y);

    Real GGX;
    if(reflect) GGX = F_g * G_g * D_g / (4 * n_dot_in);
    else {
        Real denominator = h_dot_in + h_dot_out * eta;
        GGX = (1 - F_g) * G_g * D_g * (h_dot_in * h_dot_out) / n_dot_in * eta * eta / ( denominator * denominator );
    }

    return base_clr * GGX;
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


    bool inner = dot(frame.n, dir_in) < 0;
    Real eta = inner ? bsdf.eta : Real(1) / bsdf.eta; // internal IOR t / externalIOR i (dir_in is external to material)

    Real n_dot_in = fabs(dot(frame.n, dir_in));
    Real n_dot_out = fabs(dot(frame.n, dir_out));
    
    Vector3 half_vector;
    if(reflect)half_vector = normalize(dir_in + dir_out); 
    else half_vector = normalize(-dir_in - eta * dir_out);
    half_vector = dot(frame.n, half_vector) > 0 ? half_vector : -half_vector; // sync half_vec in the same side with frame normal

    Real h_dot_in = fabs(dot(half_vector, dir_in)); // cos(theta_i)
    Real h_dot_out = fabs(dot(half_vector, dir_out)); // cos(theta_t)


    bool criticCheck = (1 - h_dot_out * h_dot_out) * eta * eta < 1;
    Real F_g = (criticCheck) ? fresnel_dielectric(h_dot_in, h_dot_out, eta) : Real(1);
    Real D_g = GTR2Aniso(to_local(frame, half_vector), alpha_x, alpha_y);

    Real pdf_reflect = D_g / (4 * h_dot_in);
    Real numerator = h_dot_in + h_dot_out * eta;
    Real pdf_refract = D_g * h_dot_out / (numerator * numerator);

    return F_g * pdf_reflect + (1 - F_g) * pdf_refract;
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

    Real cos_theta_t = sqrt(1 - sin2_theta_t);
    return normalize(eta * -i + (eta * cos_theta_i - cos_theta_t) * h);
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
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


    bool inner = dot(frame.n, dir_in) < 0;
    Real eta = inner ? bsdf.eta : Real(1) / bsdf.eta; // internal IOR t / externalIOR i (dir_in is external to material)


    Vector3 local_in = inner ? to_local(frame, -dir_in) : to_local(frame, dir_in);
    Vector3 half_vector = to_world(frame, sample_visible_normals(local_in, alpha_x, alpha_y, rnd_param_uv));
    half_vector = dot(dir_in, half_vector) > 0 ? half_vector : -half_vector;

    Real h_dot_in = fabs(dot(half_vector, dir_in)); // cos(theta_i)
    

    Vector3 dir_out = reflect(dir_in, half_vector);
    Real F_g = fresnel_dielectric(h_dot_in, eta);
    if(F_g < rnd_param_w){
        Vector3 refracted = refract(dir_in, half_vector, eta);
        if(refracted.x != 0 || refracted.y != 0 || refracted.z != 0)
            dir_out = refracted;
    }

    return BSDFSampleRecord{
        dir_out,
        eta, 
        roughness
    };
}

TextureSpectrum get_texture_op::operator()(const DisneyGlass &bsdf) const {
    return bsdf.base_color;
}
