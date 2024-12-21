#include "../microfacet.h"
#include <vector>

void weightComputation(
    Real specular_transmission,
    Real metallic,
    Real clearcoat,
    Real sheen,
    std::vector<Real>& weights
){
    weights.resize(5);
    weights[0] = (1 - metallic) * (1 - specular_transmission);
    weights[1] = 0.25 * clearcoat;
    weights[2] = (1 - specular_transmission * (1 - metallic));
    weights[3] = (1 - metallic) * specular_transmission;
    weights[4] = sheen * (1 - metallic);
}

Spectrum eval_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }

    // Homework 1: implement this!
    if(dot(vertex.geometric_normal, dir_in) <= 0){
        // only glass
        DisneyGlass glass = {bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta};
        return (*this)(glass);
    }
    else{
        Spectrum result;

        
        Spectrum base_clr = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real cc = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real sh = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real eta = bsdf.eta;
        
        Spectrum C0;
        {
            Real l = luminance(base_clr);
            Spectrum Ctint = Spectrum(1);
            if (l > 0) Ctint = base_clr / l;
            Spectrum Ks = Spectrum(1-specular_tint) + specular_tint * Ctint;
            Real r0 = (1.0 - eta) / (1.0 + eta);
            r0 = r0 * r0;
            C0 = specular * r0 * (1 - metallic) * Ks + metallic * base_clr;
        }
        Texture<Spectrum> C0_tex = make_constant_spectrum_texture(C0);

        DisneyDiffuse diffuse = {bsdf.base_color, bsdf.roughness, bsdf.subsurface};
        DisneyClearcoat clearcoat = {bsdf.clearcoat_gloss};
        DisneyMetal metal = {bsdf.base_color, bsdf.roughness, bsdf.anisotropic};
        DisneyGlass glass = {bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta};
        DisneySheen sheen = {bsdf.base_color, bsdf.sheen_tint};
        
        std::vector<Real> weights;
        weightComputation(specular_transmission, metallic, cc, sh, weights);

        result += weights[0] * (*this)(diffuse);
        result += weights[1] * (*this)(clearcoat);
        result += weights[2] * (*this)(metal);
        result += weights[3] * (*this)(glass);
        result += weights[4] * (*this)(sheen);

        return result;
    }
}

Real pdf_sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }

    // Homework 1: implement this!
    if(dot(vertex.geometric_normal, dir_in) <= 0){
        // only glass
        DisneyGlass glass = {bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta};
        return (*this)(glass);
    }
    else{
        Real result;

        Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real cc = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real sh = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);

        DisneyDiffuse diffuse = {bsdf.base_color, bsdf.roughness, bsdf.subsurface};
        DisneyClearcoat clearcoat = {bsdf.clearcoat_gloss};
        DisneyMetal metal = {bsdf.base_color, bsdf.roughness, bsdf.anisotropic};
        DisneyGlass glass = {bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta};
        DisneySheen sheen = {bsdf.base_color, bsdf.sheen_tint};

        std::vector<Real> weights;
        weightComputation(specular_transmission, metallic, cc, sh, weights);
        Real totalWeight = 0;
        for(auto w : weights){ totalWeight += w; }
        for(Real& w : weights){ w /= totalWeight; }

        result += weights[0] * (*this)(diffuse);
        result += weights[1] * (*this)(clearcoat);
        result += weights[2] * (*this)(metal);
        result += weights[3] * (*this)(glass);
        result += weights[4] * (*this)(sheen);

        return result;
    }
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }

    // Homework 1: implement this!
    if(dot(vertex.geometric_normal, dir_in) <= 0){
        // only glass
        DisneyGlass glass = {bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta};
        return (*this)(glass);
    }
    else{
        Real result;

        Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real cc = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real sh = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);

        std::vector<Real> weights;
        weightComputation(specular_transmission, metallic, cc, sh, weights);
        Real totalWeight = 0;
        for(auto w : weights){ totalWeight += w; }
        for(Real& w : weights){ w /= totalWeight; }
        for(int i = 1; i < weights.size(); i++){
            weights[i] += weights[i-1]; // CDF
        }

        // check rnd_param_w in which range
        Real rnd = rnd_param_w;
        int index = weights.size() - 1;
        for (int i = 0; i < weights.size(); i++) {
            if (rnd <= weights[i]) {
                index = i;
                break;
            }
        }
        
        DisneyDiffuse diffuse = {bsdf.base_color, bsdf.roughness, bsdf.subsurface};
        DisneyClearcoat clearcoat = {bsdf.clearcoat_gloss};
        DisneyMetal metal = {bsdf.base_color, bsdf.roughness, bsdf.anisotropic};
        DisneyGlass glass = {bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta};
        DisneySheen sheen = {bsdf.base_color, bsdf.sheen_tint};

        switch(index){
            case 0: 
                return (*this)(diffuse);
            case 1:
                return (*this)(clearcoat);
            case 2:
                return (*this)(metal);
            case 3:
                return (*this)(glass);
            case 4:
                return (*this)(sheen);
        }
    }
}

TextureSpectrum get_texture_op::operator()(const DisneyBSDF &bsdf) const {
    return bsdf.base_color;
}
