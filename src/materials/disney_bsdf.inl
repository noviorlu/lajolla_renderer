#include "../microfacet.h"
#include <vector>

void weightComputation(
    Real specular,
    Real specular_transmission,
    Real metallic,
    Real clearcoat,
    Real sheen,
    std::vector<Real>& weights
){
    // diffuse, clearcoat, metal, glass, sheen
    weights.resize(5);
    weights[0] = (1 - metallic) * (1 - specular_transmission);
    weights[1] = 0.25 * clearcoat;
    weights[2] = (1 - specular_transmission * (1 - metallic));
    weights[3] = (1 - metallic) * specular_transmission;
    weights[4] = (1 - metallic) * sheen;
}

Spectrum eval_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    // Homework 1: implement this!
    if(!reflect || dot(frame.n, dir_in) <= 0){
        // only glass
        DisneyGlass glass = {bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta};
        return (*this)(glass);
    }
    else{
        Spectrum base_clr = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real cc = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real sh = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
        Spectrum C0;
        {
            Real l = luminance(base_clr);
            Spectrum Ctint = Spectrum(1);
            if (l > 0) Ctint = base_clr / l;
            Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
            Spectrum Ks = (1 - specular_tint) + specular_tint * Ctint;
            Real r0 = (1.0 - eta) / (1.0 + eta);
            r0 = r0 * r0;
            C0 = specular * r0 * (1 - metallic) * Ks + metallic * base_clr;
        }
        Texture<Spectrum> C0_tex = make_constant_spectrum_texture(C0);

        DisneyDiffuse diffuse_bsdf = {bsdf.base_color, bsdf.roughness, bsdf.subsurface};
        DisneyClearcoat clearcoat_bsdf = {bsdf.clearcoat_gloss};
        DisneyMetal metal_bsdf = {C0_tex, bsdf.roughness, bsdf.anisotropic};
        DisneyGlass glass_bsdf = {bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta};
        DisneySheen sheen_bsdf = {bsdf.base_color, bsdf.sheen_tint};

        Spectrum f_Diffuse = make_zero_spectrum();
        Spectrum f_Clearcoat = make_zero_spectrum();
        Spectrum f_Sheen = make_zero_spectrum();
        Spectrum f_Glass = make_zero_spectrum();
        Spectrum f_Metal = make_zero_spectrum();


        f_Glass = this->operator()(glass_bsdf);
        if (dot(vertex.shading_frame.n, dir_out) > 0)
        {
            f_Metal = this->operator()(metal_bsdf);
            f_Diffuse = this->operator()(diffuse_bsdf);
            f_Clearcoat = this->operator()(clearcoat_bsdf);
            f_Sheen = this->operator()(sheen_bsdf);
        }
        std::vector<Real> weights;
        weightComputation(specular, specular_transmission, metallic, cc, sh, weights);

        Spectrum result = weights[0] * f_Diffuse
                        + weights[1] * f_Clearcoat
                        + weights[2] * f_Metal
                        + weights[3] * f_Glass
                        + weights[4] * f_Sheen;

        return result;
    }
}

Real pdf_sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    // Homework 1: implement this!
    if(!reflect || dot(frame.n, dir_in) <= 0){
        // only glass
        DisneyGlass glass = {bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta};
        return (*this)(glass);
    }
    else{
        Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
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
        weightComputation(specular, specular_transmission, metallic, cc, sh, weights);
        Real totalWeight = 0;
        for(auto w : weights){ totalWeight += w; }
        for(Real& w : weights){ w /= totalWeight; }

        Real diffuse_bsdf_pdf = (*this)(diffuse);
        Real clearcoat_bsdf_pdf = (*this)(clearcoat);
        Real metal_bsdf_pdf = (*this)(metal);
        Real glass_bsdf_pdf = (*this)(glass);
        Real sheen_bsdf_pdf = (*this)(sheen);

        Real result = Real(0);

        result += weights[0] * diffuse_bsdf_pdf;
        result += weights[1] * clearcoat_bsdf_pdf;
        result += weights[2] * metal_bsdf_pdf;
        result += weights[3] * glass_bsdf_pdf;
        result += weights[4] * sheen_bsdf_pdf;

        return result;
    }
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    // Homework 1: implement this!
    if(dot(frame.n, dir_in) <= 0){
        // only glass
        DisneyGlass glass = {bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta};
        return (*this)(glass);
    }
    else{
        Real result;

        Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real cc = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real sh = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);

        std::vector<Real> weights;
        weightComputation(specular, specular_transmission, metallic, cc, sh, weights);
        weights.pop_back();
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

        switch(index){
            case 0: 
                return (*this)(diffuse);
            case 1:
                return (*this)(clearcoat);
            case 2:
                return (*this)(metal);
            case 3:
                return (*this)(glass);
        }
    }
}

TextureSpectrum get_texture_op::operator()(const DisneyBSDF &bsdf) const {
    return bsdf.base_color;
}
