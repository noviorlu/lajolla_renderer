#pragma once

// The simplest volumetric renderer: 
// single absorption only homogeneous volume
// only handle directly visible light sources
Spectrum vol_path_tracing_1(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};

    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
    if (!vertex_) {
        // Assume Non Environment Map
        return make_zero_spectrum();
    }
    PathVertex vertex = *vertex_;

    if(!is_light(scene.shapes[vertex.shape_id])) {
        return make_zero_spectrum();
    }

    int medium_id = dot(vertex.shading_frame.n, vertex.geometric_normal) < 0 ? vertex.interior_medium_id : vertex.exterior_medium_id;
    Spectrum sigma_a = get_sigma_a(scene.media[medium_id], vertex.position);

    Real thit = distance(vertex.position, ray.org);

    return exp(-thit * sigma_a) * emission(vertex, -ray.dir, scene);
}

// The second simplest volumetric renderer: 
// single monochromatic homogeneous volume with single scattering,
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_2(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
    
    Real t_hit = infinity<Real>();
    Medium medium = scene.media[scene.camera.medium_id];
    PathVertex vertex;
    if (vertex_) {
        vertex = *vertex_;
        int medium_id = dot(vertex.shading_frame.n, vertex.geometric_normal) < 0 ? vertex.interior_medium_id : vertex.exterior_medium_id;
        medium = scene.media[medium_id];
        t_hit = distance(vertex.position, ray.org);
    }
    
    // since homogeneous medium, doesnot care about position
    Spectrum sigma_a = get_sigma_a(medium, Vector3(0));
    Spectrum sigma_s = get_sigma_s(medium, Vector3(0));
    Spectrum sigma_t = sigma_a + sigma_s;

    Real t = -log(1 - next_pcg32_real<Real>(rng)) / sigma_t[0];
    if (t >= t_hit) {
        // hit surface, account for emission
        if(!is_light(scene.shapes[vertex.shape_id])) {
            return make_zero_spectrum();
        }
        // where (transmittance / transmittance_pdf) 
        //          = exp(-sigma_t * t_hit) / exp(-sigma_t * t_hit) = 1
        return emission(vertex, -ray.dir, scene);
    }
    else{
        // hit volume, account for single scattering
        Vector3 volume_point = ray.org + t * ray.dir;

        // sample a light source
        Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
        Real light_w = next_pcg32_real<Real>(rng);
        Real shape_w = next_pcg32_real<Real>(rng);
        int light_id = sample_light(scene, light_w);
        const Light &light = scene.lights[light_id];
        PointAndNormal point_on_light =
            sample_point_on_light(light, volume_point, light_uv, shape_w, scene);
        

        Spectrum L_distanceSampling = make_zero_spectrum();
        Real p_distanceSampling = 0;
        {
            // compute G, p1, L, light_transmittance
            Real t_light = distance(point_on_light.position, volume_point);
            Vector3 dir_light = normalize(point_on_light.position - volume_point);
            
            Real G = 0;
            Ray shadow_ray{volume_point, dir_light, get_shadow_epsilon(scene), t_light - get_shadow_epsilon(scene)};
            if (!occluded(scene, shadow_ray))
                G = max(abs(dot(dir_light, point_on_light.normal)), Real(0)) / (t_light * t_light);

            Real p_light = light_pmf(scene, light_id) *
                pdf_point_on_light(light, point_on_light, volume_point, scene);

            Spectrum L = Spectrum(0);
            if(p_light > 0 && G > 0)
                L = emission(light, -dir_light, Real(0), point_on_light, scene);

            PhaseFunction& phase_function = get_phase_function(medium);
            Spectrum rho = eval(phase_function, -ray.dir, dir_light);
            
            Spectrum light_transmittance = exp(-t_light * sigma_t);

            // monochromatic medium
            Real sigma_ts = sigma_t[0];

            L_distanceSampling = exp(-t * sigma_ts)* rho * G * light_transmittance * L;
            p_distanceSampling = p_light * exp(-t * sigma_ts) * sigma_ts;
        }
        // return sigma_s * L_distanceSampling / p_distanceSampling;

        // MIS with equiangular sampling
        Spectrum L_equiangularSampling = make_zero_spectrum();
        Real p_equiangularSampling = 0;
        {
            Vector3 dir_a = point_on_light.position - ray.org;
    
            Real a = dot(dir_a, ray.dir) / length_squared(ray.dir);
            Real b = t_hit - a;
            Real D = length(ray.dir * a + ray.org - point_on_light.position);

            Real theta_a = atan(a / D);
            Real theta_b = atan(b / D);

            Real u_rand = next_pcg32_real<Real>(rng);
            t = a + D * tan((theta_a  + theta_b) * u_rand - theta_a);

            volume_point = ray.org + t * ray.dir;

            // compute G, p, L, light_transmittance
            Real t_light = distance(point_on_light.position, volume_point);
            Vector3 dir_light = normalize(point_on_light.position - volume_point);
            
            Real G = 0;
            Ray shadow_ray{volume_point, dir_light, get_shadow_epsilon(scene), t_light - get_shadow_epsilon(scene)};
            if (!occluded(scene, shadow_ray))
                G = max(abs(dot(dir_light, point_on_light.normal)), Real(0)) / (t_light * t_light);

            Real p_light = light_pmf(scene, light_id) *
                pdf_point_on_light(light, point_on_light, volume_point, scene);

            Spectrum L = Spectrum(0);
            if(p_light > 0 && G > 0)
                L = emission(light, -dir_light, Real(0), point_on_light, scene);

            PhaseFunction& phase_function = get_phase_function(medium);
            Spectrum rho = eval(phase_function, -ray.dir, dir_light);
            
            Spectrum light_transmittance = exp(-t_light * sigma_t);

            // monochromatic medium
            Real sigma_ts = sigma_t[0];

            L_equiangularSampling = exp(-t * sigma_ts)* rho * G * light_transmittance * L;
            p_equiangularSampling = p_light * D / ((theta_a + theta_b) * (D*D + t_light * t_light));
        }
        // return sigma_s * L_equiangularSampling / p_equiangularSampling;

        Real p_sum = p_distanceSampling + p_equiangularSampling;
        if(p_sum == 0)
            return make_zero_spectrum();
        
        return sigma_s * (L_distanceSampling + L_equiangularSampling) / p_sum;
    }
}

int update_medium(const PathVertex& isect, const Ray& ray, int current_medium_id){
    int medium_id = current_medium_id;
    if(isect.interior_medium_id != isect.exterior_medium_id){
        if(dot(ray.dir, isect.geometric_normal) > 0){
            medium_id = isect.exterior_medium_id;
        }
        else{
            medium_id = isect.interior_medium_id;
        }
    }
    return medium_id;
}

// The third volumetric renderer (not so simple anymore): 
// multiple monochromatic homogeneous volumes with multiple scattering
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_3(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
    
    int current_medium_id = scene.camera.medium_id;

    Spectrum current_path_throughput = fromRGB(Vector3{1, 1, 1});
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;
    
    const int max_depth = scene.options.max_depth;
    const int rr_depth = scene.options.rr_depth;

    while(true){
        bool scatter = false;
        std::optional<PathVertex> isect_ = intersect(scene, ray, ray_diff);
        
        Real transmittance = Real(1);
        Real trans_pdf = Real(1);
        Real sigma_s = Real(0);

        if(current_medium_id != -1){
            const Medium& current_medium = scene.media[current_medium_id];
            Real sigma_a = get_sigma_a(current_medium, ray.org)[0];
            sigma_s = get_sigma_s(current_medium, ray.org)[0];
            Real sigma_t = sigma_a + sigma_s;

            Real t = -log(1 - next_pcg32_real<Real>(rng)) / sigma_t;

            Real t_hit = infinity<Real>();
            if(isect_){
                t_hit = distance((*isect_).position, ray.org);
            }

            if(t < t_hit){
                scatter = true;
                transmittance = exp(-t * sigma_t);
                trans_pdf = sigma_t * exp(-sigma_t * t);
            }
            else{
                // see zhihu section 4.2
                transmittance = exp(-sigma_t * t_hit);
                trans_pdf = exp(-sigma_t * t_hit);
            }

            ray.org += t * ray.dir;
        }

        current_path_throughput *= (transmittance / trans_pdf);

        if(!scatter && isect_ && is_light(scene.shapes[isect_->shape_id])){
            radiance += current_path_throughput * emission(*isect_, -ray.dir, scene);
        }

        // Reach max bounces
        if(bounces == max_depth - 1 && max_depth != -1) break; 

        // Index-matched surface, pass through
        if(!scatter && isect_ && (*isect_).material_id == -1){
            current_medium_id = update_medium(*isect_, ray, current_medium_id);
            bounces++;
            continue;
        }

        // Scatter YES, sample dir & update path throughput
        if(scatter){
            PhaseFunction& phase_function = get_phase_function(scene.media[current_medium_id]);
            Vector2 phase_func_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};

            std::optional<Vector3> next_dir_ = sample_phase_function(
                phase_function, 
                -ray.dir, 
                phase_func_param_uv
            );
            
            current_path_throughput *= 
                eval(phase_function, -ray.dir, *next_dir_) / 
                pdf_sample_phase(phase_function, -ray.dir, *next_dir_) * 
                sigma_s;

            ray.dir = *next_dir_;
        }
        else{
            break;
        }

        // Russian roulette
        Real rr_prob = 1;
        if(bounces >= rr_depth){
            rr_prob = min(max(current_path_throughput), Real(0.95));
            if(next_pcg32_real<Real>(rng) > rr_prob){
                break;
            }else{
                current_path_throughput /= rr_prob;
            }
        }
        bounces++;
    }

    return radiance;
}


struct LightSampleInfo{
public:
    int light_id;
    PointAndNormal point_on_light;
    Vector3 dir_light;
    Real pdf_light = Real(0);
    Real dist_light = Real(0); // distance from point p to light
    Spectrum L = make_zero_spectrum();
    Real G = Real(0);

    LightSampleInfo(const Scene &scene, pcg32_state &rng, const Vector3& p){
        Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
        Real light_w = next_pcg32_real<Real>(rng);
        Real shape_w = next_pcg32_real<Real>(rng);
        light_id = sample_light(scene, light_w);
        const Light &light = scene.lights[light_id];
        
        point_on_light = sample_point_on_light(light, p, light_uv, shape_w, scene);
        pdf_light = light_pmf(scene, light_id) * pdf_point_on_light(light, point_on_light, p, scene);
        
        dir_light = normalize(point_on_light.position - p);
        dist_light = distance(point_on_light.position, p);
        
        L = emission(light, -dir_light, Real(0), point_on_light, scene);
        G = max(abs(dot(dir_light, point_on_light.normal)), Real(0)) / (dist_light * dist_light);
    }
    
    void Recompute(const Scene &scene, const Vector3& p){
        pdf_light = light_pmf(scene, light_id) * pdf_point_on_light(scene.lights[light_id], point_on_light, p, scene);
        
        dir_light = normalize(point_on_light.position - p);
        dist_light = distance(point_on_light.position, p);

        L = emission(scene.lights[light_id], -dir_light, Real(0), point_on_light, scene);
        G = max(abs(dot(dir_light, point_on_light.normal)), Real(0)) / (dist_light * dist_light);
    }

    void EquiangularSampling(
        pcg32_state &rng, const Vector3& p, const Vector3& dir, const Real& t_hit, const Real& sigma_t,
        Real& t, Real& transmittance, Real& trans_pdf
    ){
        Vector3 dir_a = point_on_light.position - p;

        Real a = dot(dir_a, dir) / length_squared(dir);
        Real b = t_hit - a;
        Real D = length(dir * a + p - point_on_light.position);

        Real theta_a = atan(a / D);
        Real theta_b = atan(b / D);

        Real u_rand = next_pcg32_real<Real>(rng);
        t = a + D * tan((theta_a  + theta_b) * u_rand - theta_a);

        transmittance = exp(-sigma_t * t);
        trans_pdf = D / ((theta_a + theta_b) * (D * D + t * t));
    }
};

inline Spectrum next_event_estimation(
    const Scene &scene, pcg32_state &rng, LightSampleInfo* l_info,
    Vector3 p, const Vector3& dir_in, 
    int current_medium_id, int bounces, bool is_volume, 
    const PathVertex* isect = nullptr // 如果是 surface，需要传入 isect
) {
    if(l_info == nullptr) l_info = new LightSampleInfo(scene, rng, p);
    if(max(l_info->L) == 0 || l_info->L == 0) return make_zero_spectrum();

    /* BSDF or Phase Function Related */
    Spectrum f_or_rho = Spectrum(1);
    Real pdf_f_or_phase = 1;
    
    if (is_volume) {
        // Volume scattering
        if (current_medium_id != -1) {
            const Medium& current_medium = scene.media[current_medium_id];
            PhaseFunction& phase_function = get_phase_function(current_medium);
            f_or_rho = eval(phase_function, dir_in, -l_info->dir_light);
            pdf_f_or_phase = pdf_sample_phase(phase_function, dir_in, -l_info->dir_light);
        }
    } else {
        // Surface scattering
        if (!isect) return make_zero_spectrum(); // 防止错误调用
        const Material& mat = scene.materials[isect->material_id];
        f_or_rho = eval(mat, dir_in, l_info->dir_light, *isect, scene.texture_pool);
        pdf_f_or_phase = pdf_sample_bsdf(mat, dir_in, l_info->dir_light, *isect, scene.texture_pool);
    }
    
    if(max(f_or_rho) == 0) return make_zero_spectrum();

    /* Transmittance Calculation */
    Spectrum T_light = make_const_spectrum(1); // Transmittance from p to light
    Real p_trans_dir = 1; // Probability of multiple transmittance
    
    int shadow_medium_id = current_medium_id;
    int shadow_bounces = 0;

    Ray shadow_ray{p, l_info->dir_light, get_shadow_epsilon(scene), l_info->dist_light - get_shadow_epsilon(scene)};
    RayDifferential shadow_ray_diff = RayDifferential{Real(0), Real(0)};

    while (true) {
        Real next_t = distance(l_info->point_on_light.position, p);
        shadow_ray.org = p;
        shadow_ray.tfar = next_t - get_shadow_epsilon(scene);

        std::optional<PathVertex> shadow_isect_ = intersect(scene, shadow_ray, shadow_ray_diff);
        if (shadow_isect_) {
            const PathVertex& shadow_isect = *shadow_isect_;
            next_t = distance(shadow_isect.position, p);
        }

        if (shadow_medium_id != -1) {
            const Medium& shadow_medium = scene.media[shadow_medium_id];

            Real sigma_a = get_sigma_a(shadow_medium, p)[0];
            Real sigma_s = get_sigma_s(shadow_medium, p)[0];
            Real sigma_t = sigma_a + sigma_s;

            Real trans = exp(-sigma_t * next_t);
            T_light *= trans;
            p_trans_dir *= trans;
        }

        if (!shadow_isect_) break; // Reaches light source
        else { // Blocked by an object
            const PathVertex& shadow_isect = *shadow_isect_;
            if (shadow_isect.material_id >= 0) {
                const Material& mat = scene.materials[shadow_isect.material_id];
                Spectrum bsdf = eval(mat, l_info->dir_light, -l_info->dir_light, shadow_isect, scene.texture_pool);
                Real pdf_bsdf = pdf_sample_bsdf(mat, l_info->dir_light, -l_info->dir_light, shadow_isect, scene.texture_pool);
                T_light *= bsdf / pdf_bsdf;
            }

            shadow_bounces++;
            if (scene.options.max_depth != -1 && bounces + shadow_bounces + 1 >= scene.options.max_depth) 
                return make_zero_spectrum();

            shadow_medium_id = update_medium(shadow_isect, shadow_ray, shadow_medium_id);
            p += shadow_ray.dir * next_t;
        }
    }

    /* Contribution Calculation */
    if (max(T_light) > 0) {
        Spectrum contrib = T_light * l_info->G * f_or_rho * l_info->L / l_info->pdf_light;
        pdf_f_or_phase *= p_trans_dir * l_info->G;

        Real w = (l_info->pdf_light * l_info->pdf_light) / (l_info->pdf_light * l_info->pdf_light + pdf_f_or_phase * pdf_f_or_phase);
        return w * contrib;
    }
    return make_zero_spectrum();
}

// The fourth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// still no surface lighting
Spectrum vol_path_tracing_4(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
                                

    int current_medium_id = scene.camera.medium_id;

    Spectrum current_path_throughput = fromRGB(Vector3{1, 1, 1});
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;

    // Caching information for MIS
    Real dir_pdf = Real(0); // in solid angle measure
    Vector3 nee_p_cache = Vector3(0);
    Real multi_trans_pdf = Real(1);

    bool never_scatter = true;

    const int max_depth = scene.options.max_depth;
    const int rr_depth = scene.options.rr_depth;

    while(true){
        bool scatter = false;
        std::optional<PathVertex> isect_ = intersect(scene, ray, ray_diff);
        
/* Recursive Transmittance into path-throughput from Ray origin to dist-sampled point */
        Real sigma_s = Real(0);
        {
            Real transmittance = Real(1);
            Real trans_pdf = Real(1);
            if(current_medium_id != -1){
                const Medium& current_medium = scene.media[current_medium_id];
                Real sigma_a = get_sigma_a(current_medium, ray.org)[0];
                sigma_s = get_sigma_s(current_medium, ray.org)[0];
                Real sigma_t = sigma_a + sigma_s;

                Real t = -log(1 - next_pcg32_real<Real>(rng)) / sigma_t;

                Real t_hit = infinity<Real>();
                if(isect_){
                    t_hit = distance((*isect_).position, ray.org);
                }

                if(t < t_hit){
                    scatter = true;
                    transmittance = exp(-t * sigma_t);
                    trans_pdf = sigma_t * transmittance;
                }
                else{
                    // see zhihu section 4.2
                    transmittance = exp(-sigma_t * t_hit);
                    trans_pdf = transmittance;
                    t = t_hit + get_intersection_epsilon(scene);
                }

                ray.org += t * ray.dir;
            } else if (isect_) { // 场景globally没有medium
                ray.org = isect_->position + ray.dir * get_intersection_epsilon(scene);
            } else { // 场景globally没有medium, 也没有hit到任何物体, 空白区域
                break;
            }
            current_path_throughput *= (transmittance / trans_pdf);
            multi_trans_pdf *= trans_pdf;
        }

/* If hit Emissive Object, Collect Color with MIS_phaseFunction_weighted */
        if(!scatter && isect_ && is_light(scene.shapes[isect_->shape_id])){
            if(never_scatter) radiance += current_path_throughput * emission(*isect_, -ray.dir, scene);
            else{
                int light_id = get_area_light_id(scene.shapes[isect_->shape_id]);
                const Light& light = scene.lights[light_id];
                PointAndNormal point_on_light{(*isect_).position, (*isect_).geometric_normal};
                // why nee_p_cache? NEE potentailly many bounces ago due to Index-matched surface
                Real pdf_light = light_pmf(scene, light_id) * pdf_point_on_light(scene.lights[light_id], point_on_light, nee_p_cache, scene);

                Vector3 dir_light = normalize(isect_->position - nee_p_cache);
                Real dist_light = distance(nee_p_cache, point_on_light.position);
                
                Real G = max(abs(dot(dir_light, point_on_light.normal)), Real(0)) / (dist_light * dist_light);
                Real pdf_phase = dir_pdf * multi_trans_pdf * G;

                Real w = (pdf_phase * pdf_phase) / (pdf_light * pdf_light + pdf_phase * pdf_phase);
                radiance += w * current_path_throughput * emission(*isect_, -ray.dir, scene);
            }
        }

/* Reach max bounces */
        if(bounces == max_depth - 1 && max_depth != -1) break; 

/* Index-matched surface, pass through */
        if (!scatter && isect_ && (*isect_).material_id == -1) {
            current_medium_id = update_medium(*isect_, ray, current_medium_id);
            bounces++;
            continue;
        }

        // Scatter YES, sample dir & update path throughput
        if(scatter){
            never_scatter = false;
            Spectrum nee_contrib = next_event_estimation(scene, rng, nullptr, ray.org, -ray.dir, current_medium_id, bounces, true);
            radiance += current_path_throughput * sigma_s* nee_contrib;

            const PhaseFunction& phase_function = get_phase_function(scene.media[current_medium_id]);
            Vector2 phase_func_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};

            std::optional<Vector3> next_dir_ = sample_phase_function(
                phase_function, 
                -ray.dir, 
                phase_func_param_uv
            );
            
            dir_pdf = pdf_sample_phase(phase_function, -ray.dir, *next_dir_);
            current_path_throughput *= eval(phase_function, -ray.dir, *next_dir_) / dir_pdf * sigma_s;

            ray.dir = *next_dir_;

            nee_p_cache = ray.org;
            multi_trans_pdf = 1;
        }
        else{
            break;
        }

        // Russian roulette
        Real rr_prob = 1;
        if(bounces >= rr_depth){
            rr_prob = min(max(current_path_throughput), Real(0.95));
            if(next_pcg32_real<Real>(rng) > rr_prob){
                break;
            }else{
                current_path_throughput /= rr_prob;
            }
        }
        bounces++;
    }

    return radiance;
}

// The fifth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing_5(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                        (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    ray.tnear = get_intersection_epsilon(scene);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
                                

    int current_medium_id = scene.camera.medium_id;

    Spectrum current_path_throughput = fromRGB(Vector3{1, 1, 1});
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;

// Caching information for MIS
    Real dir_pdf = Real(0); // in solid angle measure
    Vector3 nee_p_cache = Vector3(0);
    Real multi_trans_pdf = Real(1);

    bool never_sampl_dir = true;

    const int max_depth = scene.options.max_depth;
    const int rr_depth = scene.options.rr_depth;

    while(true){
        bool scatter = false;
        std::optional<PathVertex> isect_ = intersect(scene, ray, ray_diff);
        
/* Recursive Transmittance into path-throughput from Ray origin to dist-sampled point */
        LightSampleInfo* l_info = nullptr;
        Real sigma_s = Real(0);
        Real sigma_t = Real(0);
        Real t = 0;
        Real t_hit = infinity<Real>();
        {
            Real transmittance = Real(1);
            Real trans_pdf = Real(1);
            if(current_medium_id != -1){
                const Medium& current_medium = scene.media[current_medium_id];
                Real sigma_a = get_sigma_a(current_medium, ray.org)[0];
                sigma_s = get_sigma_s(current_medium, ray.org)[0];
                sigma_t = sigma_a + sigma_s;

                Real t = -log(1 - next_pcg32_real<Real>(rng)) / sigma_t;
                if(isect_){
                    t_hit = distance((*isect_).position, ray.org);
                }

                if(t < t_hit){
                    scatter = true;
                    transmittance = exp(-t * sigma_t);
                    trans_pdf = sigma_t * transmittance;
                    ray.org += t * ray.dir;
                    // l_info = new LightSampleInfo(scene, rng, ray.org);
                    // l_info->EquiangularSampling(rng, ray.org, ray.dir, t_hit, sigma_t, t, transmittance, trans_pdf);
                    // ray.org += t * ray.dir;
                    // l_info->Recompute(scene, ray.org);
                }
                else{
                    // see zhihu section 4.2
                    transmittance = exp(-sigma_t * t_hit);
                    trans_pdf = transmittance;
                    t = t_hit;
                    ray.org += t * ray.dir;
                }
            } else if (isect_) { // 场景globally没有medium
                ray.org = isect_->position;
            } else { // 场景globally没有medium, 也没有hit到任何物体, 空白区域
                break;
            }
            current_path_throughput *= (transmittance / trans_pdf);
            multi_trans_pdf *= trans_pdf;
        }

/* If hit Emissive Object, Collect Color with MIS_phaseFunction_weighted */
        if(!scatter && isect_ && is_light(scene.shapes[isect_->shape_id])){
            if(never_sampl_dir) radiance += current_path_throughput * emission(*isect_, -ray.dir, scene);
            else{
                int light_id = get_area_light_id(scene.shapes[isect_->shape_id]);
                const Light& light = scene.lights[light_id];
                PointAndNormal point_on_light{(*isect_).position, (*isect_).geometric_normal};
                // why nee_p_cache? NEE potentailly many bounces ago due to Index-matched surface
                Real pdf_light = light_pmf(scene, light_id) * pdf_point_on_light(scene.lights[light_id], point_on_light, nee_p_cache, scene);

                Vector3 dir_light = normalize(isect_->position - nee_p_cache);
                Real dist_light = distance(nee_p_cache, point_on_light.position);
                
                Real G = max(abs(dot(dir_light, point_on_light.normal)), Real(0)) / (dist_light * dist_light);
                Real pdf_phase = dir_pdf * multi_trans_pdf * G;

                Real w = (pdf_phase * pdf_phase) / (pdf_light * pdf_light + pdf_phase * pdf_phase);
                radiance += w * current_path_throughput * emission(*isect_, -ray.dir, scene);
            }
        }

/* Reach max bounces */
        if(bounces == max_depth - 1 && max_depth != -1) break; 

/* Index-matched surface, pass through */
        if (!scatter && isect_ && (*isect_).material_id == -1) {
            current_medium_id = update_medium(*isect_, ray, current_medium_id);
            bounces++;
            continue;
        }

// Scatter YES, sample dir & update path throughput
        never_sampl_dir = false;
        if(scatter){
            Spectrum nee_contrib = next_event_estimation(scene, rng, l_info, ray.org, -ray.dir, current_medium_id, bounces, true);
            radiance += current_path_throughput * sigma_s * nee_contrib;

            const PhaseFunction& phase_function = get_phase_function(scene.media[current_medium_id]);
            Vector2 phase_func_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};

            std::optional<Vector3> next_dir_ = sample_phase_function(
                phase_function, 
                -ray.dir, 
                phase_func_param_uv
            );
            
            dir_pdf = pdf_sample_phase(phase_function, -ray.dir, *next_dir_);
            current_path_throughput *= eval(phase_function, -ray.dir, *next_dir_) / dir_pdf * sigma_s;

            ray.dir = *next_dir_;

            nee_p_cache = ray.org;
            multi_trans_pdf = 1;
        }
        else{
            const Material &mat = scene.materials[isect_->material_id];
            Spectrum nee_contrib = next_event_estimation(scene, rng, nullptr, ray.org, -ray.dir, current_medium_id, bounces, false, &(*isect_));
            radiance += current_path_throughput * nee_contrib;

            Vector2 bsdf_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            Real bsdf_rnd_param_w = next_pcg32_real<Real>(rng);

            std::optional<BSDFSampleRecord> bsdf_sample_ =
                sample_bsdf(
                    mat, -ray.dir, *isect_, scene.texture_pool,
                    bsdf_rnd_param_uv, bsdf_rnd_param_w
                );
            
            if (!bsdf_sample_) break;
            const BSDFSampleRecord &bsdf_sample = *bsdf_sample_;
            
            dir_pdf = pdf_sample_bsdf(mat, -ray.dir, bsdf_sample.dir_out, *isect_, scene.texture_pool);
            current_path_throughput *= eval(mat, -ray.dir, bsdf_sample.dir_out, *isect_, scene.texture_pool) / dir_pdf;
            
            // refracted, update medium
            if(bsdf_sample.eta != 0) current_medium_id = update_medium(*isect_, ray, current_medium_id);

            ray.dir = bsdf_sample.dir_out;

            nee_p_cache = ray.org;
            multi_trans_pdf = 1;
        }

// Russian roulette
        Real rr_prob = 1;
        if(bounces >= rr_depth){
            rr_prob = min(max(current_path_throughput), Real(0.95));
            if(next_pcg32_real<Real>(rng) > rr_prob){
                break;
            }else{
                current_path_throughput /= rr_prob;
            }
        }
        
        bounces++;
    }

    return radiance;
}






// The final volumetric renderer: 
// multiple chromatic heterogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing(const Scene &scene,
                          int x, int y, /* pixel coordinates */
                          pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
}
