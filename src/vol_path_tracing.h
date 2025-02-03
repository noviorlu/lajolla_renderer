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

int update_medium(const PathVertex& isect, const Ray& ray){
    int medium_id = isect.interior_medium_id;
    if(isect.interior_medium_id != isect.exterior_medium_id){
        medium_id = dot(isect.shading_frame.n, isect.geometric_normal) < 0 ? isect.interior_medium_id : isect.exterior_medium_id;
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
    const int max_bounces = scene.options.max_depth;
    const int rr_depth = scene.options.rr_depth;

    while(true){
        bool isScatter = false;
        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        PathVertex vertex = *vertex_;
        Real t_hit = infinity<Real>();
        if(vertex_){
            t_hit = distance(vertex.position, ray.org);
        }

        Spectrum sigma_a = make_zero_spectrum();
        Spectrum sigma_s = make_zero_spectrum();
        Spectrum sigma_t = make_zero_spectrum();

        if(current_medium_id != -1){
            Medium current_medium = scene.media[current_medium_id];
            sigma_a = get_sigma_a(current_medium, Vector3(0));
            sigma_s = get_sigma_s(current_medium, Vector3(0));
            sigma_t = sigma_a + sigma_s;

            Real t = -log(1 - next_pcg32_real<Real>(rng)) / sigma_t[0];

            if(t < t_hit){
                isScatter = true;
            }
            else{
                t = t_hit;
            }

            // where (transmittance / transmittance_pdf) 
            //      = exp(-sigma_t * t) / (exp(-sigma_t * t) * sigma_t) = 1 / sigma_t 
            current_path_throughput /= sigma_t[0];

            ray.org += t * ray.dir;
        }

        // Reach Surface/ EnvMap, include emission
        if(isScatter == false){
            radiance += current_path_throughput * emission(vertex, -ray.dir, scene);
        }

        // Reach max bounces
        if(bounces == max_bounces - 1 && max_bounces != -1) break; 

        // Index-matched surface, pass through
        if(isScatter == false && !vertex_ && vertex.material_id == -1){
            current_medium_id = update_medium(vertex, ray);
            bounces++;
            continue;
        }

        // Scatter YES, sample dir & update path throughput
        if(isScatter){
            Medium current_medium = scene.media[current_medium_id];

            PhaseFunction& phase_function = get_phase_function(current_medium);

            Vector3 next_dir = sample_phase_function(
                phase_function, 
                -ray.dir, 
                Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng))
            );

            current_path_throughput *= eval(phase_function, -ray.dir, next_dir) / pdf_sample_phase(phase_function, -ray.dir, next_dir) * sigma_s[0];

            ray.dir = next_dir;
        }
        else{
            break;
        }

        // Russian roulette
        Real rr_prob = 1;
        if(bounces >= rr_depth){
            rr_prob = min(current_path_throughput[0], Real(0.95));
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

// The fourth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// still no surface lighting
Spectrum vol_path_tracing_4(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
}

// The fifth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing_5(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
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
