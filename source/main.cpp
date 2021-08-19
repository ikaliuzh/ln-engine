#include <limits>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <vector>
#include "geometry.h"

const size_t MAX_RAY_CAST_RECURSION_DEPTH = 4;

struct Material{
	Material(float r, const vec<4, float>& alb, const vec<3, float>& color, const float& spec)
		: refractive_index(r),
		  albedo(alb),
		  diffuse_color(color),
		  specular_exponent(spec){}
	Material() : albedo({1, 0, 0, 0}), diffuse_color(), specular_exponent(){}

	vec<4, float> albedo; // diffuse, specular, reflection
	vec<3, float> diffuse_color;
	float specular_exponent, refractive_index;

};

struct Light{
	Light(const vec<3, float>& p, float i) : position(p), intensity(i) {}


	vec<3, float> position;
	float intensity;
};

vec<3, float> reflect(const vec<3, float>& L, const vec<3, float>& N){
	return L - 2.f*(L*N)*N;
}

vec<3, float> refract(const vec<3, float>& L, const vec<3, float>& N, const float& refractive_index){
	float cosl = - std::max(-1.f, std::min(1.f, L*N));
	float etal = 1;
	float etar = refractive_index;
	vec<3, float> n = N;
	if (cosl < 0){
		cosl = -cosl;
		std::swap(etal, etar);
		n = -N;
	}
    float eta = etal / etar;
    float k = 1 - eta*eta*(1 - cosl*cosl);

    return k < 0 ? vec<3, float>(0,0,0) : L*eta + n*(eta * cosl - sqrtf(k));
}

struct Sphere {
	vec<3, float> center;
	float radius;
	Material material;

	Sphere(const vec<3, float>& c, const float& r) : center(c), radius(r){}
	Sphere(const vec<3, float>& c, const float& r, const Material& m) : center(c), radius(r), material(m){}


    bool ray_intersect(const vec<3, float> &orig, const vec<3, float> &dir, float &t0) const {
        // assert(dir.norm() == 1.0);
    	vec<3, float> L = center - orig;
        float tca = L*dir;
        float d2 = L*L - tca*tca;
        if (d2 > radius*radius) return false;
        float thc = sqrtf(radius*radius - d2);
        t0 = tca - thc; // distance to the first intersection with the sphere
        float t1 = tca + thc; // distance to the second intersection
        if (t0 < 0) t0 = t1;
        if (t0 < 0) return false; // sphere is behind the scene
        return true;
    }
};

bool scene_intersect(
		const vec<3, float>& orig,
		const vec<3, float>& dir,
		const std::vector<Sphere>& spheres,
		vec<3, float>& hit, // hit point
		vec<3, float>& N, // normal to be calculated from sphere\ray position
		Material& material){ // material to be updated from sphere material

	float spheres_dist = std::numeric_limits<float>::max(); // distance to the closest sphere
	for (size_t i = 0; i < spheres.size(); ++i){
		float dist_i; // distance to the intersection along the ray
		if (spheres[i].ray_intersect(orig, dir, dist_i) && dist_i < spheres_dist){
			spheres_dist = dist_i;
			hit = orig + dir*dist_i;
			N = (hit - spheres[i].center).normalize();
			material = spheres[i].material;
		}
	}
	// finally hit, N, material correspond to those of the closest sphere
	return spheres_dist < 1000;
}

vec<3, float> cast_ray(
		const vec<3, float> &orig,
		const vec<3, float> &dir,
		const std::vector<Sphere>& spheres,
		const std::vector<Light>& lights,
		size_t recursion_depth = 0) {
	vec<3, float> point, N;
	Material material;

	if (recursion_depth>MAX_RAY_CAST_RECURSION_DEPTH || !scene_intersect(orig, dir, spheres, point, N, material)) {
        return vec<3, float>(0.2, 0.7, 0.8);
    }

	float diffuse_light_intensity = 0;
	float specular_light_intensity = 0;

	for (size_t i = 0; i < lights.size(); ++i){
		vec<3, float> light_dir = (lights[i].position - point).normalize();
		float light_dist = (lights[i].position - point).norm();

		vec<3, float> shadow_orig = (light_dir*N < 0)? point - N*1e-3 : point + N*1e-3;
        vec<3, float> shadow_pt, shadow_N;
        Material tmpmaterial;
        if (scene_intersect(shadow_orig, light_dir, spheres, shadow_pt, shadow_N, tmpmaterial) && (shadow_pt-shadow_orig).norm() < light_dist)
            continue;

		diffuse_light_intensity += lights[i].intensity * std::max<float>(0.0, light_dir * N);
		specular_light_intensity += lights[i].intensity * powf(std::max(0.f, reflect(light_dir, N)*dir), material.specular_exponent);
	}

	vec<3, float> reflected_dir = reflect(dir, N);
	vec<3, float> reflected_orig = reflected_dir*N < 0 ? point - N*1e-3 : point + N*1e-3;
	vec<3, float> reflected_color = cast_ray(reflected_orig, reflected_dir, spheres, lights, recursion_depth+1);

	vec<3, float> refracted_dir = refract(dir, N, material.refractive_index).normalize();
	vec<3, float> refracted_orig = refracted_dir*N < 0 ? point - N*1e-3 : point + N*1e-3;
	vec<3, float> refracted_color = cast_ray(refracted_orig, refracted_dir, spheres, lights, recursion_depth+1);

	return material.diffuse_color * diffuse_light_intensity * material.albedo[0] +
			vec<3,float>(1., 1., 1.)*specular_light_intensity * material.albedo[1] +
			reflected_color * material.albedo[2] +
			refracted_color * material.albedo[3];
}



void render(const std::vector<Sphere>& spheres, const std::vector<Light>& lights) {
    const int width    = 1024;
    const int height   = 768;
    const int fov      = M_PI/2.;
    std::vector<vec<3, float>> framebuffer(width*height);

	// #pragma omp parallel for
    for (size_t j = 0; j<height; j++) {
        for (size_t i = 0; i<width; i++) {
            float x =  (2*(i + 0.5)/(float)width  - 1)*tan(fov/2.)*width/(float)height;
            float y = -(2*(j + 0.5)/(float)height - 1)*tan(fov/2.);
            vec<3, float> dir = vec<3, float>(x, y, -1).normalize();
            framebuffer[i+j*width] = cast_ray(vec<3, float>(0,0,0), dir, spheres, lights);
        }
    }

    std::ofstream ofs; // save the framebuffer to file
    ofs.open("./out.ppm");
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (size_t i = 0; i < height*width; ++i) {
        for (size_t j = 0; j<3; j++) {
            vec<3, float> &c = framebuffer[i];
            float max = std::max(c[0], std::max(c[1], c[2]));
            if (max>1) c = c*(1./max);
            ofs << (char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }
    ofs.close();
}

int main() {
    Material      ivory(1.0, vec<4, float>({0.6,  0.3, 0.1, 0.0}), vec<3, float>(0.4, 0.4, 0.3),   50.);
    Material      glass(1.5, vec<4, float>({0.0,  0.5, 0.1, 0.8}), vec<3, float>(0.6, 0.7, 0.8),  125.);
    Material red_rubber(1.0, vec<4, float>({0.9,  0.1, 0.0, 0.0}), vec<3, float>(0.3, 0.1, 0.1),   10.);
    Material     mirror(1.0, vec<4, float>({0.0, 10.0, 0.8, 0.0}), vec<3, float>(1.0, 1.0, 1.0), 1425.);

    std::vector<Sphere> spheres;
    spheres.push_back(Sphere(vec<3, float>(-3,    0,   -16), 2,      ivory));
    spheres.push_back(Sphere(vec<3, float>(-1.0, -1.5, -12), 2, glass));
    spheres.push_back(Sphere(vec<3, float>( 1.5, -0.5, -18), 3, red_rubber));
    spheres.push_back(Sphere(vec<3, float>( 3,   10,   -18), 4,      mirror));

     std::vector<Light>  lights;
     lights.push_back(Light(vec<3, float>(-20, 20,  20), 1.5));
     lights.push_back(Light(vec<3, float>( 30, 50, -25), 1.8));
     lights.push_back(Light(vec<3, float>( 30, 20,  30), 1.7));

    render(spheres, lights);
    return 0;
}
