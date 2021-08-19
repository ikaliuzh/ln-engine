#include <limits>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <vector>
#include "geometry.h"

struct Material{
	Material(const vec<3, float>& color) : diffuse_color(color) {}
	Material() : diffuse_color(){}

	vec<3, float> diffuse_color;
};

struct Light{
	Light(const vec<3, float>& p, float i) : position(p), intensity(i) {}


	vec<3, float> position;
	float intensity;
};

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
	return spheres_dist < 100;
}

vec<3, float> cast_ray(
		const vec<3, float> &orig,
		const vec<3, float> &dir,
		const std::vector<Sphere>& spheres,
		const std::vector<Light>& lights) {
	vec<3, float> point, N;
	Material material;

	if (!scene_intersect(orig, dir, spheres, point, N, material)) {
        return vec<3, float>(0.2, 0.7, 0.8);
    }

	float diffuse_light_intensity = 0;

	for (size_t i = 0; i < lights.size(); ++i){
		vec<3, float> light_dir = (lights[i].position - point).normalize();
		diffuse_light_intensity += lights[i].intensity * std::max<float>(0.0, light_dir * N);
	}

    return material.diffuse_color * diffuse_light_intensity;
}



void render(const std::vector<Sphere>& spheres, const std::vector<Light>& lights) {
    const int width    = 1024;
    const int height   = 768;
    const int fov      = M_PI/2.;
    std::vector<vec<3, float>> framebuffer(width*height);

	#pragma omp parallel for
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
            ofs << (char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }
    ofs.close();
}

int main() {
    Material      ivory(vec<3, float>(0.4, 0.4, 0.3));
    Material red_rubber(vec<3, float>(0.3, 0.1, 0.1));

    std::vector<Sphere> spheres;
    spheres.push_back(Sphere(vec<3, float>(-3,    0,   -16), 2,      ivory));
    spheres.push_back(Sphere(vec<3, float>(-1.0, -1.5, -12), 2, red_rubber));
    spheres.push_back(Sphere(vec<3, float>( 1.5, -0.5, -18), 3, red_rubber));
    spheres.push_back(Sphere(vec<3, float>( 7,    5,   -18), 4,      ivory));

    std::vector<Light> lights;
    lights.push_back(Light(vec<3, float>(-10, 10, -2), 2.0));

    render(spheres, lights);
    return 0;
}
