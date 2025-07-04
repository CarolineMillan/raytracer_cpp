/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2018.
 *
 * Do what you like with this code as long as you retain this comment.
 */

#include "scene.h"
#include "transform.h"
#include "polymesh.h"
#include "phong.h"
#include "sphere.h"
#include "directional_light.h"
#include "point_light.h"
#include "photon.h"
#include "kdtree.h"

#include <string>
#include <filesystem>

#include <cmath>

//#pragma message ("RESOURCE_DIR => " RESOURCE_DIR)



Scene::Scene()
{
	object_list = 0;
	light_list = 0;
    causticTree   = nullptr;
    globalTree    = nullptr;
}

Scene::~Scene() {
    // 1) Delete all objects in your linked-list:
    Object *o = object_list;
    while (o) {
        Object *next = o->next;
        delete o;
        o = next;
    }

    // 2) Delete all lights in your linked-list:
    PointLight *pl = light_list;
    while (pl) {
        PointLight *next = pl->next;
        delete pl;
        pl = next;
    }

    // 3) Delete your KD-trees:
    delete causticTree;   // safe if nullptr
    delete globalTree;    // safe if nullptr

    // 4) (Optional) clear vectors if you still use them:
    causticPhotons.clear();
    globalPhotons.clear();
}

// TODO shorten the teapot_box() and test() method
/// A NICE SCENE -- TEAPOT AND SPHERE INSIDE AN OPEN BOX
void Scene::teapot_box() {

	string resourceDir = RESOURCE_DIR;
	
	// The following transform allows 4D homogeneous coordinates to be transformed. It moves the supplied teapot model to somewhere visible.
	Transform *transform = new Transform(1.0f, 0.0f, 0.0f,  0.0f,
			0.0f, 0.0f, 1.0f, -2.7f,
			0.0f, 1.0f, 0.0f, 5.0f,
			0.0f, 0.0f, 0.0f, 1.0f);

	//  Read in the teapot model.
	//PolyMesh *pm = new PolyMesh((char *)"teapot_smaller.ply", transform);

	std::string teapotPath = resourceDir + "/teapot_smaller.ply";
    PolyMesh *pm = new PolyMesh((char *)teapotPath.c_str(), transform);

	//creates Phong surface illumination model for polymesh

	float scaling = 1.0/M_PI;

	// rgb(244, 250, 252)
	//Phong glass; 
	glass.ambient = Colour(244.0/255.0, 250.0/255.0, 252.0/255.0, 255.0/255.0);
	glass.diffuse = Colour(244.0/255.0, 250.0/255.0, 252.0/255.0, 255.0/255.0);
	glass.specular = Colour(255.0/255.0, 255.0/255.0, 255.0/255.0, 255.0/255.0);
	glass.BRDF_d = glass.diffuse;
	glass.BRDF_d.scale(scaling);
	glass.power = 40.0f;

    //Phong bp1; 
	// rgb(211, 141, 255)
	mat_pm.ambient = Colour(207.0/255.0, 207.0/255.0, 207.0/255.0, 255.0/255.0);
	mat_pm.diffuse = Colour(207.0/255.0, 207.0/255.0, 207.0/255.0, 255.0/255.0);
	mat_pm.specular = Colour(255.0/255.0, 255.0/255.0, 255.0/255.0, 255.0/255.0);
	mat_pm.BRDF_d = mat_pm.diffuse;
	mat_pm.BRDF_d.scale(scaling);
	mat_pm.power = 40.0f;

	// rgb(255, 252, 230)
	//Phong bp2;
	mat_wall2.ambient = Colour(255.0/255.0, 252.0/255.0, 230.0/255.0, 255.0/255.0);
	mat_wall2.diffuse = Colour(255.0/255.0, 252.0/255.0, 230.0/255.0, 255.0/255.0);
	mat_wall2.specular = Colour(255.0/255.0, 255.0/255.0, 255.0/255.0, 255.0/255.0);
	mat_wall2.BRDF_d = mat_wall2.diffuse;
	mat_wall2.BRDF_d.scale(scaling);
	mat_wall2.power = 40.0f;

	// rgb(247, 198, 198)
	//Phong bp3;
	mat_wall3.ambient = Colour(247.0/255.0, 198.0/255.0, 198.0/255.0, 255.0/255.0);
	mat_wall3.diffuse = Colour(247.0/255.0, 198.0/255.0, 198.0/255.0, 255.0/255.0);
	mat_wall3.specular = Colour(255.0/255.0, 255.0/255.0, 255.0/255.0, 255.0/255.0);
	mat_wall3.BRDF_d = mat_wall3.diffuse;
	mat_wall3.BRDF_d.scale(scaling);
	mat_wall3.power = 40.0f;

	// rgb(199, 247, 198)
	//Phong bp4;
	mat_wall4.ambient = Colour(199.0/255.0, 247.0/255.0, 198.0/255.0, 255.0/255.0);
	mat_wall4.diffuse = Colour(199.0/255.0, 247.0/255.0, 198.0/255.0, 255.0/255.0);
	mat_wall4.specular = Colour(255.0/255.0, 255.0/255.0, 255.0/255.0, 255.0/255.0);
	mat_wall4.BRDF_d = mat_wall4.diffuse;
	mat_wall4.BRDF_d.scale(scaling);
	mat_wall4.power = 40.0f;

	// rgb(198, 242, 247)
	//Phong bp5;
	mat_wall5.ambient = Colour(198.0/255.0, 242.0/255.0, 247.0/255.0, 255.0/255.0);
	mat_wall5.diffuse = Colour(198.0/255.0, 242.0/255.0, 247.0/255.0, 255.0/255.0);
	mat_wall5.specular = Colour(255.0/255.0, 255.0/255.0, 255.0/255.0, 255.0/255.0);
	mat_wall5.BRDF_d = mat_wall5.diffuse;
	mat_wall5.BRDF_d.scale(scaling);
	mat_wall5.power = 40.0f;

	// rgb(255, 176, 249)
	//Phong bp6;
	mat_wall6.ambient = Colour(255.0/255.0, 176.0/255.0, 249.0/255.0, 255.0/255.0);
	mat_wall6.diffuse = Colour(255.0/255.0, 176.0/255.0, 249.0/255.0, 255.0/255.0);
	mat_wall6.specular = Colour(255.0/255.0, 255.0/255.0, 255.0/255.0, 255.0/255.0);
	mat_wall6.BRDF_d = mat_wall6.diffuse;
	mat_wall6.BRDF_d.scale(scaling);
	mat_wall6.power = 40.0f;

	// rgb(180, 169, 245)
	//Phong bp7;
	mat_wall7.ambient = Colour(180.0/255.0, 169.0/255.0, 245.0/255.0, 255.0/255.0);
	mat_wall7.diffuse = Colour(180.0/255.0, 169.0/255.0, 245.0/255.0, 255.0/255.0);
	mat_wall7.specular = Colour(255.0/255.0, 255.0/255.0, 255.0/255.0, 255.0/255.0);
	mat_wall7.BRDF_d = mat_wall7.diffuse;
	mat_wall7.BRDF_d.scale(scaling);
	mat_wall7.power = 40.0f;

    pm->material = &mat_pm;
    pm->material->transparent = false;
    pm->material->reflective = true;
	pm->material->eta = 1.5; //glass refractive index
	pm->material->kr = 0.4;
	pm->material->kt = 0.95;

    Vertex v;
	v.x = 1.0f;
	v.y = 1.0f;
	v.z = 3.0f;
    Sphere *sphere = new Sphere(v, 1.0f);

    sphere->material = &glass;
    sphere->material->transparent = true;
    sphere->material->reflective = true;
	sphere->material->eta = 1.5; //glass refractive index
	sphere->material->kr = 0.4;

    pm->next = sphere;

	/// put the scene in a box 

	string fl_path = resourceDir + "/square.ply";
	PolyMesh *fl = new PolyMesh((char *)fl_path.c_str());

	string ce_path = resourceDir + "/ceiling.ply";
	PolyMesh *ce = new PolyMesh((char *)ce_path.c_str());

	string w1_path = resourceDir + "/wall1.ply";
	PolyMesh *w1 = new PolyMesh((char *)w1_path.c_str());

	string w2_path = resourceDir + "/wall2.ply";
	PolyMesh *w2 = new PolyMesh((char *)w2_path.c_str());

	string w3_path = resourceDir + "/wall3.ply";
	PolyMesh *w3 = new PolyMesh((char *)w3_path.c_str());

	string w4_path = resourceDir + "/wall4.ply";
	PolyMesh *w4 = new PolyMesh((char *)w4_path.c_str());

	
	fl->material = &mat_wall2;
    fl->material->transparent = false;
    fl->material->reflective = false;

	ce->material = &mat_wall3;
    ce->material->transparent = false;
    ce->material->reflective = false;

	w1->material = &mat_wall4;
    w1->material->transparent = false;
    w1->material->reflective = false;

	w2->material = &mat_wall5;
    w2->material->transparent = false;
    w2->material->reflective = false;
	
	w3->material = &mat_wall6;
    w3->material->transparent = false;
    w3->material->reflective = false;

	w4->material = &mat_wall7;
    w4->material->transparent = false;
    w4->material->reflective = false;

	sphere->next = fl;
	fl->next = ce;
	ce->next = w1;
	w1->next = w2;
	w2->next = w3;
	w3->next = w4;
	w4->next = nullptr;

	//creates a light source
	Vertex v1 = Vertex(-1.0, 1.0, -1.0);
	Colour c = Colour(1.0f, 1.0f, 1.0f, 1.0f);
	Vector d = Vector(0.01f, 0.01f, 1.0f);

	Vertex sampling = v;
	sampling.sub(v1);

	cout << "sampling: " << sampling.x << ", " << sampling.y << ", " << sampling.z;



	PointLight *pl = new PointLight(v1, c, d);

	// define object_list and light_list
	//object_list = pm;
	//object_list->next = sphere;
	object_list = pm;
	object_list->next = sphere;
	light_list = pl;
	light_list->next = nullptr;
}


void Scene::test() {
	
	/// A SIMPLE SCENE -- SPHERE, BACK WALL, ONE SIDE WALL
	// rgb(211, 141, 255)
	// rgb(244, 250, 252)
	//Phong bp1; 
	mat_sphere = Phong();
	mat_sphere.ambient = Colour(244.0/255.0, 250.0/255.0, 252.0/255.0, 255.0/255.0);
	mat_sphere.diffuse = Colour(244.0/255.0, 250.0/255.0, 252.0/255.0, 255.0/255.0);
	mat_sphere.specular = Colour(255.0/255.0, 255.0/255.0, 255.0/255.0, 255.0/255.0);
	mat_sphere.power = 40.0f;
	mat_sphere.transparent = true;
    mat_sphere.reflective = true;
	mat_sphere.eta = 1.5; //glass refractive index
	mat_sphere.kr = 0.4;
	//mat_sphere.BRDF_s = Colour(0.0, 0.0, 0.0, 0.0);
	mat_sphere.BRDF_d = Colour(0.0, 0.0, 0.0, 0.0);

	//p1 = bp1;

	// rgb(199, 247, 198)
	//Phong bp4;
	mat_wall1.ambient = Colour(199.0/255.0, 247.0/255.0, 198.0/255.0, 255.0/255.0);
	mat_wall1.diffuse = Colour(199.0/255.0, 247.0/255.0, 198.0/255.0, 255.0/255.0);
	mat_wall1.specular = Colour(255.0/255.0, 255.0/255.0, 255.0/255.0, 255.0/255.0);
	mat_wall1.power = 40.0f;

	//p2 = bp4;

	// rgb(255, 176, 249)
	//Phong bp6;
	mat_wall3.ambient = Colour(255.0/255.0, 176.0/255.0, 249.0/255.0, 255.0/255.0);
	mat_wall3.diffuse = Colour(255.0/255.0, 176.0/255.0, 249.0/255.0, 255.0/255.0);
	mat_wall3.specular = Colour(255.0/255.0, 255.0/255.0, 255.0/255.0, 255.0/255.0);
	mat_wall3.power = 40.0f;

	//p3 = bp6;

	Vertex v;
	v.x = -1.0f;
	v.y = -1.0f;
	v.z = 5.0f;
    Sphere *sphere = new Sphere(v, 1.0f);

    sphere->material = &mat_sphere;
	/*
    sphere->material->transparent = true;
    sphere->material->reflective = true;
	sphere->material->eta = 1.5; //glass refractive index
	sphere->material->kr = 0.4;
	sphere->material->BRDF_s = Colour(0.0, 0.0, 0.0, 0.0);
	sphere->material->BRDF_d = Colour(0.0, 0.0, 0.0, 0.0);
*/
	PolyMesh *w1 = new PolyMesh((char *)"wall1.ply");

	w1->material = &mat_wall1;
    w1->material->transparent = true;
    w1->material->reflective = true;
	//w1->material->BRDF_s = Colour(0.0, 0.0, 0.0, 0.0);
	w1->material->BRDF_d = Colour(1.0, 1.0, 1.0, 1.0);
	w1->material->eta = 1.5; //glass refractive index
	w1->material->kr = 0.4;

	PolyMesh *w3 = new PolyMesh((char *)"wall3.ply");
	
	w3->material = &mat_wall3;
    w3->material->transparent = false;
    w3->material->reflective = false;
	//w3->material->BRDF_s = Colour(0.0, 0.0, 0.0, 0.0);
	w3->material->BRDF_d = Colour(1.0, 1.0, 1.0, 1.0);

	object_list = sphere;
	object_list->next = w1;
	w1->next = w3;
	w3->next = nullptr;

	std::cout << "sphere: " << object_list << std::endl;
	std::cout << "w1: " << object_list->next << std::endl;
	std::cout << "w3: " << object_list->next->next << std::endl;

	//creates a light source
    DirectionalLight *dl = new DirectionalLight(Vector(0.01f, 0.01f, 1.0f),Colour(1.0f, 1.0f, 1.0f, 1.0f)); 

	Vertex pt = Vertex(-1.0, -1.0, -1.0);
	Colour c = Colour(1.0, 1.0, 1.0, 1.0);
	Vector dir = Vector(0.1, 0.1, 0.1);
	light_list = new PointLight(pt, c, dir);
	light_list->next = nullptr; 

	std::cout << "light_list->direction: " << light_list->direction.x << ", " << light_list->direction.y << ", " << light_list->direction.z << std::endl;
			
}


Colour Scene::background_colour(float depth) {
	Colour colour = Colour(163.0/255.0, 249.0/255.0, 255.0/255.0, 1.0f);
	// I like rgb(163, 249, 255)
	depth = 7.0f;
	return colour;
}


/// best_hit returns the closest intersection between the ray and the list of objects
void Scene::object_intersection(Ray ray, Hit &best_hit) {
	Object *obj = object_list;
	best_hit.flag = false;

	while(obj != 0) {

		Hit obj_hit;
		obj_hit.flag=false;

		obj->intersection(ray, obj_hit);
		
		// if we have an intersection and it's in front of the camera
		if (obj_hit.flag && obj_hit.t > 0.0f) {
			// if this is the first hit or closest hit
			if (best_hit.flag == false || obj_hit.t < best_hit.t) best_hit = obj_hit;
		}
		obj = obj->next;
	}
	return;
}


void Scene::point_light_intersection(Ray ray, PointLight*& pl, float &depth, bool &flag) {
	
	PointLight *light = light_list;
	pl = nullptr;
	flag = false;

	while(light != 0) {

		Vector toLight;
		toLight.x = light->point.x - ray.position.x;
		toLight.y = light->point.y - ray.position.y;
		toLight.z = light->point.z - ray.position.z;
		float t = toLight.dot(ray.direction);

		if (t > 0) {

			Vector closest;
			closest.x = ray.at(t).x - light->point.x;
			closest.y = ray.at(t).y - light->point.y;
			closest.z = ray.at(t).z - light->point.z;

			// 1e-4f is squared because we square the thing we're comparing it to
			if (closest.len_sqr() < (1e-4f)*(1e-4f) && t < depth) {
				depth = t;
				pl = light;
				flag = true;
			}
		}

		light = light->next;
	}
	return;
}



Ray Scene::get_shadow_ray(Vector ldir, Hit &best_hit) {
	Ray shadow_ray;
	shadow_ray.direction.x = -ldir.x;
	shadow_ray.direction.y = -ldir.y;
	shadow_ray.direction.z = -ldir.z;

	shadow_ray.position.x = best_hit.position.x + (0.0001f * shadow_ray.direction.x);
	shadow_ray.position.y = best_hit.position.y + (0.0001f * shadow_ray.direction.y);
	shadow_ray.position.z = best_hit.position.z + (0.0001f * shadow_ray.direction.z);
	return shadow_ray;
}


// takes a ray and a hit and outputs a refracted ray from the hit point. Also sets hit material's kr and kt values.
void Scene::refract_ray(const Ray &incoming, Hit &hit, Ray &refracted, bool &total_internal_reflection) {

	float eta = hit.what->material->eta;

	Vector n = hit.normal;
	n.normalise();
	Vector d = incoming.direction;
	d.normalise();
	float cos_theta_i = fabs(n.dot(-1*d));

	// test here for total_internal_reflection
	float cos_theta_t = sqrt(1.0 - ( 1.0 / (pow(eta,2.0) ) ) * ( 1.0 - pow(cos_theta_i,2.0) ) );

	if(cos_theta_t >= 0) {

		refracted.direction = (1/eta)*incoming.direction - (cos_theta_t - (1/eta)*cos_theta_i)*hit.normal;

		refracted.position.x = hit.position.x + 0.001*refracted.direction.x;
		refracted.position.y = hit.position.y + 0.001*refracted.direction.y;
		refracted.position.z = hit.position.z + 0.001*refracted.direction.z;

		//fresnel equations
		float rpar = (eta*cos_theta_i - cos_theta_t)/(eta*cos_theta_i + cos_theta_t);
		float rper = (cos_theta_i - eta*cos_theta_t)/(cos_theta_i + eta*cos_theta_t);
		hit.what->material->kr = 0.5*(rpar*rpar + rper*rper); 
		hit.what->material->kt = 1.0-hit.what->material->kr;
	}
	else {
		total_internal_reflection = true;
	}
}


// takes a ray and a hit and outputs a reflected ray from the hit point
void Scene::reflect_ray(const Ray &incoming, Hit &hit, Ray &reflected) {

	Vector n = hit.normal;
	Vector dir = incoming.direction;

	//calculate reflected eye ray
	reflected.direction = dir - 2.0*(n.dot(dir))*n;

	reflected.position.x = hit.position.x + 0.001*reflected.direction.x;
	reflected.position.y = hit.position.y + 0.001*reflected.direction.y;
	reflected.position.z = hit.position.z + 0.001*reflected.direction.z;
}


Colour Scene::get_shadow_colour(Ray ray, Hit best_hit, int ref_limit) {
	
	// get a handle on the first light -- check this
	PointLight *light = light_list;
	Hit shadow_hit;
	shadow_hit.flag = false;	
	Colour colour = Colour();

	// check all the light sources for shadow rays
	// this gives us direct lighting L_d
	while (light != (Light *)0) {
		// check for a hit between best_hit and light source
		Vector viewer;
		Vector ldir;

		//viewing ray position
		viewer.x = -best_hit.position.x; 
		viewer.y = -best_hit.position.y;
		viewer.z = -best_hit.position.z;
		viewer.normalise();
		//is best_hit hit by the light?
		bool lit;
		lit = light->get_direction(best_hit.position, ldir); 

		if(ldir.dot(best_hit.normal)>0)	lit=false; //light is facing wrong way.

		if(lit) {
			// create a "shadow ray" from best_hit to the light source
			// if it intersects anything, then the point is in shadow
			
			Ray shadow_ray = get_shadow_ray(ldir, best_hit);

			auto temp = best_hit.position;
			float dist = temp.distance(light->point);

			object_intersection(shadow_ray, shadow_hit);

			//there's a shadow so no lighting, if  the hit is closer than the point light
			if(shadow_hit.flag && (shadow_hit.t < dist)) lit = false;
		}

		//do colour
		if (lit) {
			Colour intensity;
			Colour scaling;

			light->get_intensity(best_hit.position, scaling);
			best_hit.what->material->compute_light_colour(viewer, best_hit.normal, ldir, intensity);
			intensity.scale(scaling);

			//add intensity to the colour
			colour.add(intensity);
			}
		
		// move on to the next light
		light = light->next;
	}
	return colour;
}


Colour Scene::get_refraction_colour(Ray ray, Hit hit, int ref_limit) {

	Colour colour = Colour();
	ref_limit -= 1;
	if(ref_limit<0)	return colour;

	Ray refracted;
	bool total_internal_reflection = false;
	

	refract_ray(ray, hit, refracted, total_internal_reflection);

	if (!total_internal_reflection) {

		float trans_depth = hit.t;
		Hit h = Hit();

		//raytrace the refracted ray
		raytrace(refracted, colour, trans_depth, ref_limit, h);

		//scale colour by kt
		colour.scale(hit.what->material->kt);
	}
	return colour;
}


Colour Scene::get_reflection_colour(Ray ray, Hit hit, int ref_limit) {

	Colour colour = Colour();
	ref_limit -= 1;
	if(ref_limit<0)	return colour;

	Ray reflected;
	Hit rhit;

	//tests if there is an object to be reflected in the surface
	reflect_ray(ray, hit, reflected);
	object_intersection(reflected, rhit); 

	if(rhit.flag) {

		float ref_depth = hit.t;
		Hit h = Hit();

		//get colour of reflection 
		raytrace(reflected, colour, ref_depth, ref_limit, h); 

		//scale colour by kr
		colour.scale(hit.what->material->kr);
	}
	return colour;
}

// 
Colour Scene::gather_diffuse(const Hit hit, const vector<Photon*> globalNeighbours) {

	Colour diffuse = Colour();
	float r=0, r_min = 0, r_max = 0; //get radius of the sphere from max distance between best_hit.position and rest of nodes
	float d=0, dA=0; 

	for (Photon* p : globalNeighbours) {

		// clamp to avoid negative colour values
		float  f = -1.0*p->direction.dot(hit.normal);
		float dot = max(0.0f, f);

		Colour temp = p->BRDF_d;
		temp.scale(p->intensity);
		temp.scale(dot);
		diffuse.add(temp);

		//find max distance
		d = hit.position.distance(p->position);
		if (d > r) r = d;
	}

	if (r>0) {
		float r_clamped = max(0.05f, min(r, 2.0f));
		dA = M_PI * r_clamped * r_clamped; //area of disc associated with sphere
		//if (dA > 190.0) cout << "dA: " << dA;
		float scaling = (1.0/dA);
		diffuse.scale(scaling);
	}

	//cout << "diffuse: " << diffuse.r << ", " << diffuse.g << ", " << diffuse.b << endl;
				
	return diffuse;
}


/// TODO shorten this method 
Colour Scene::gather_diffuse_reflection(Ray ray, Hit best_hit, vector<Photon*> globalNeighbours) {
	// get rays -- TODO look into emplace_back() instead
	vector<Ray> current_rays;
	Vertex origin = best_hit.position;
	
	for (Photon* p : globalNeighbours) {
		Ray reflected;
		reflect_ray(ray, best_hit, reflected);
		//Ray new_r = Ray(best_hit.position, p->direction);
		current_rays.push_back(reflected);
	}
	//store the intensity at this point, multiply by kr of the object

	//use importance sampling to decide which to follow and then resample at the next hit point (reflection)
	Colour diffuse_reflection = Colour();
	//bounce limit for photons
	int bounce_limit = 3;
	while (bounce_limit > 0 && !current_rays.empty()) {

		bounce_limit--;
		
		vector<Ray> nxt_rays;

		//ideally for imp_rays
		for (Ray &bounce_ray : current_rays) {
		
			Hit nxt_hit;
			object_intersection(bounce_ray, nxt_hit); //get next hit point of each ray

			if(!nxt_hit.flag) continue;

			//do k nearest neighbours sampling again, but scale each reflection by kr (BRDF?) and add them all together
			vector<Photon*> globalNeighbours;
			globalTree->kNearest(nxt_hit.position, 5, globalNeighbours);

			Colour d = gather_diffuse(nxt_hit, globalNeighbours);
			// FIXME somewhere the kr value is being set to -2.something, I've hardcoded a kr value here for now
			//cout << "kr: " << nxt_hit.what->material->kr;
			float f = 0.4;
			d.scale(f); //nxt_hit.what->material->kr);
			diffuse_reflection.add(d);

			// get next rays, if you loop over all globalNeighbours the number of rays increases exponentially
			// so only add one here
			for (Photon* p : globalNeighbours) {
				Ray reflected;
				reflect_ray(bounce_ray, best_hit, reflected);	
				Ray new_r = Ray(nxt_hit.position, p->direction);
				nxt_rays.push_back(reflected);
				break;
			}
		}

		if (!current_rays.empty()) {
			//cout << "no rays: " << current_rays.size() << endl;
			float no_rays = float(current_rays.size());
			float inv = 1.0f / no_rays;
			diffuse_reflection.scale(inv);
		}
		//colour.add(diffuse_reflection);
		current_rays.swap(nxt_rays);
	}
	return diffuse_reflection;
}


Colour Scene::compute_colour(Ray ray, Hit best_hit, float &depth, int ref_limit) {

	Colour colour = Colour();

	// get colour of the material we've hit
	best_hit.what->material->compute_base_colour(colour);

	// get the depth of the object
	depth = best_hit.t;

	// get shadow colour here
	Colour shadow_colour = get_shadow_colour(ray, best_hit, ref_limit);
	colour.add(shadow_colour);

	// this is L_s, specular
	if(best_hit.what->material->transparent)
	{
		Colour refraction = get_refraction_colour(ray, best_hit, ref_limit);
		colour.add(refraction);
	}

	if(best_hit.what->material->reflective)	{
		Colour reflection = get_reflection_colour(ray, best_hit, ref_limit);
		colour.add(reflection);
	}

	// this is L_c, caustic
	if (causticTree) {
		vector<Photon*> causticNeighbours;
		causticTree->kNearest(best_hit.position, 50, causticNeighbours);
		Colour caustic = gather_diffuse(best_hit, causticNeighbours);
		colour.add(caustic);
	}
	
	// this is L_d, diffuse
	if (globalTree) {

		vector<Photon*> globalNeighbours;	
		globalTree->kNearest(best_hit.position, 5, globalNeighbours);

		Colour diffuse = gather_diffuse(best_hit, globalNeighbours);
		colour.add(diffuse);

		Colour diffuse_reflection = gather_diffuse_reflection(ray, best_hit, globalNeighbours);
		colour.add(diffuse_reflection);
	}
	return colour;
}


/// returns the depth and colour of the closest intersection between ray and the list of objects via classic Whitted-style ray tracing with global and caustic photon mapping 
void Scene::raytrace(Ray ray, Colour &colour, float &depth, int ref_limit, Hit &best_hit) {

	if (ref_limit < 0) return;

	// check for point light intersection, gives us L_e value 
	// pl returns the exact point light hit, so you do need it as an input
	PointLight*pl;
	bool flag = false;

	point_light_intersection(ray, pl, depth, flag);

	if (flag) {
		// we've intersected a light source, so return the light colour
		colour = pl->intensity;
		return;
	}

	object_intersection(ray, best_hit);

	// if we found a primitive then compute the colour we should see
	if(best_hit.flag) colour = compute_colour(ray, best_hit, depth, ref_limit); 
	else colour = background_colour(depth);
	return;
}


//doesn't emit towards specular surfaces, it emits in random directions then checks if the intersection is specular
//traces a single photon through the scene, saves them in causticPhotons and globalPhotons
void Scene::photon_trace(Photon *photon, int ref_limit) {

	bool saw_specular = false;
	Ray photon_ray;
	photon->ray(photon_ray);

	while (ref_limit > 0) {
		
		// put the contents of this loop into a sub function
		ref_limit--;
		Hit best_hit;
		object_intersection(photon_ray, best_hit);

		if (!best_hit.flag) break;

		photon->position = best_hit.position;
		bool is_specular = best_hit.what->material->reflective || best_hit.what->material->transparent;

		if (!is_specular) {
			// add to global and check for caustic
			photon->position = best_hit.position;
			photon->BRDF_d = best_hit.what->material->BRDF_d;
			globalPhotons.push_back(photon);

			if (saw_specular) {
				causticPhotons.push_back(photon);
				// reset the caustic photon
				saw_specular = false;
			}
			//photon->g_russian_roulette(best_hit);
			break;
		} else {
			saw_specular = true;
			photon->c_russian_roulette(best_hit);
			Ray new_ray;

			if (photon->reflected) {
				// change the intensity -- here or in c_russian_roulette?
				// don't want to kill a specular photon in russian roulette, 
				// it only terminates if it hits a diffuse surface, or the bounce limit runs out
				reflect_ray(photon_ray, best_hit, new_ray);
				photon_ray = new_ray;
			}
			if (photon->transmitted) {
				bool tir = false;
                refract_ray(photon_ray, best_hit, new_ray, tir);
				photon_ray = new_ray;
				if (tir) {
                    // if total internal reflection, treat as mirror
                    reflect_ray(photon_ray, best_hit, new_ray);
					photon_ray = new_ray;
                }
			}
		}
	}
}


/// creates caustic and global photon maps, saves them in causticTree and globalTree
void Scene::create_photon_maps() {

	causticPhotons.clear();
	globalPhotons.clear();

	PointLight *light = light_list;
	int no_of_photons = 100000;
	for (int n = 0; n < no_of_photons; n++) {
		while (light != (PointLight*) 0) {
			//trace a photon through the scene and store all the hits in a kd tree - how to go through specular objects only?
			Photon *photon = new Photon(*light, no_of_photons);
			if (n < no_of_photons/5) {
				//Vector temp = Vector(0.8, 0.8, 2.8);
				Vector temp = Vector(2.0, 0.0, 4.0);
				temp.normalise();
				photon->direction = temp;
			}
			photon_trace(photon, 5);
			light = light->next;
		}
		light = light_list;
		if (n % (no_of_photons/1000) == 0) { 
			std::cerr << "+" << flush;
		}
	}

	std::cout << "causticPhotons.size: " << causticPhotons.size() << std::endl;
	std::cout << "globalPhotons.size: " << globalPhotons.size() << std::endl;

	if (!causticPhotons.empty()) {
		delete causticTree;                        // in case it existed
		causticTree = new KDTree(causticPhotons);  // build once
	}

	if (!globalPhotons.empty()) {
		delete globalTree;                        // in case it existed
		globalTree = new KDTree(globalPhotons);  // build once
	}
}

