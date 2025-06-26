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

#include <cmath>


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

/*
void Scene::teapot_box() {
	
	/// A NICE SCENE -- TEAPOT AND SPHERE INSIDE AN OPEN BOX

	// The following transform allows 4D homogeneous coordinates to be transformed. It moves the supplied teapot model to somewhere visible.
	Transform *transform = new Transform(1.0f, 0.0f, 0.0f,  0.0f,
			0.0f, 0.0f, 1.0f, -2.7f,
			0.0f, 1.0f, 0.0f, 5.0f,
			0.0f, 0.0f, 0.0f, 1.0f);

	//  Read in the teapot model.
	PolyMesh *pm = new PolyMesh((char *)"teapot_smaller.ply", transform);
	//objects->next = pm;

	//creates Phong surface illumination model for polymesh

	// rgb(244, 250, 252)
	//Phong glass; 
	glass.ambient = Colour(244.0/255.0, 250.0/255.0, 252.0/255.0, 255.0/255.0);
	glass.diffuse = Colour(244.0/255.0, 250.0/255.0, 252.0/255.0, 255.0/255.0);
	glass.specular = Colour(255.0/255.0, 255.0/255.0, 255.0/255.0, 255.0/255.0);
	glass.power = 40.0f;

    //Phong bp1; 
	// rgb(211, 141, 255)
	mat_pm.ambient = Colour(211.0/255.0, 141.0/255.0, 255.0/255.0, 255.0/255.0);
	mat_pm.diffuse = Colour(211.0/255.0, 141.0/255.0, 255.0/255.0, 255.0/255.0);
	mat_pm.specular = Colour(255.0/255.0, 255.0/255.0, 255.0/255.0, 255.0/255.0);
	mat_pm.power = 40.0f;

	// rgb(255, 252, 230)
	//Phong bp2;
	mat_wall2.ambient = Colour(255.0/255.0, 252.0/255.0, 230.0/255.0, 255.0/255.0);
	mat_wall2.diffuse = Colour(255.0/255.0, 252.0/255.0, 230.0/255.0, 255.0/255.0);
	mat_wall2.specular = Colour(255.0/255.0, 255.0/255.0, 255.0/255.0, 255.0/255.0);
	mat_wall2.power = 40.0f;

	// rgb(247, 198, 198)
	//Phong bp3;
	mat_wall3.ambient = Colour(247.0/255.0, 198.0/255.0, 198.0/255.0, 255.0/255.0);
	mat_wall3.diffuse = Colour(247.0/255.0, 198.0/255.0, 198.0/255.0, 255.0/255.0);
	mat_wall3.specular = Colour(255.0/255.0, 255.0/255.0, 255.0/255.0, 255.0/255.0);
	mat_wall3.power = 40.0f;

	// rgb(199, 247, 198)
	//Phong bp4;
	mat_wall4.ambient = Colour(199.0/255.0, 247.0/255.0, 198.0/255.0, 255.0/255.0);
	mat_wall4.diffuse = Colour(199.0/255.0, 247.0/255.0, 198.0/255.0, 255.0/255.0);
	mat_wall4.specular = Colour(255.0/255.0, 255.0/255.0, 255.0/255.0, 255.0/255.0);
	mat_wall4.power = 40.0f;

	// rgb(198, 242, 247)
	//Phong bp5;
	mat_wall5.ambient = Colour(198.0/255.0, 242.0/255.0, 247.0/255.0, 255.0/255.0);
	mat_wall5.diffuse = Colour(198.0/255.0, 242.0/255.0, 247.0/255.0, 255.0/255.0);
	mat_wall5.specular = Colour(255.0/255.0, 255.0/255.0, 255.0/255.0, 255.0/255.0);
	mat_wall5.power = 40.0f;

	// rgb(255, 176, 249)
	//Phong bp6;
	mat_wall6.ambient = Colour(255.0/255.0, 176.0/255.0, 249.0/255.0, 255.0/255.0);
	mat_wall6.diffuse = Colour(255.0/255.0, 176.0/255.0, 249.0/255.0, 255.0/255.0);
	mat_wall6.specular = Colour(255.0/255.0, 255.0/255.0, 255.0/255.0, 255.0/255.0);
	mat_wall6.power = 40.0f;

	// rgb(180, 169, 245)
	//Phong bp7;
	mat_wall7.ambient = Colour(180.0/255.0, 169.0/255.0, 245.0/255.0, 255.0/255.0);
	mat_wall7.diffuse = Colour(180.0/255.0, 169.0/255.0, 245.0/255.0, 255.0/255.0);
	mat_wall7.specular = Colour(255.0/255.0, 255.0/255.0, 255.0/255.0, 255.0/255.0);
	mat_wall7.power = 40.0f;

	float scaling = 1.0/M_PI;

    pm->material = &glass;
    pm->material->transparent = true;
    pm->material->reflective = true;
	pm->material->eta = 1.5; //glass refractive index
	pm->material->kr = 0.4;
	pm->material->BRDF_d = glass.diffuse;
	pm->material->BRDF_d.scale(scaling);

    Vertex v;
	v.x = 1.0f;
	v.y = 1.0f;
	v.z = 3.0f;
    Sphere *sphere = new Sphere(v, 1.0f);

    sphere->material = &mat_pm;
    sphere->material->transparent = false;
    sphere->material->reflective = false;
	sphere->material->eta = 1.5; //glass refractive index
	sphere->material->kr = 0.4;
	sphere->material->BRDF_d = mat_pm.diffuse;
	sphere->material->BRDF_d.scale(scaling); // mat_sphere.diffuse;

    
	//sphere->next = nullptr;

	/// put the scene in a box -- maybe only have 4/5 walls for light, see how it goes
	// rewrite this as a loop when you incorporate scene.cpp

	PolyMesh *fl = new PolyMesh((char *)"square.ply");
	PolyMesh *ce = new PolyMesh((char *)"ceiling.ply");
	PolyMesh *w1 = new PolyMesh((char *)"wall1.ply");
	PolyMesh *w2 = new PolyMesh((char *)"wall2.ply");
	PolyMesh *w3 = new PolyMesh((char *)"wall3.ply");
	PolyMesh *w4 = new PolyMesh((char *)"wall4.ply");

	
	fl->material = &mat_wall2;
    fl->material->transparent = false;
    fl->material->reflective = false;
	fl->material->BRDF_d = mat_wall2.diffuse;
	fl->material->BRDF_d.scale(scaling);
	fl->material->kr = 0.4;

	ce->material = &mat_wall3;
    ce->material->transparent = false;
    ce->material->reflective = false;
	ce->material->BRDF_d = mat_wall3.diffuse;
	ce->material->BRDF_d.scale(scaling);
	ce->material->kr = 0.4;

	w1->material = &mat_wall4;
    w1->material->transparent = false;
    w1->material->reflective = false;
	w1->material->BRDF_d = mat_wall4.diffuse;
	w1->material->BRDF_d.scale(scaling);
	w1->material->kr = 0.4;

	w2->material = &mat_wall5;
    w2->material->transparent = false;
    w2->material->reflective = false;
	w2->material->BRDF_d = mat_wall5.diffuse;
	w2->material->BRDF_d.scale(scaling);
	w2->material->kr = 0.4;
	
	w3->material = &mat_wall6;
    w3->material->transparent = false;
    w3->material->reflective = false;
	w3->material->BRDF_d = mat_wall6.diffuse;
	w3->material->BRDF_d.scale(scaling);
	w3->material->kr = 0.4;

	w4->material = &mat_wall7;
    w4->material->transparent = false;
    w4->material->reflective = false;
	w4->material->BRDF_d = mat_wall7.diffuse;
	w4->material->BRDF_d.scale(scaling);
	w4->material->kr = 0.4;

	//object_list = pm;
	object_list = pm;
	object_list->next = sphere;
	sphere->next = fl;
	fl->next = ce;
	ce->next = w1;
	w1->next = w2;
	w2->next = w3;
	w3->next = w4;
	w4->next = nullptr;

	// creates a light source
	Vertex v1 = Vertex(-1.0, -1.0, -1.0);
	Colour c = Colour(1.0, 1.0, 1.0, 1.0);
	Vector d = Vector(0.1, 0.1, 0.1);
	PointLight *pl = new PointLight(v1, c, d);
	pl->next = nullptr;

	// FIXME -- define object_list and light_list
	
	light_list = pl;
	light_list->next = nullptr;

}
*/


void Scene::teapot_box() {
	
	/// A NICE SCENE -- TEAPOT AND SPHERE INSIDE AN OPEN BOX

	// The following transform allows 4D homogeneous coordinates to be transformed. It moves the supplied teapot model to somewhere visible.
	Transform *transform = new Transform(1.0f, 0.0f, 0.0f,  0.0f,
			0.0f, 0.0f, 1.0f, -2.7f,
			0.0f, 1.0f, 0.0f, 5.0f,
			0.0f, 0.0f, 0.0f, 1.0f);

	//  Read in the teapot model.
	PolyMesh *pm = new PolyMesh((char *)"teapot_smaller.ply", transform);
	//objects->next = pm;

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
	pm->material->kt = 0.4;

    Vertex v;
	v.x = 1.0f;
	v.y = 1.0f;
	v.z = 5.0f;
    Sphere *sphere = new Sphere(v, 1.0f);

    sphere->material = &glass;
    sphere->material->transparent = true;
    sphere->material->reflective = true;
	sphere->material->eta = 1.5; //glass refractive index
	sphere->material->kr = 0.4;

    pm->next = sphere;
	//sphere->next = nullptr;

	/// put the scene in a box -- maybe only have 4/5 walls for light, see how it goes
	// rewrite this as a loop when you incorporate scene.cpp

	Transform *t = new Transform();

	PolyMesh *fl = new PolyMesh((char *)"square.ply", t);
	PolyMesh *ce = new PolyMesh((char *)"ceiling.ply", t);
	PolyMesh *w1 = new PolyMesh((char *)"wall1.ply", t);
	PolyMesh *w2 = new PolyMesh((char *)"wall2.ply", t);
	PolyMesh *w3 = new PolyMesh((char *)"wall3.ply", t);
	PolyMesh *w4 = new PolyMesh((char *)"wall4.ply", t);

	
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
    //DirectionalLight *dl = new DirectionalLight(Vector(0.01f, 0.01f, 1.0f),Colour(1.0f, 1.0f, 1.0f, 1.0f)); 
	//lights->next = dl;
	//dl->next = nullptr; 

	Vertex v1 = Vertex(0.2, 0.2, 0.2);
	Colour c = Colour(1.0f, 1.0f, 1.0f, 1.0f);
	Vector d = Vector(0.01f, 0.01f, 1.0f);


	PointLight *pl = new PointLight(v1, c, d);

	// FIXME -- define object_list and light_list
	object_list = pm;
	object_list->next = sphere;
	light_list = pl;
	light_list->next = nullptr;

	std::cout << "pm: " << pm << std::endl;
	std::cout << "sphere: " << sphere << std::endl;
	std::cout << "fl: " << fl << std::endl;
	std::cout << "ce: " << ce << std::endl;
	std::cout << "w1: " << w1 << std::endl;
	std::cout << "w2: " << w2 << std::endl;
	std::cout << "w3: " << w3 << std::endl;
	std::cout << "w4: " << w4 << std::endl;
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
	//cos_theta_t = sqrt(1.0f - sin2_t);

	if(cos_theta_t >= 0) {

		// should I use n insteaed of hit.normal?
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
	//this is the value needed for reflection in this surface at this point
	// it is used in the reflection section
	//hit.what->material->kr = kr;
	//hit.what->material->kt = kt;
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


Colour Scene::gather_diffuse(const Hit best_hit, const vector<Photon*> globalNeighbours) {

	Colour diffuse = Colour();
	float r=0; //get radius of the sphere from max distance between best_hit.position and rest of nodes
	float d=0, dA=0; 

	//find where the photons came from
	vector<Ray> initial_rays;

	for (Photon* p : globalNeighbours) {

		float dot = -1.0* (p->direction).dot(best_hit.normal);

		diffuse.r += p->BRDF_d.r * (p->intensity.r) * dot;
		diffuse.g += p->BRDF_d.g * (p->intensity.g) * dot;
		diffuse.b += p->BRDF_d.b * (p->intensity.b) * dot;

		//find max distance
		d = best_hit.position.distance(p->position);
		if (d > r) r = d;
	}

	if (r>0) {
		dA = M_PI * r * r; //area of disc associated with sphere
		float scaling = (1.0/dA);
		diffuse.scale(scaling);
	}

	if (diffuse.r != 0 || diffuse.g != 0 || diffuse.b != 0) {
		//std::cout << "diffuse: " << diffuse.r << ", " << diffuse.g << ", " << diffuse.b << std::endl;
	}

	return diffuse;
}



/// returns the depth and colour of the closest intersection between ray and the list of objects via classic Whitted-style ray tracing with global and caustic photon mapping 
void Scene::raytrace(Ray ray, Colour &colour, float &depth, int ref_limit, Hit &best_hit) {

	if (ref_limit < 0) return;

	// check for point light intersection, gives us L_e value 
	PointLight* pl;
	bool flag = false;

	point_light_intersection(ray, pl, depth, flag);

	if (flag) {
		// we've intersected a light source, so return the light colour
		colour = pl->intensity;
		return;
	}

	//Hit best_hit;
	Hit shadow_hit;
	shadow_hit.flag = false;
	object_intersection(ray, best_hit);

	// if we found a primitive then compute the colour we should see
	if(best_hit.flag) {

		//std::cout << "colour 1: " << colour.r << ", " << colour.g << ", " << colour.b << std::endl;
		// get colour of the material we've hit
		best_hit.what->material->compute_base_colour(colour);

		//std::cout << "colour 2: " << colour.r << ", " << colour.g << ", " << colour.b << std::endl;
		// get the depth of the object
		depth = best_hit.t;

		// get a handle on the first light -- check this
		PointLight *light = light_list;
		// check all the light sources for shadow rays
		// tis gives us direct lighting L_d
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
			//std::cout << "lit: " << lit << std::endl;

			if(ldir.dot(best_hit.normal)>0)	{
				lit=false;//light is facing wrong way.
			}
			if(lit) {
				// create a "shadow ray" from best_hit to the light source
				// if it intersects anything, then the point is in shadow
				
				Ray shadow_ray;

				// you want this to be the vector from the hit point to the point light
				// ie pl - best_hit.position

				//shadow_ray.direction.x = light->point.x - best_hit.position.x;
				//shadow_ray.direction.y = light->point.y - best_hit.position.y;
				//shadow_ray.direction.z = light->point.z - best_hit.position.z;

				shadow_ray.direction.x = -ldir.x;
				shadow_ray.direction.y = -ldir.y;
				shadow_ray.direction.z = -ldir.z;

				shadow_ray.position.x = best_hit.position.x + (0.0001f * shadow_ray.direction.x);
				shadow_ray.position.y = best_hit.position.y + (0.0001f * shadow_ray.direction.y);
				shadow_ray.position.z = best_hit.position.z + (0.0001f * shadow_ray.direction.z);

				// debug: print where the shadow ray goes, its t‐limit, and what it hits
				/* std::cout
				<< " ShadowRay from (" << best_hit.position.x << ","
										<< best_hit.position.y << ","
										<< best_hit.position.z << ") "
				<< "dir (" << shadow_ray.direction.x << ","
							<< shadow_ray.direction.y << ","
							<< shadow_ray.direction.z << ") "
				<< "\n"; */

				auto temp = best_hit.position;
				float dist = temp.distance(light->point);

				/* std::cerr
					<< " Shadow test: origin=("
					<< shadow_ray.position.x << ","
					<< shadow_ray.position.y << ","
					<< shadow_ray.position.z << ") dir=("
					<< shadow_ray.direction.x << ","
					<< shadow_ray.direction.y << ","
					<< shadow_ray.direction.z << ") "
					<< "dist-to-light=" << dist
					<< "  t_hit=" << shadow_hit.t
					<< "ldir: "
					<< ldir.x << ", "
					<< ldir.y << ", "
					<< ldir.z << ", "
					<< "\n"; */



				object_intersection(shadow_ray, shadow_hit);

				//std::cout << "shadow_hit.what: " << shadow_hit.what << std::endl;
				//std::cout << "shadow_hit.position: " << shadow_hit.position.x << ", " << shadow_hit.position.y << ", " << shadow_hit.position.z  << std::endl;
				//std::cout << "shadow_hit.t: " << shadow_hit.t << std::endl;

				

				//std::cout << "dist: " << dist << std::endl;



				//there's a shadow so no lighting, if  the hit is closer than the point light
				if(shadow_hit.flag && (shadow_hit.t < dist)) { //1000000000.0f) {
					lit = false;
					//printf("0");
				}
				else {
					//printf("1");
				}

			}
			//do colour
			if (lit) {

				//printf("we are lit!!");
				Colour intensity;
				Colour scaling;

				light->get_intensity(best_hit.position, scaling);
				best_hit.what->material->compute_light_colour(viewer, best_hit.normal, ldir, intensity);
				
				//scale by scaling to get intensity
				intensity.scale(scaling);

				//add intensity to the colour
				colour.add(intensity);
				//std::cout << "colour 3: " << colour.r << ", " << intensti.g << ", " << intensity.b << std::endl;
				}
			// move on to the next light
			light = light->next;
		}
		
		//std::cout << "colour 3: " << colour.r << ", " << colour.g << ", " << colour.b << std::endl;
		// generate secondary ray (reflection and/or refraction)
		// ref_limit stops infinite recursion
		// this is some L_s, specular
		if(best_hit.what->material->transparent)
		{
			//std::cout << "in transparent" << std::endl;
			Ray refracted;
			bool total_internal_reflection = false;

			refract_ray(ray, best_hit, refracted, total_internal_reflection);

			// I dont think you need to save these values here, but just in case something goes wrong try uncommenting this
			//float kt = best_hit.what->material->kt;
			//float kr = best_hit.what->material->kr;
			
			if (!total_internal_reflection) {
				Colour trans_colour;
				float trans_depth;
				Hit h = Hit();

				ref_limit -= 1;

				//raytrace the refracted ray
				raytrace(refracted, trans_colour, trans_depth, ref_limit, h);

				//scale trans_colour by kt
				trans_colour.scale(best_hit.what->material->kt); 

				//add the refracted colour to the colour
				colour.add(trans_colour);
			}

		}

		//std::cout << "colour 4: " << colour.r << ", " << colour.g << ", " << colour.b << std::endl;

		// compute reflection ray if material supports it.
		// this is more L_s, specular
		if(best_hit.what->material->reflective)	{
			Ray reflected;
			Hit rhit;

			reflect_ray(ray, best_hit, reflected);

			//tests if there is an object to be reflected in the surface
			object_intersection(reflected, rhit); 

			ref_limit -= 1;

			if(ref_limit<0)	return;

			if(rhit.flag) {

				float ref_depth;
				Colour ref_colour;
				Hit h = Hit();

				//get colour of reflection 
				raytrace(reflected, ref_colour, ref_depth, ref_limit, h); 

				//kr = best_hit.what->material->kr;

				//scale by kr
				ref_colour.scale(best_hit.what->material->kr);

				//add reflective colour to the colour
				colour.add(ref_colour);
			}
		}
		
		//std::cout << "colour 5: " << colour.r << ", " << colour.g << ", " << colour.b << std::endl;
		// add photon map checks here

		// this is L_c, caustic

		if (causticTree) {
			vector<Photon*> causticNeighbours;
			causticTree->kNearest(best_hit.position, 50, causticNeighbours);

			//std::cout << "causticNeighbours.size: " << causticNeighbours.size() << std::endl;
			//printf("out of knearest");
			Colour caustic = Colour();

			float r = 0.0f; //get radius of the sphere from max distance between best_hit.position and rest of nodes
			float d=0.0f; //, max_d = 0;

			//each photon in sphere
			for (Photon* p : causticNeighbours) {

				//find colour
				float dot = -1.0f * p->direction.dot(best_hit.normal);
				
				// need temp because we're in a for loop and don't want to overwrite caustic each iteration
				Colour temp = p->BRDF_d;
				temp.scale(p->intensity);
				temp.scale(dot);

				caustic.add(temp);

				//find max distance 
				d = best_hit.position.distance(p->position);
				//std::cout<< "d: " << d << std::endl;
				if (d > r) r = d;
			}
			//std::cout<< "r: " << r << std::endl;
			float dA = M_PI * r * r; //area of disc associated with sphere of radius r
			//normalise by the disc area (check it's positive maybe? FIXME)
			if (r>0) { 
				float scaling = 1.0/dA;
				caustic.scale(scaling); 
			}
			float photonBoost = 100.0f;             // try 10, 100, 1000…
			Colour boosted = caustic;
			boosted.scale(photonBoost);
			colour.add(boosted);

			//colour.add(caustic);
			if (caustic.r != 0 || caustic.g != 0 || caustic.b != 0) {
				//std::cout << "caustic: " << caustic.r << ", " << caustic.g << ", " << caustic.b << std::endl;
			}
		}

		//std::cout << "colour 6: " << colour.r << ", " << colour.g << ", " << colour.b << std::endl;
		// this is L_d, diffuse

		if (globalTree) {

			vector<Photon*> globalNeighbours;	
			globalTree->kNearest(best_hit.position, 5, globalNeighbours);

			Colour diffuse = gather_diffuse(best_hit, globalNeighbours);
			float photonBoost = 100.0f;             // try 10, 100, 1000…
			Colour boosted = diffuse;
			boosted.scale(photonBoost);
			colour.add(boosted);

			//colour.add(diffuse);

			// get rays -- TODO look into emplace_back() instead
			vector<Ray> current_rays;
			for (Photon* p : globalNeighbours) {
				Ray new_r = Ray(p->position, p->direction);
				current_rays.push_back(new_r);
			}
			// ----------------------------------------------

			//store the intensity at this point, multiply by kr of the object

			//use importance sampling to decide which to follow and then resample at the next hit point (reflection)

			//bounce limit for photons
			int bounce_limit = 3;
			

			while (bounce_limit > 0 && !current_rays.empty()) {

				bounce_limit--;
				Colour diffuse_reflection = Colour();
				vector<Ray> nxt_rays;

				//ideally for imp_rays
				for (Ray &ray : current_rays) {
				
					Hit nxt_hit;
					object_intersection(ray, nxt_hit); //get next hit point of each ray

					if(!nxt_hit.flag) continue;

					//do k nearest neighbours sampling again, but scale each reflection by kr (BRDF?) and add them all together
					vector<Photon*> globalNeighbours;
					globalTree->kNearest(nxt_hit.position, 5, globalNeighbours);

					Colour d = gather_diffuse(nxt_hit, globalNeighbours);
					d.scale(nxt_hit.what->material->kr);
					diffuse_reflection.add(d);

					// get next rays
					for (Photon* p : globalNeighbours) {
						Ray new_r = Ray(p->position, p->direction);
						nxt_rays.push_back(new_r);
					}
				}
				float photonBoost = 100.0f;             // try 10, 100, 1000…
				Colour boosted = diffuse_reflection;
				boosted.scale(photonBoost);
				colour.add(boosted);

				//colour.add(diffuse_reflection);
				current_rays.swap(nxt_rays);
				if (diffuse_reflection.r != 0 || diffuse_reflection.g != 0 ||diffuse_reflection.b != 0 ) {
					//std::cout << "diffuse_reflection: " << diffuse_reflection.r << ", " << diffuse_reflection.g << ", " << diffuse_reflection.b << std::endl;
				}
			}
			if (diffuse.r != 0 || diffuse.g != 0 || diffuse.b != 0) {
				//std::cout << "diffuse: " << diffuse.r << ", " << diffuse.g << ", " << diffuse.b << std::endl;
			}
		}
	
		//std::cout << "colour 7: " << colour.r << ", " << colour.g << ", " << colour.b << std::endl;
	} else {
		// colour = background_colour
		// I like rgb(163, 249, 255)
		colour.r = 163.0/255.0;
		colour.g = 249.0/255.0;
		colour.b = 255.0/255.0;
		depth = 7.0f;
	}
	
	//std::cout << "colour 8: " << colour.r << ", " << colour.g << ", " << colour.b << std::endl;
	return;
}


//doesn't emit towards specular surfaces, it emits in random directions then checks if the intersection is specular
//traces a single photon through the scene, saves them in causticPhotons and globalPhotons
void Scene::photon_trace(Photon *photon, int ref_limit) {

	//printf("1");

	bool saw_specular = false;
	Ray photon_ray;
	photon->ray(photon_ray);

	//printf("2");
	// loop starts here
	// do you use ref_limit or photon->absorbed and russian roulette?
	// what is ref_limit used for here?
	while (ref_limit > 0) {
		
		/*if (ref_limit < 5) {
			std::cout << "ref_limit: " << ref_limit << std::endl;
			std::cout << "first saw_specular: " << saw_specular << std::endl;
			
		}*/

		
		ref_limit--;
		//printf("3a");
		Hit best_hit;
		object_intersection(photon_ray, best_hit); //where does this photon hit first? maybe add in the shadow photons here


		//std::cout << "photon ray: " << photon_ray << std::endl;
		//printf("3b");
		if (!best_hit.flag) break;

		//printf("4");
		photon->position = best_hit.position;

		//std::cout << "hit object: " << best_hit.what << std::endl;

		bool is_specular = best_hit.what->material->reflective || best_hit.what->material->transparent;

		if (!is_specular) {
			//printf("5");
			// add to global and check for caustic
			photon->position = best_hit.position;
			photon->BRDF_d = best_hit.what->material->BRDF_d;
			//photon->BRDF_s = best_hit.what->material->BRDF_s;
			globalPhotons.push_back(photon);

			//std::cout << "globalPhotons.size: " << globalPhotons.size() << std::endl;

			if (saw_specular) {
				//printf("6");
				causticPhotons.push_back(photon);
				// reset the caustic photon
				saw_specular = false;
			}
			//photon->g_russian_roulette(best_hit);
			// FIXME bounce function to get the next ray
			// actually, start with just one diffuse hit, inc for global photons
			break;
		} else {
			//printf("7");
			saw_specular = true;
			photon->c_russian_roulette(best_hit);
			Ray new_ray;

			//std::cout << "transmitted? " << photon->transmitted << std::endl;
			//std::cout << "reflected? " << photon->reflected << std::endl;

			if (photon->reflected) {
				//printf("8");
				// change the intensity -- here or in c_russian_roulette?
				// don't want to kill a specular photon in russian roulette, 
				// it only terminates if it hits a diffuse surface, or the bounce limit runs out
				reflect_ray(photon_ray, best_hit, new_ray);
				photon_ray.position = new_ray.position;
				photon_ray.direction = new_ray.direction;
			}
			if (photon->transmitted) {
				bool tir = false;
				//std::cout << "photon_ray: " << photon_ray << std::endl;
                refract_ray(photon_ray, best_hit, new_ray, tir);
				photon_ray.position = new_ray.position;
				photon_ray.direction = new_ray.direction;
				//std::cout << "photon_ray: " << photon_ray << std::endl;
				if (tir) {
                    // if total internal reflection, treat as mirror
                    reflect_ray(photon_ray, best_hit, new_ray);
					photon_ray.position = new_ray.position;
					photon_ray.direction = new_ray.direction;
					//std::cout << "photon_ray: " << photon_ray << std::endl;
                }
			}
		}
		//printf("10");
		//std::cout << "saw specular: " << saw_specular << std::endl;
	}
	//printf("11");
}


/// creates caustic and global photon maps, saves them in causticTree and globalTree
void Scene::create_photon_maps() {
	//caustic_photon_map(5);
	//global_photon_map(5);

	causticPhotons.clear();

	globalPhotons.clear();

	PointLight *light = light_list;
	int no_of_photons = 1000000;

	for (int n = 0; n < no_of_photons; n++) {
		while (light != (PointLight*) 0) {
			//trace a photon through the scene and store all the hits in a kd tree - how to go through specular objects only?
			Photon *photon = new Photon(*light, no_of_photons);
			if (n < 100000) {
			// First 20k: aim at sphere center
			Vector temp = Vector(-0.2, -0.2, 1.8);
			temp.normalise();
			photon->direction = temp;
			} else if (n < 200000) {
				Vector temp = Vector(0.8, 0.8, 2.8);
				temp.normalise();
			}
			photon_trace(photon, 5);
			light = light->next;
		}
		light = light_list;
		if (n % 1000 == 0) { 
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

