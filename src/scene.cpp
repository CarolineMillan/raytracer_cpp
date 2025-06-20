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


Scene::Scene()
{
	object_list = 0;
	light_list = 0;
}


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

    pm->material = &mat_pm;
    pm->material->transparent = false;
    pm->material->reflective = false;
	pm->material->eta = 1.5; //glass refractive index
	pm->material->kr = 0.4;

    Vertex v;
	v.x = 0.0f;
	v.y = 0.0f;
	v.z = 0.5f;
    Sphere *sphere = new Sphere(v, 0.25f);

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
	w3->next = nullptr;
	//w4->next = nullptr;

	//creates a light source
    DirectionalLight *dl = new DirectionalLight(Vector(0.01f, 0.01f, 1.0f),Colour(1.0f, 1.0f, 1.0f, 1.0f)); 
	//lights->next = dl;
	dl->next = nullptr; 

	// FIXME -- define object_list and light_list
	object_list = pm;
	light_list = dl;

}


void Scene::test() {
	
	/// A SIMPLE SCENE -- SPHERE, BACK WALL, ONE SIDE WALL
	// rgb(211, 141, 255)
	// rgb(244, 250, 252)
	//Phong bp1; 
	mat_sphere.ambient = Colour(244.0/255.0, 250.0/255.0, 252.0/255.0, 255.0/255.0);
	mat_sphere.diffuse = Colour(244.0/255.0, 250.0/255.0, 252.0/255.0, 255.0/255.0);
	mat_sphere.specular = Colour(255.0/255.0, 255.0/255.0, 255.0/255.0, 255.0/255.0);
	mat_sphere.power = 40.0f;

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
	v.x = 1.0f;
	v.y = 1.0f;
	v.z = 5.0f;
    Sphere *sphere = new Sphere(v, 1.0f);

    sphere->material = &mat_sphere;
    sphere->material->transparent = true;
    sphere->material->reflective = true;
	sphere->material->eta = 1.5; //glass refractive index
	sphere->material->kr = 0.4;

	PolyMesh *w1 = new PolyMesh((char *)"wall1.ply");

	w1->material = &mat_wall1;
    w1->material->transparent = false;
    w1->material->reflective = false;

	PolyMesh *w3 = new PolyMesh((char *)"wall3.ply");
	
	w3->material = &mat_wall3;
    w3->material->transparent = false;
    w3->material->reflective = false;

	object_list = sphere;
	object_list->next = w1;
	w1->next = w3;
	w3->next = nullptr;

	//creates a light source
    DirectionalLight *dl = new DirectionalLight(Vector(0.01f, 0.01f, 1.0f),Colour(1.0f, 1.0f, 1.0f, 1.0f)); 
	light_list = dl;
	dl->next = nullptr; 
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


/// returns the depth and colour of the closest intersection between ray and the list of objects via classic Whitted-style ray tracing
void Scene::raytrace(Ray ray, Colour &colour, float &depth, int ref_limit) {

	if (ref_limit < 0) return;

	Hit best_hit;
	Hit shadow_hit;
	object_intersection(ray, best_hit);

	// if we found a primitive then compute the colour we should see
	if(best_hit.flag) {
		// get colour of the material we've hit
		best_hit.what->material->compute_base_colour(colour);
		// get the depth of the object
		depth = best_hit.t;

		// get a handle on the first light -- check this
		Light *light = light_list;
		// check all the light sources for shadow rays
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
			if(ldir.dot(best_hit.normal)>0)	lit=false;//light is facing wrong way.
			if(lit) {
				// create a "shadow ray" from best_hit to the light source
				// if it intersects anything, then the point is in shadow
				Ray shadow_ray;

				shadow_ray.direction.x = -ldir.x;
				shadow_ray.direction.y = -ldir.y;
				shadow_ray.direction.z = -ldir.z;
				shadow_ray.position.x = best_hit.position.x + (0.0001f * shadow_ray.direction.x);
				shadow_ray.position.y = best_hit.position.y + (0.0001f * shadow_ray.direction.y);
				shadow_ray.position.z = best_hit.position.z + (0.0001f * shadow_ray.direction.z);

				object_intersection(shadow_ray, shadow_hit);

				//there's a shadow so no lighting, if realistically close
				if(shadow_hit.flag && shadow_hit.t < 1000000000.0f) lit = false;
			}
			//do colour
			if (lit) {
					Colour intensity;
					Colour scaling;

					light->get_intensity(best_hit.position, scaling);
					best_hit.what->material->compute_light_colour(viewer, best_hit.normal, ldir, intensity);
					
					//scale by scaling to get intensity
					intensity.scale(scaling);

					//add intensity to the colour
					colour.add(intensity);
				}
			// move on to the next light
			light = light->next;
		}
		// generate secondary ray (reflection and/or refraction)
		// ref_limit stops infinite recursion
		if(best_hit.what->material->transparent)
		{
			Ray T;
			Colour trans_colour, refl_colour;
			float trans_depth, eta, cos_theta_t, cos_theta_i, rpar, rper, kr, kt;

			eta = best_hit.what->material->eta;

			Vector n;
			n =  best_hit.normal;
			n.normalise();
			Vector d;
			d = ray.direction;
			d.normalise();
			cos_theta_i = fabs(n.dot(-1*d));

			if (1.0 - ( 1.0 / (pow(eta,2) ) >=0.0 )) {
				cos_theta_t = sqrt(1.0 - ( 1.0 / (pow(eta,2.0) ) ) * ( 1.0 - pow(cos_theta_i,2.0) ) );
			}

			float sin2_t = (1.0f / (eta * eta)) * (1.0f - cos_theta_i * cos_theta_i);
			
			// why do we test sin2_t and not if(cos_theta_t >= 0) ?
			
			if (sin2_t <= 1.0f) {
				cos_theta_t = sqrt(1.0 - ( 1.0 / (pow(eta,2.0) ) ) * ( 1.0 - pow(cos_theta_i,2.0) ) );
				//cos_theta_t = sqrt(1.0f - sin2_t);

				T.direction = (1/eta)*ray.direction - (cos_theta_t - (1/eta)*cos_theta_i)*best_hit.normal;

				T.position.x = best_hit.position.x + 0.001*T.direction.x;
				T.position.y = best_hit.position.y + 0.001*T.direction.y;
				T.position.z = best_hit.position.z + 0.001*T.direction.z;

				ref_limit -= 1;

				//raytrace the refracted ray
				raytrace(T, trans_colour, trans_depth, ref_limit);

				//fresnel equations
				rpar = (eta*cos_theta_i - cos_theta_t)/(eta*cos_theta_i + cos_theta_t);
				rper = (cos_theta_i - eta*cos_theta_t)/(cos_theta_i + eta*cos_theta_t);
				kr = 0.5*( pow(rpar,2) + pow(rper,2) ); 
				kt = 1-kr;

				//this is the value needed for reflection in this surface at this point
				// it is used in the reflection section
				best_hit.what->material->kr = kr;
				best_hit.what->material->kt = kt;

				if (kr < 0 || kr > 1) std::cout << "kr: " << kr << std::endl;
				if (kt < 0 || kt > 1) std::cout << "kt: " << kt << std::endl;

				//scale trans_colour by kt
				trans_colour.r *= kt; //use scale function?
				trans_colour.g *= kt;
				trans_colour.b *= kt;

				//add the refracted colour to the colour
				colour.add(trans_colour);
			}

		}
		
		// compute reflection ray if material supports it.
		if(best_hit.what->material->reflective)	{
			Ray R;
			float kr, ref_depth;
			Hit rhit;
			Colour ref_colour;

			//calculate reflected eye ray
			R.direction = ray.direction - 2*((best_hit.normal).dot(ray.direction))*best_hit.normal;

			R.position.x = best_hit.position.x + 0.001*R.direction.x;
			R.position.y = best_hit.position.y + 0.001*R.direction.y;
			R.position.z = best_hit.position.z + 0.001*R.direction.z;

			//tests if there is an object to be reflected in the surface
			object_intersection(R, rhit); 

			ref_limit -= 1;

			if(ref_limit<0)	return;

			if(rhit.flag) {
				//get colour of reflection 

				raytrace(R, ref_colour, ref_depth, ref_limit); 

				kr = best_hit.what->material->kr;

				//scale by kr
				ref_colour.r *= kr; 
				ref_colour.g *= kr;
				ref_colour.b *= kr;

				//add reflective colour to the colour
				colour.add(ref_colour);
			}
		}
	} else {
		// colour = background_colour
		// I like rgb(163, 249, 255)
		colour.r = 163.0/255.0;
		colour.g = 249.0/255.0;
		colour.b = 255.0/255.0;
		depth = 7.0f;
	}
	return;
}



/* PHOTON MAPPING STUFF

void g_trace(Photon *photon, Object *objects, PointLight *pl, int ref_limit, Kdtree::KdNodeVector g_nodes) //traces a single photon through the scene, saves them as a vector of nodes
		{
	Ray photon_ray;
	photon->ray(photon_ray);

	Hit best_hit;

	object_test(photon_ray, objects, best_hit); //where does this photon hit first? maybe add in the shadow photons here

	photon->position = best_hit.position;

	Kdtree::KdNode p;

	std::vector<double> pos;

	pos[0] = photon->position.x;
	pos[1] = photon->position.y;
	pos[2] = photon->position.z;

	p.point = pos; //what point should i give it?
	p.data = photon;

	g_nodes.push_back(p);

	int new_w = 0;

	photon->g_russian_roulette(best_hit, photon->w, new_w);

	photon->w = new_w; //gets the new weight for russian_roulette

	while (!photon->absorbed) //russian roulette threshhold/weighting should mean the photon is eventually absorbed
	{
		if (photon->reflected) //calculate BDRF for reflection coefficient, use rendering equation later
		{

			Ray R;
			Hit rhit;

			//use BDRF to determine new intensity
			photon->BRDF_s = best_hit.what->material->BRDF_s;
			photon->BRDF_d = best_hit.what->material->BRDF_d;

			//add this photon to kd tree
			Kdtree::KdNode p;
			std::vector<double> pos;

			pos[0] = photon->position.x;
			pos[1] = photon->position.y;
			pos[2] = photon->position.z;

			p.point = pos;
			p.data = photon;

			g_nodes.push_back(p);
		}

		photon->g_russian_roulette(best_hit, photon->w, new_w); //need to do russian roulette again every iteration of the while loop

		photon->w = new_w;
	}

	//store last absorbed photon of the while loop in the tree - if statement?
	pos[0] = photon->position.x;
	pos[1] = photon->position.y;
	pos[2] = photon->position.z;

	p.point = pos;
	;
	p.data = photon;

	g_nodes.push_back(p);

	//continue path and store shadow photons
	while (objects != 0) {
		Hit obj_hit;
		obj_hit.flag = false;

		//does the photon intersect with this object?
		objects->intersection(photon_ray, obj_hit);

		if (obj_hit.flag) {
			if (obj_hit.position.x == best_hit.position.x
					&& obj_hit.position.y == best_hit.position.y
					&& obj_hit.position.y == best_hit.position.y) {
				photon->position = obj_hit.position;
				photon->shadow = true;

				Kdtree::KdNode p;
				std::vector<double> pos;

				pos[0] = photon->position.x;
				pos[1] = photon->position.y;
				pos[2] = photon->position.z;

				p.point = pos;
				p.data = photon;

				g_nodes.push_back(p);
			}
		}

	}
}


//doesn't emit towards specular surfaces, it emits in random directions then checks if the intersection is specular
void c_trace(Photon *photon, Object *objects, int ref_limit, Kdtree::KdNodeVector c_nodes) //traces a single photon through the scene, saves them as a vector of nodes
		{
	Ray photon_ray;
	photon->ray(photon_ray);

	Hit best_hit;

	object_test(photon_ray, objects, best_hit); //where does this photon hit first? maybe add in the shadow photons here

	if (best_hit.what->material->transparent) //only does it with photons that directly hit a transparent surface, what if they are reflected onto a transparent surface?
	{
		photon->position = best_hit.position;

		Kdtree::KdNode p;

		std::vector<double> pos;

		pos[0] = photon->position.x;
		pos[1] = photon->position.y;
		pos[2] = photon->position.z;

		p.point = pos;
		; //what point should i give it?
		p.data = photon;

		c_nodes.push_back(p);

		//continue path and store shadow photons - not in caustic

		int new_w;

		photon->c_russian_roulette(best_hit, photon->w, new_w);

		photon->w = new_w;

		while (!photon->absorbed) {
			if (photon->transmitted) {
				//direction is a function type of transmission - new direction is
				//if its transmitted we need a new direction, use raytracing transparency stuff
				Kdtree::KdNode p;

				std::vector<double> pos;

				pos[0] = photon->position.x;
				pos[1] = photon->position.y;
				pos[2] = photon->position.z;

				p.point = pos;
				p.data = photon;
				c_nodes.push_back(p);
				photon->c_russian_roulette(best_hit, photon->w, new_w);
				photon->w = new_w;
			}
		}

		//stores the last absorbed photon of the while loop

		pos[0] = photon->position.x;
		pos[1] = photon->position.y;
		pos[2] = photon->position.z;

		p.point = pos;
		p.data = photon;

		c_nodes.push_back(p);
	}
}


void g_photon_map(Object *objects, PointLight *lights, int ref_limit, Kdtree::KdNodeVector g_nodes) //creates a global photon map
		{

	Kdtree::KdNodeVector g_nodes1;

	for (int n = 0; n == 100000; n++) {
		while (lights != (Light*) 0) {
			Photon *photon = new Photon(*lights);

			//trace a photon through the scene and store all the hits in a kd tree
			g_trace(photon, objects, lights, 6, g_nodes1);

			for (int i = 0; i == g_nodes1.size(); i++) //add all these photons to the tree
					{
				Kdtree::KdNode kdnode;

				kdnode.data = g_nodes1[i].data;
				kdnode.point = g_nodes1[i].point;

				g_nodes.push_back(kdnode);
			}

			lights = lights->next; //move on to next light
		}
	}
		}


void c_photon_map(Object *objects, PointLight *lights, int ref_limit, Kdtree::KdNodeVector c_nodes) //creates a caustic photon map, only deals with transparency, don't need depth?
		{
	Kdtree::KdNodeVector c_nodes1;

	for (int n = 0; n == 100000; n++) //1 million photons?
			{
		while (lights != (Light*) 0) {
			Photon *photon = new Photon(*lights);

			//trace a photon through the scene and store all the hits in a kd tree - how to go through specular objects only?
			c_trace(photon, objects, 6, c_nodes1);

			for (int i = 0; i == c_nodes1.size(); i++) {
				c_nodes.push_back(c_nodes1[i]); //adds all photon hits to a tree
			}
			lights = lights->next;
		}
	}
}


/// returns the depth and colour of the closest intersection between ray and the list of objects via classic Whitted-style ray tracing with global and caustic photon mapping 
void render(Ray ray, Object *objects, Light *lights, Colour &colour, float &depth, int ref_limit) {
				//Raytrace to get closest object

			Colour L_l, L_s, L_c, L_d, L_0, L_e;
			Hit best_hit;

			Object *obj = objects;
			Light *light = lights;

			object_intersection(ray, obj, best_hit);

			Vector L_i, V;

			L_e.r = 0;
			L_e.g = 0;
			L_e.b = 0;

			//is the hit point a light source?
			while (light != (Light*) 0) {

				if (light->point.x == best_hit.position.x
						&& light->point.y == best_hit.position.y
						&& light->point.z == best_hit.position.z) //if the hit position is on a light, then it emits radiance
								{
					L_e.r += light->intensity.r; //addition in case two lights are in the same position (unlikely)
					L_e.g += light->intensity.g;
					L_e.b += light->intensity.b;
				}
				light = light->next;
			}

			L_0 = L_e;

			//DIRECT - L_l
			raytrace(ray, obj, light, L_l, depth, best_hit);

			L_l.r *= (best_hit.what->material->BRDF_s.r + best_hit.what->material->BRDF_d.r) * (-1* (pl->direction).dot(best_hit.normal));
			L_l.g *= (best_hit.what->material->BRDF_s.g + best_hit.what->material->BRDF_d.g) * (-1* (pl->direction).dot(best_hit.normal));
			L_l.b *= (best_hit.what->material->BRDF_s.b + best_hit.what->material->BRDF_d.b) * (-1* (pl->direction).dot(best_hit.normal));

			//SPECULAR - L_s
			if (best_hit.what->material->reflective) //material supports reflection
			{
				Ray R;
				float kr, ref_depth;
				Hit rhit;

				//calculate reflected eye ray
				R.direction = ray.direction
						- 2 * ((best_hit.normal).dot(ray.direction))
								* best_hit.normal;

				R.position.x = best_hit.position.x + 0.001 * R.direction.x;
				R.position.y = best_hit.position.y + 0.001 * R.direction.y;
				R.position.z = best_hit.position.z + 0.001 * R.direction.z;

				object_intersection(R, obj, rhit); //tests if there is an object to be reflected in the surface

				if (rhit.flag)
				{
					Hit ref_hit;
					raytrace(R, obj, pl, L_s, ref_depth, ref_hit);

					L_s.r *= best_hit.what->material->BRDF_s.r * (-1* (pl->direction).dot(best_hit.normal));
					L_s.g *= best_hit.what->material->BRDF_s.g * (-1* (pl->direction).dot(best_hit.normal));
					L_s.b *= best_hit.what->material->BRDF_s.b * (-1* (pl->direction).dot(best_hit.normal));
				}
			} else {
				L_s.r = 0;
				L_s.g = 0;
				L_s.b = 0;
			}

			//CAUSTICS - L_c

			//create a photon map
			Kdtree::KdNodeVector nodes;
			c_photon_map(sphere, pl, 6, nodes);

			//create kd tree
			Kdtree::KdTree c_tree(&nodes);

			Kdtree::KdNodeVector c_nodes_knn;
			std::vector<double> v;

			v[0] = best_hit.position.x;
			v[1] = best_hit.position.y;
			v[2] = best_hit.position.z;

			const std::vector<double> v0 = v;

			c_tree.k_nearest_neighbors(v0, 50, c_nodes_knn);

			float r; //get radius of the sphere from max distance between best_hit.position and rest of nodes
			float a, b, c, d, max_d = 0;

			for (int i = 0; i == c_nodes_knn.size(); i++) //finds max distance
					{
				a = best_hit.position.x - c_nodes_knn[i].point[0];
				b = best_hit.position.y - c_nodes_knn[i].point[1];
				c = best_hit.position.z - c_nodes_knn[i].point[2];
				d = sqrt(pow(a, 2) + pow(b, 2) + pow(c, 2));

				if (d > max_d)
				{
					max_d = d;
				}
			}

			r = max_d;

			float dA = M_PI * pow(r, 2); //area of disc associated with sphere of radius r

			for (int i = 0; i == 49; i++) //each photon in sphere
					{
				Photon *p = new Photon();
				p = c_nodes_knn[i].data;

				//divide flux that photon represents by area of s and multiply by photon's BDRF to get radiance here
				L_c.r += p->BRDF_d.r * (p->intensity.r) * (-1* (p->direction).dot(best_hit.normal));
				L_c.g += p->BRDF_d.g * (p->intensity.g) * (-1* (p->direction).dot(best_hit.normal));
				L_c.b += p->BRDF_d.b * (p->intensity.b) * (-1* (p->direction).dot(best_hit.normal));
			}
			L_c.r *= 1/dA;
			L_c.g *= 1/dA;
			L_c.b *= 1/dA;

			//SOFT INDIRECT - L_d

			//create a photon map
			Kdtree::KdNodeVector nodes1;
			g_photon_map(sphere, pl, 6, nodes1);

			//create kd tree
			Kdtree::KdTree g_tree(&nodes1);

			Kdtree::KdNodeVector g_nodes_knn;

			g_tree.k_nearest_neighbors(v0, 50, g_nodes_knn);

			for (int i = 0; i == g_nodes_knn.size(); i++) //finds max distance
					{
				a = best_hit.position.x - g_nodes_knn[i].point[0];
				b = best_hit.position.y - g_nodes_knn[i].point[1];
				c = best_hit.position.z - g_nodes_knn[i].point[2];
				d = sqrt(pow(a, 2) + pow(b, 2) + pow(c, 2));

				if (d > max_d) {
					max_d = d;
				}
			}

			r = max_d;

			dA = M_PI * pow(r, 2); //area of disc associated with sphere

			//find where the photons came from
			vector<Ray> nxt_rays;

			for (int i = 0; i == 49; i++) {
				Photon *p;

				p = g_nodes_knn[i].data;

				L_d.r += p->BRDF_d.r * (p->intensity.r) * (-1* (p->direction).dot(best_hit.normal));
				L_d.g += p->BRDF_d.g * (p->intensity.g) * (-1* (p->direction).dot(best_hit.normal));
				L_d.b += p->BRDF_d.b * (p->intensity.b) * (-1* (p->direction).dot(best_hit.normal));

				nxt_rays[i].position.x = g_nodes_knn[i].point[0];
				nxt_rays[i].position.y = g_nodes_knn[i].point[1];
				nxt_rays[i].position.z = g_nodes_knn[i].point[2];
				nxt_rays[i].direction = p->direction; //direction the photon came from
			}

			L_d.r *= 1/dA;
			L_d.g *= 1/dA;
			L_d.b *= 1/dA;

			vector<float> kr;

			for (int k = 0; k == nxt_rays.size(); k++) {
				kr[k] = best_hit.what->material->kr; //BRDF?
			}
			//store the intensity at this point, multiply by kr of the object

			//use importance sampling to decide which to follow and then resample at the next hit point (reflection)

			//bounce limit for photons
			int ref_limit = 3;

			Colour col, L_g_ref; //global reflections colour, will then be added onto L_d.
			float g_depth;
			Hit nxt_hit;

			while (ref_limit > 0) {
				for (int k = 0; k == nxt_rays.size(); k++) //ideally for imp_rays
				{
					object_intersection(nxt_rays[k], obj, nxt_hit); //get next hit point of each ray

					//do k nearest neighbours sampling again, but scale each reflection by kr (BRDF?) and add them all together
					Kdtree::KdNodeVector g_nodes_knn_1;

					std::vector<double> v1;

					v1[0] = nxt_hit.position.x;
					v1[1] = nxt_hit.position.y;
					v1[2] = nxt_hit.position.z;

					const std::vector<double> v2;

					g_tree.k_nearest_neighbors(v2, 50, g_nodes_knn_1);

					float r; //get radius of the sphere from max distance between best_hit.position and rest of nodes
					float a, b, c, d, max_d = 0;

					for (int i = 0; i == g_nodes_knn_1.size(); i++) //finds max distance
					{

						a = best_hit.position.x - g_nodes_knn_1[i].point[0];
						b = best_hit.position.y - g_nodes_knn_1[i].point[1];
						c = best_hit.position.z - g_nodes_knn_1[i].point[2];
						d = sqrt(pow(a, 2) + pow(b, 2) + pow(c, 2));

						if (d > max_d)
						{
							max_d = d;
						}
					}

					r = max_d;

					float dA;
					dA = M_PI * pow(r, 2); //area of disc associated with sphere

					//now work out average colour
					for (int i = 0; i == 49; i++)
					{
						Photon *p;

						p = g_nodes_knn_1[i].data;

						col.r += kr[k] * p->BRDF * (p->intensity.r) * (-1* (p->direction).dot(best_hit.normal));
						col.g += kr[k] * p->BRDF * (p->intensity.g) * (-1* (p->direction).dot(best_hit.normal));
						col.b += kr[k] * p->BRDF * (p->intensity.b) * (-1* (p->direction).dot(best_hit.normal));

						nxt_rays[i].position.x = g_nodes_knn_1[i].point[0];
						nxt_rays[i].position.y = g_nodes_knn_1[i].point[1];
						nxt_rays[i].position.z = g_nodes_knn_1[i].point[2];
						nxt_rays[i].direction = p->direction; //direction the photon came from
					}

					col.r *= 1/dA;
					col.g *= 1/dA;
					col.b *= 1/dA;

					//get next kr value
					kr[k] = nxt_hit.what->material->kr;

					ref_limit -= 1;
				}
				L_d.r += col.r;
				L_d.g += col.g;
				L_d.b += col.b;
			}
			//TOTAL
			colour.r = L_e.r + L_l.r + L_s.r + L_c.r + L_d.r;
			colour.g = L_e.g + L_l.g + L_s.g + L_c.g + L_d.g;
			colour.b = L_e.b + L_l.b + L_s.b + L_c.b + L_d.b;
}
*/
