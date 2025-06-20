// g++ -g -fsanitize=address -o mainexecutable main.cpp framebuffer.cpp polymesh.cpp directional_light.cpp sphere.cpp phong.cpp -lm
// ./mainexecutable

#include <iostream>

#include "hit.h"
#include "ray.h"
#include "framebuffer.h"
#include "transform.h"
#include "polymesh.h"
#include "directional_light.h"
#include "phong.h"
#include "sphere.h"

// object intersection
void object_intersection(Ray ray, Object *objects, Hit &best_hit)
{
	Object *obj = objects;

	best_hit.flag = false;

	// in object_intersection, before entering the loop:
	//std::cerr << "Object list head: " << objects << "\n";


    //printf("1");
	while(obj != 0) {
		// inside the loop:
		//std::cerr << "  visiting obj=" << obj << " at next=" << obj->next << "\n";	
        //printf("2");
		Hit obj_hit;
		obj_hit.flag=false;

		obj->intersection(ray, obj_hit);
		

		if (obj_hit.flag) {
			//std::cout << "(in obj int) obj_hit.normal: " << obj_hit.normal.x << ", " << obj_hit.normal.y << ", " << obj_hit.normal.z << std::endl;
            //printf("3");
			if (obj_hit.t > 0.0f) {
                //printf("4");
				if (best_hit.flag == false) {
                    //printf("5");
					best_hit.flag = obj_hit.flag;
					best_hit.t = obj_hit.t;
					best_hit.what = obj_hit.what;
					best_hit.position.x = obj_hit.position.x;
					best_hit.position.y = obj_hit.position.y;
					best_hit.position.z = obj_hit.position.z;
					best_hit.normal.x = obj_hit.normal.x;
					best_hit.normal.y = obj_hit.normal.y;
					best_hit.normal.z = obj_hit.normal.z;
					//printf("5b");
					//std::cout << "(in obj int) best_hit.normal: " << best_hit.normal.x << ", " << best_hit.normal.y << ", " << best_hit.normal.z << std::endl;
					/*
						bool flag;
						float t;
						Object *what;
						Vertex position;
						Vector normal;
					*/
				} else if (obj_hit.t < best_hit.t) {
                    //printf("6");
					//best_hit = obj_hit;
					best_hit.flag = obj_hit.flag;
					best_hit.t = obj_hit.t;
					best_hit.what = obj_hit.what;
					best_hit.position.x = obj_hit.position.x;
					best_hit.position.y = obj_hit.position.y;
					best_hit.position.z = obj_hit.position.z;
					best_hit.normal.x = obj_hit.normal.x;
					best_hit.normal.y = obj_hit.normal.y;
					best_hit.normal.z = obj_hit.normal.z;
					//std::cout << "(in obj int) best_hit.normal: " << best_hit.normal.x << ", " << best_hit.normal.y << ", " << best_hit.normal.z << std::endl;
				}
			}
		}
		//printf("5c");
		obj = obj->next;
		//printf("5d");
	}
	//printf("finished");
	//std::cout << "(in obj int) best_hit.normal: " << best_hit.normal.x << ", " << best_hit.normal.y << ", " << best_hit.normal.z << std::endl;
	return;
}

/*
void raytrace(Ray ray, Object *objects, Light *lights, Colour &colour, float &depth, int ref_limit) //gives depth and colour of hit point
{
	// first step, find the closest primitive

	Hit shadow_hit;
	Hit best_hit;
	printf("in 1");
	object_intersection(ray, objects, best_hit);
	printf("out 1");

	// if we found a primitive then compute the colour we should see
	if(best_hit.flag)
	{
		//std::cout << "best_hit.normal: " << best_hit.normal.x << ", " << best_hit.normal.y << ", " << best_hit.normal.z << std::endl;
		//printf("we have a hit!");
		best_hit.what->material->compute_base_colour(colour); //this gives 0?
		depth = best_hit.t;
		Light *light = lights;
		int no_of_reflections = 5; //depth counter to stop infinite recursion

		//std::cerr << "Light list head: " << lights << "\n";

		while (light != (Light *)0) //there is a light source
		{
			//std::cerr << "  visiting light=" << light << " at next=" << light->next << "\n";

            //printf("1");
			Vector viewer;
			Vector ldir;

			viewer.x = -best_hit.position.x; //viewing ray position
			viewer.y = -best_hit.position.y;
			viewer.z = -best_hit.position.z;
			viewer.normalise();

			bool lit;
			lit = light->get_direction(best_hit.position, ldir); //is it hit by the light?

			if(ldir.dot(best_hit.normal)>0)	lit=false;//light is facing wrong way.

			// make sure ldir is the vector *from* point *to* light
			// (or flip its sign consistently in get_direction)
			//if (best_hit.normal.dot(ldir) <= 0) lit = false;

			if(lit)
			{
                //printf("3");
				Ray shadow_ray;

				shadow_ray.direction.x = -ldir.x;
				shadow_ray.direction.y = -ldir.y;
				shadow_ray.direction.z = -ldir.z;
				shadow_ray.position.x = best_hit.position.x + (0.0001f * shadow_ray.direction.x);
				shadow_ray.position.y = best_hit.position.y + (0.0001f * shadow_ray.direction.y);
				shadow_ray.position.z = best_hit.position.z + (0.0001f * shadow_ray.direction.z);


				//printf("in 2");
				object_intersection(shadow_ray, objects, shadow_hit);
				//printf("out 2");
				// if (shadow_hit.flag) {std::cout << "shadow_hit.normal: " << shadow_hit.normal.x << ", " << shadow_hit.normal.y << ", " << shadow_hit.normal.z << std::endl;}
				if(shadow_hit.flag==true)
				{
                    //printf("4");
					if (shadow_hit.t < 1000000000.0f)
					{
                        //printf("4b");
						//std::cout << "shadow: " << best_hit.what->material->kr << std::endl;

						lit = false; //there's a shadow so no lighting, if realistically close
					}
				}
			}

			if (lit) //do colour
			{
                //printf("5");
				Colour intensity;
				Colour scaling;

				light->get_intensity(best_hit.position, scaling);

				best_hit.what->material->compute_light_colour(viewer, best_hit.normal, ldir, intensity);

				//intensity.scale(scaling); //scale by scaling to get intensity
				intensity.r *= scaling.r;
				intensity.g *= scaling.g;
				intensity.b *= scaling.b;
				intensity.a *= scaling.a;

				//std::cout << "Intensity, r: " << intensity.r << ", g: " << intensity.g << ", b: " << intensity.b << std::endl;
				//std::cout << "Colour 5 before, r: " << colour.r << ", g: " << colour.g << ", b: " << colour.b << std::endl;

				//colour.add(intensity); //add intensity to the colour
				colour.r += intensity.r;
				colour.g += intensity.g;
				colour.b += intensity.b;

                //std::cout << "Colour 5 after, r: " << colour.r << ", g: " << colour.g << ", b: " << colour.b << std::endl;
			}

			light = light->next; //move on to next light
		}

		// TODO: compute refraction ray if material supports it.
		if(best_hit.what->material->transparent)
		{
            //printf("6");
			Ray T;
			Colour trans_colour;
			float trans_depth, eta, cos_theta_t, cos_theta_i, rpar, rper, kr, kt;

			eta = best_hit.what->material->eta;

			Vector n;
			n.x =  best_hit.normal.x;
			n.y =  best_hit.normal.y;
			n.z =  best_hit.normal.z;
			cos_theta_i = n.dot({-1*ray.direction.x, -1*ray.direction.y, -1*ray.direction.z});
			cos_theta_t = sqrt(1 - ( 1 / (pow(eta,2) ) ) * ( 1 - pow(cos_theta_i,2) ) );

			/*
			//std::cout << "eta: " << eta << std::endl;
			//std::cout << "best_hit.normal: " << best_hit.normal.x << ", " << best_hit.normal.y << ", " << best_hit.normal.z << std::endl;
			//std::cout << "n: " << n.x << ", " << n.y << ", " << n.z << std::endl;
			//std::cout << "cos_theta_i: " << cos_theta_i << std::endl;
			//std::cout << "cos_theta_t: " << cos_theta_i << std::endl;
			
			//tests for total internal reflection
			if(cos_theta_t >= 0) 
			{
                //printf("THIS 7");
				//get the refracted ray
				T.direction.x = (1/eta)*ray.direction.x - (cos_theta_t - (1/eta)*cos_theta_i)*best_hit.normal.x;
				T.direction.y = (1/eta)*ray.direction.y - (cos_theta_t - (1/eta)*cos_theta_i)*best_hit.normal.y;
				T.direction.z = (1/eta)*ray.direction.z - (cos_theta_t - (1/eta)*cos_theta_i)*best_hit.normal.z;

				T.position.x = best_hit.position.x + 0.001*T.direction.x;
				T.position.y = best_hit.position.y + 0.001*T.direction.y;
				T.position.z = best_hit.position.z + 0.001*T.direction.z;

				//raytrace the refracted ray
				raytrace(T, objects, lights, trans_colour, trans_depth, ref_limit);

				//fresnel equations
				rpar = (eta*cos_theta_i - cos_theta_t)/(eta*cos_theta_i + cos_theta_t);
				rper = (cos_theta_i - eta*cos_theta_t)/(cos_theta_i + eta*cos_theta_t);
				kr = 0.5*( pow(rpar,2) + pow(rper,2) ); //do we use this in reflection?
				kt = 1-kr;

				best_hit.what->material->kr = kr; //this is the value needed for reflection in this surface at this point
				best_hit.what->material->kt = kt;

				//scale trans_colour by kt
				trans_colour.r *= kt; //use scale function?
				trans_colour.g *= kt;
				trans_colour.b *= kt;

				//std::cout << "Trans Colour, r: " << trans_colour.r << ", g: " << trans_colour.g << ", b: " << trans_colour.b << std::endl;

				//std::cout << "Colour 7 before, r: " << colour.r << ", g: " << colour.g << ", b: " << colour.b << std::endl;
				//add the refracted colour to the colour
				//colour.add(trans_colour);
				colour.r += trans_colour.r;
				colour.g += trans_colour.g;
				colour.b += trans_colour.b;

                //std::cout << "Colour 7 after , r: " << colour.r << ", g: " << colour.g << ", b: " << colour.b << std::endl;
			}

		}
		// TODO: compute reflection ray if material supports it.
		if(best_hit.what->material->reflective) //material supports reflection
		{
            //printf("8");
			Ray R;
			float kr, ref_depth;
			Hit rhit;
			Colour ref_colour;

			//calculate reflected eye ray
			R.direction.x = ray.direction.x - 2*((best_hit.normal).dot(ray.direction))*best_hit.normal.x;
			R.direction.y = ray.direction.y - 2*((best_hit.normal).dot(ray.direction))*best_hit.normal.y;
			R.direction.z = ray.direction.z - 2*((best_hit.normal).dot(ray.direction))*best_hit.normal.z;

			R.position.x = best_hit.position.x + 0.001*R.direction.x;
			R.position.y = best_hit.position.y + 0.001*R.direction.y;
			R.position.z = best_hit.position.z + 0.001*R.direction.z;

			//printf("in 3");
			//tests if there is an object to be reflected in the surface
			object_intersection(R, objects, rhit); 
			//printf("out 3");

			//if (rhit.flag) {std::cout << "raytrace reflective intersection rhit.flag: " << rhit.flag << ", rhit.normal: " << rhit.normal.x << ", " << rhit.normal.y << ", " << rhit.normal.z << std::endl;}
			

			ref_limit -=1;

			if(ref_limit<0)
			{
                //printf("9");
				return;
			}

			if(rhit.flag)
			{
                //printf("10");
				//get colour of reflection 
				raytrace(R, objects, lights, ref_colour, ref_depth, ref_limit); 

				kr = best_hit.what->material->kr;

				//scale by kr
				ref_colour.r *= kr; 
				ref_colour.g *= kr;
				ref_colour.b *= kr;

				//std::cout << "rhit Colour, r: " << colour.r << ", g: " << colour.g << ", b: " << colour.b << std::endl;

				//add reflective colour to the colour
				colour.add(ref_colour);

                //std::cout << "rhit2 Colour, r: " << colour.r << ", g: " << colour.g << ", b: " << colour.b << std::endl;
			}
		}
	}
	else
	{
		depth = 7.0f;
		colour.r = 0.0f;
		colour.g = 0.0f;
		colour.b = 0.0f;
	}
}
*/

// gives depth and colour of hit pt
void raytrace(Ray ray, Object *objects, Light *lights, Colour &colour, float &depth, int ref_limit) {

	if (ref_limit < 0) return;

	Hit best_hit;
	Hit shadow_hit;
	object_intersection(ray, objects, best_hit);

	// if we found a primitive then compute the colour we should see
	if(best_hit.flag) {

		// get colour of the material we've hit
		// compute base colour is a method on the materials class

		best_hit.what->material->compute_base_colour(colour);

		// get the depth of the object
		depth = best_hit.t;

		// get a handle on the first light -- check this
		Light *light = lights;


		// add shadow, refraction and reflection

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

				object_intersection(shadow_ray, objects, shadow_hit);

				//there's a shadow so no lighting, if realistically close
				if(shadow_hit.flag && shadow_hit.t < 1000000000.0f) lit = false;
			}

			//do colour
			if (lit) 
				{
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
		// TODO: compute refraction ray if material supports it.

		//ref_limit--;
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
			/*
			//tests for total internal reflection
			if(cos_theta_t >= 0) {
				
				//printf("cos_theta_t >= 0");
				//get the refracted ray
				T.direction = (1/eta)*ray.direction - (cos_theta_t - (1/eta)*cos_theta_i)*best_hit.normal;

				T.position.x = best_hit.position.x + 0.001*T.direction.x;
				T.position.y = best_hit.position.y + 0.001*T.direction.y;
				T.position.z = best_hit.position.z + 0.001*T.direction.z;

				ref_limit -= 1;

				//std::cout << "ref_limit2: " << ref_limit << ". " << std::endl;

				//raytrace the refracted ray
				raytrace(T, objects, lights, trans_colour, trans_depth, ref_limit);

				//fresnel equations
				rpar = (eta*cos_theta_i - cos_theta_t)/(eta*cos_theta_i + cos_theta_t);
				rper = (cos_theta_i - eta*cos_theta_t)/(cos_theta_i + eta*cos_theta_t);
				kr = 0.5*( pow(rpar,2) + pow(rper,2) ); //do we use this in reflection?
				kt = 1-kr;

				//best_hit.what->material->kr = kr; //this is the value needed for reflection in this surface at this point
				//best_hit.what->material->kt = kt;

				//scale trans_colour by kt
				trans_colour.r *= kt; //use scale function?
				trans_colour.g *= kt;
				trans_colour.b *= kt;

				//add the refracted colour to the colour
				colour.add(trans_colour);
			}
			*/
			if (sin2_t <= 1.0f) {
				cos_theta_t = sqrt(1.0 - ( 1.0 / (pow(eta,2.0) ) ) * ( 1.0 - pow(cos_theta_i,2.0) ) );
				//cos_theta_t = sqrt(1.0f - sin2_t);

				T.direction = (1/eta)*ray.direction - (cos_theta_t - (1/eta)*cos_theta_i)*best_hit.normal;

				T.position.x = best_hit.position.x + 0.001*T.direction.x;
				T.position.y = best_hit.position.y + 0.001*T.direction.y;
				T.position.z = best_hit.position.z + 0.001*T.direction.z;

				ref_limit -= 1;

				//std::cout << "ref_limit2: " << ref_limit << ". " << std::endl;

				//raytrace the refracted ray
				raytrace(T, objects, lights, trans_colour, trans_depth, ref_limit);

				//fresnel equations
				rpar = (eta*cos_theta_i - cos_theta_t)/(eta*cos_theta_i + cos_theta_t);
				rper = (cos_theta_i - eta*cos_theta_t)/(cos_theta_i + eta*cos_theta_t);
				kr = 0.5*( pow(rpar,2) + pow(rper,2) ); //do we use this in reflection?
				kt = 1-kr;

				best_hit.what->material->kr = kr; //this is the value needed for reflection in this surface at this point
				best_hit.what->material->kt = kt;


				//std::cout << "n: " << n.x << ", " << n.y << ", " << n.z << std::endl;
				//std::cout << "ray.direction: " << ray.direction.x << ", " << ray.direction.y << ", " << ray.direction.z <<std::endl;
				//std::cout << "cos_theta_i: " << cos_theta_i << std::endl;
				//std::cout << "cos_theta_t: " << cos_theta_t << std::endl;
				//std::cout << "rpar: " << rpar << std::endl;
				//std::cout << "rper: " << rper << std::endl;
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
		if(best_hit.what->material->reflective) //material supports reflection
		{
			Ray R;
			float kr, ref_depth;
			Hit rhit;
			Colour ref_colour;

			//calculate reflected eye ray
			R.direction = ray.direction - 2*((best_hit.normal).dot(ray.direction))*best_hit.normal;

			R.position.x = best_hit.position.x + 0.001*R.direction.x;
			R.position.y = best_hit.position.y + 0.001*R.direction.y;
			R.position.z = best_hit.position.z + 0.001*R.direction.z;

			//std::cout << "ref_limit3: " << ref_limit << ". " << std::endl;

			//tests if there is an object to be reflected in the surface
			object_intersection(R, objects, rhit); 

			ref_limit -= 1;

			if(ref_limit<0)
			{
				return;
			}

			if(rhit.flag)
			{

				//std::cout << "ref_limit4: " << ref_limit << ". " << std::endl;
				//get colour of reflection 

				raytrace(R, objects, lights, ref_colour, ref_depth, ref_limit); 

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
		// get these between 0 and 1 
		//printf("in background colour bit");
		
		//colour.r = 163.0/255.0;
		//colour.g = 249.0/255.0;
		//colour.b = 255.0/255.0;
		depth = 7.0f;
		
	}

	return;
}

		/*
			◦ If (there is a hit) then
			◦ Generate a shadow ray
			◦ If (there is no hit between the shadow ray and a light) then
			◦ c = c + shading()
			◦ Generate a secondary ray (reflection or refraction ray) // increase the ray depth +1
			◦ Else
			◦ c = c + background color
			◦ Set the pixel color with c
		*/

		// while we have a light source

        // if the object is illuminated by a directional light source
            // create a shadow ray
            // do colour for the hit point

        // compute refraction ray if material supports it

        // compute reflection ray if material supports it

/*

• For each pixel 
◦ Color c = (0, 0, 0)
◦ Generate a primary ray (with depth 0)
◦ While (depth < d) 
◦ Find the closest intersection point between the ray and objects
◦ If (there is a hit) then
◦ Generate a shadow ray
◦ If (there is no hit between the shadow ray and a light) then
◦ c = c + shading()
◦ Generate a secondary ray (reflection or refraction ray) // increase the ray depth +1
◦ Else
◦ c = c + background color
◦ Set the pixel color with c

*/


int main() {

    // create a framebuffer
    int width = 1280;
	int height = 1280;
    FrameBuffer *fb = new FrameBuffer(width,height);

	//creates a light source
    DirectionalLight *dl = new DirectionalLight(Vector(0.01f, 0.01f, 1.0f),Colour(1.0f, 1.0f, 1.0f, 1.0f)); 
	//lights->next = dl;
	dl->next = nullptr;  

	// create scene

	//A NICE SCENE -- TEAPOT AND SPHERE INSIDE AN OPEN BOX

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
	Phong glass; 
	glass.ambient = Colour(244.0/255.0, 250.0/255.0, 252.0/255.0, 255.0/255.0);
	glass.diffuse = Colour(244.0/255.0, 250.0/255.0, 252.0/255.0, 255.0/255.0);
	glass.specular = Colour(255.0/255.0, 255.0/255.0, 255.0/255.0, 255.0/255.0);
	glass.power = 40.0f;

    Phong bp1; 
	// rgb(211, 141, 255)
	bp1.ambient = Colour(211.0/255.0, 141.0/255.0, 255.0/255.0, 255.0/255.0);
	bp1.diffuse = Colour(211.0/255.0, 141.0/255.0, 255.0/255.0, 255.0/255.0);
	bp1.specular = Colour(255.0/255.0, 255.0/255.0, 255.0/255.0, 255.0/255.0);
	bp1.power = 40.0f;

	// rgb(255, 252, 230)
	Phong bp2;
	bp2.ambient = Colour(255.0/255.0, 252.0/255.0, 230.0/255.0, 255.0/255.0);
	bp2.diffuse = Colour(255.0/255.0, 252.0/255.0, 230.0/255.0, 255.0/255.0);
	bp2.specular = Colour(255.0/255.0, 255.0/255.0, 255.0/255.0, 255.0/255.0);
	bp2.power = 40.0f;

	// rgb(247, 198, 198)
	Phong bp3;
	bp3.ambient = Colour(247.0/255.0, 198.0/255.0, 198.0/255.0, 255.0/255.0);
	bp3.diffuse = Colour(247.0/255.0, 198.0/255.0, 198.0/255.0, 255.0/255.0);
	bp3.specular = Colour(255.0/255.0, 255.0/255.0, 255.0/255.0, 255.0/255.0);
	bp3.power = 40.0f;

	// rgb(199, 247, 198)
	Phong bp4;
	bp4.ambient = Colour(199.0/255.0, 247.0/255.0, 198.0/255.0, 255.0/255.0);
	bp4.diffuse = Colour(199.0/255.0, 247.0/255.0, 198.0/255.0, 255.0/255.0);
	bp4.specular = Colour(255.0/255.0, 255.0/255.0, 255.0/255.0, 255.0/255.0);
	bp4.power = 40.0f;

	// rgb(198, 242, 247)
	Phong bp5;
	bp5.ambient = Colour(198.0/255.0, 242.0/255.0, 247.0/255.0, 255.0/255.0);
	bp5.diffuse = Colour(198.0/255.0, 242.0/255.0, 247.0/255.0, 255.0/255.0);
	bp5.specular = Colour(255.0/255.0, 255.0/255.0, 255.0/255.0, 255.0/255.0);
	bp5.power = 40.0f;

	// rgb(255, 176, 249)
	Phong bp6;
	bp6.ambient = Colour(255.0/255.0, 176.0/255.0, 249.0/255.0, 255.0/255.0);
	bp6.diffuse = Colour(255.0/255.0, 176.0/255.0, 249.0/255.0, 255.0/255.0);
	bp6.specular = Colour(255.0/255.0, 255.0/255.0, 255.0/255.0, 255.0/255.0);
	bp6.power = 40.0f;

	// rgb(180, 169, 245)
	Phong bp7;
	bp7.ambient = Colour(180.0/255.0, 169.0/255.0, 245.0/255.0, 255.0/255.0);
	bp7.diffuse = Colour(180.0/255.0, 169.0/255.0, 245.0/255.0, 255.0/255.0);
	bp7.specular = Colour(255.0/255.0, 255.0/255.0, 255.0/255.0, 255.0/255.0);
	bp7.power = 40.0f;

    pm->material = &bp1;
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

	
	fl->material = &bp2;
    fl->material->transparent = false;
    fl->material->reflective = false;

	ce->material = &bp3;
    ce->material->transparent = false;
    ce->material->reflective = false;

	w1->material = &bp4;
    w1->material->transparent = false;
    w1->material->reflective = false;

	w2->material = &bp5;
    w2->material->transparent = false;
    w2->material->reflective = false;
	
	w3->material = &bp6;
    w3->material->transparent = false;
    w3->material->reflective = false;

	w4->material = &bp7;
    w4->material->transparent = false;
    w4->material->reflective = false;

	sphere->next = fl;
	fl->next = ce;
	ce->next = w1;
	w1->next = w2;
	w2->next = w3;
	w3->next = nullptr;
	//w4->next = nullptr;


    // create a ray starting at (0,0,0) to use for the camera
	Ray ray;
	ray.position.x = 0.0001f;
	ray.position.y = 0.0f;
	ray.position.z = 0.0f;

    // loop through every pixel in the screen
    for (int y = 0; y < height; y += 1) {
		for (int x = 0; x < width; x += 1) {

			// fire a ray into the scene
            float fx = (float)x/(float)width;
			float fy = (float)y/(float)height;

			Vector direction;
			ray.direction.x = (fx-0.5f);
			ray.direction.y = (0.5f-fy);
			ray.direction.z = 0.5f;
			ray.direction.normalise();

			Colour colour;
			float depth = 0;

            // do a raytrace
			raytrace(ray, pm, dl, colour, depth, 4);

            // plot it in the framebuffer
			fb->plotPixel(x, y, colour.r, colour.g, colour.b);
			fb->plotDepth(x,y, depth);
		}
		std::cerr << "*" << flush;
	}
	// write framebuffer to a ppm file
	fb->writeRGBFile((char *)"test.ppm");
	printf("finished");
	return 0;
}

// raytrace function

    // find the closest primitive/object

    // if we found a primitive then compute the colour we should see

        // while we have a light source

        // if the object is illuminated by a directional light source
            // create a shadow ray
            // do colour for the hit point

        // compute refraction ray if material supports it

        // compute reflection ray if material supports it


/* PHONG INTERPOLATION

Surface rendering is done with the help of Phong shading in the following manner:

Determine the average unit normal vector at each polygon vertex.
              So, for n polygons => summation of Ni/ | summation of Ni | (where i is initialized from 1 to N)

Linearly interpolate the vertex normal over the surfaces of the polygon.   N= (y-y2)/(y1-y2) N1 + (y1-y)/(y1-y2) N2

Interpolation of the surface normal
                                                      

By applying the illumination model along each scan we have to determine the projected pixel intensities of the surface points.

*/





	/*
	// rgb(211, 141, 255)
	// rgb(244, 250, 252)
	Phong bp1; 
	bp1.ambient = Colour(244.0/255.0, 250.0/255.0, 252.0/255.0, 255.0/255.0);
	bp1.diffuse = Colour(244.0/255.0, 250.0/255.0, 252.0/255.0, 255.0/255.0);
	bp1.specular = Colour(255.0/255.0, 255.0/255.0, 255.0/255.0, 255.0/255.0);
	bp1.power = 40.0f;

	// rgb(199, 247, 198)
	Phong bp4;
	bp4.ambient = Colour(199.0/255.0, 247.0/255.0, 198.0/255.0, 255.0/255.0);
	bp4.diffuse = Colour(199.0/255.0, 247.0/255.0, 198.0/255.0, 255.0/255.0);
	bp4.specular = Colour(255.0/255.0, 255.0/255.0, 255.0/255.0, 255.0/255.0);
	bp4.power = 40.0f;

	// rgb(255, 176, 249)
	Phong bp6;
	bp6.ambient = Colour(255.0/255.0, 176.0/255.0, 249.0/255.0, 255.0/255.0);
	bp6.diffuse = Colour(255.0/255.0, 176.0/255.0, 249.0/255.0, 255.0/255.0);
	bp6.specular = Colour(255.0/255.0, 255.0/255.0, 255.0/255.0, 255.0/255.0);
	bp6.power = 40.0f;

	Vertex v;
	v.x = 0.0f;
	v.y = 0.0f;
	v.z = 5.0f;
    Sphere *sphere = new Sphere(v, 1.0f);

    sphere->material = &bp1;
    sphere->material->transparent = true;
    sphere->material->reflective = true;
	sphere->material->eta = 1.5; //glass refractive index
	sphere->material->kr = 0.4;

	PolyMesh *w1 = new PolyMesh((char *)"wall1.ply");

	w1->material = &bp4;
    w1->material->transparent = false;
    w1->material->reflective = false;

	PolyMesh *w3 = new PolyMesh((char *)"wall3.ply");
	
	w3->material = &bp6;
    w3->material->transparent = false;
    w3->material->reflective = false;

	sphere->next = w1;
	w1->next = w3;
	w3->next = nullptr;
	*/