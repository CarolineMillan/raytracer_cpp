/*
 * point_light.cpp
 *
 *  Created on: 26 Nov 2019
 *      Author: carolinemillan
 */

#include "point_light.h"
#include "vector.h"
#include <random>
#include <ctime>

PointLight::PointLight()
{
	Light();
}

PointLight::PointLight(Vertex pt, Colour col, Vector dir)
{
	Light();

	direction = dir;
	direction.normalise();
	intensity = col;
	point = pt;
}

PointLight::PointLight(Vertex pt, Colour col)
{
	Light();

	this->intensity = col;
	this->point = pt;

	//srand (time(NULL)); //initialises random seed
	this->direction.x = (rand() % 100000)/100000;
	this->direction.y = (rand() % 100000)/100000;
	this->direction.z = (rand() % 100000)/100000;


}

bool PointLight::get_direction(Vertex &surface, Vector &dir) //unsure about this
{
	/// dir is the vector from pointlight to surface, returns true

	dir.x = surface.x - point.x;
	dir.y = surface.y - point.y;
	dir.z = surface.z - point.z;

	dir.normalise();

	//dir = direction;

	// I think viewing direction (dir) is irrelevant here

	// need to check to see whether there is a clear line between surface and self.point

	// this won't be hard, just check for a closer intersection between pt and surface

	// get light dir

	//Vector surface_vec = {surface.x, surface.y, surface.z};
	//Vector pt_vec = (point.x, point.y, point.z);
	//surface_vec.sub(pt_vec);

	// now surface_vec is the vector from the point light and the surface

	// so check to see if there is an object intersection with anything closer

	// unless I just use the point light for photons, and use direction for the general ray tracing
	// do this for now


	//dir = direction;

	return true;
}

void PointLight::get_intensity(Vertex D, Colour &level)
{
	level = intensity;
	/*
	if (direction.x == 0 && direction.y == 0 && direction.z == 0)
	{
		level = intensity;
	}
	else
	{
		//a function of the angle between light direction and light to object vector
		//I_l= I_i*(L.D)^n

		//get angle between two vectors
		float	cos_theta, d;
		cos_theta = direction.dot(D);
		//cos_theta.normalise();

		float n = 40;

		if(cos_theta >= 90) //is this the right test? radians?
		{
			level.r = 0;
			level.g = 0;
			level.b = 0;
		}
		else
		{
			level.r = intensity.r * pow(direction.dot(D), n);
			level.g = intensity.g * pow(direction.dot(D), n);
			level.b = intensity.b * pow(direction.dot(D), n);
		}
	}
		*/
}
