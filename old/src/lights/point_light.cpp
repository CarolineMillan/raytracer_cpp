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

	return true;
}

void PointLight::get_intensity(Vertex D, Colour &level)
{
	level = intensity;//.scale(0.1f);
}
