/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2018.
 *
 * Do what you like with this code as long as you retain this comment.
 */

#pragma once

#include "colour.h"
#include "ray.h"
#include "object.h"
#include "light.h"
#include "hit.h"
#include "phong.h"

// Scene is a class that is used to build a scene database of objects
// and lights and then trace a ray through it.

class Scene {
public:

	Object *object_list;
	Light *light_list;
	Phong glass;
	Phong  mat_sphere;
  	Phong  mat_wall1;
  	Phong  mat_wall3;	
	Phong mat_pm;
	Phong mat_wall2;
	Phong mat_wall4;
	Phong mat_wall5;
	Phong mat_wall6;
	Phong mat_wall7;

	Scene();

	void teapot_box();

	void test();

	// Trace a ray through the scene and find the closest if any object
	// intersection in front of the ray.
	void object_intersection(Ray ray, Hit &hit);
	
	// Trace a ray through the scene and return its colour. This function
	// is the one that should recurse down the reflection/refraction tree.
	void raytrace(Ray ray, Colour &colour, float &depth, int ref_limit);

};
