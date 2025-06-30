/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

// Object is the base class for objects.
#ifndef _OBJECT_H_
#define _OBJECT_H_

#include "ray.h"
#include "hit.h"
#include "material.h"

class Object {
public:

	Object *next;
	Object *c_next; //to create a list for caustic photon map
	Material *material;

	Object()
	{
		next = (Object *)0;
		c_next = (Object *)0;
		material = 0;
	}
	
	virtual void intersection(Ray ray, Hit &hit)
	{

	}
};

#endif
