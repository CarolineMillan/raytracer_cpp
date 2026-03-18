/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2018.
 *
 * Do what you like with this code as long as you retain this comment.
 */

// Light is the base class for lights.

#pragma once

#include "../../utils/vertex.h"
#include "../../utils/vector.h"
#include "../../core/colour.h"

class Light {
public:
	Light *next;

	Light()
	{
		next = (Light *)0;
	}

    virtual ~Light(){};

	// Get the direction towards the light at the point on the surface
	// should return true if the surface is in front of the light
	// and false if it's behind and not illuminated.
	virtual bool get_direction(Vertex &surface, Vector &dir)
	{
		return false;
	}

	// Get the intensity of the light in the direction of the surface
	virtual void get_intensity(Vertex &surface, Colour &intensity)
	{

	}
};
