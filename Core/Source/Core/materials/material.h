/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2018.
 *
 * Do what you like with this code as long as you retain this comment.
 */

// Material is the base class for materials.

#pragma once

#include "../utils/vector.h"
#include "../core/colour.h"
#include <cstdio>

#include <iostream>
using namespace std;

class Material {
public:

    Material(){};

    virtual ~Material(){};

	virtual void compute_base_colour(Colour &result)
	{
		result.r = 0.0f;
		result.g = 0.0f;
		result.b = 0.0f;
	}
	virtual void compute_light_colour(Vector &viewer, Vector &normal, Vector &ldir, Colour &result)
	{
		result.r = 0.0f;
		result.g = 0.0f;
		result.b = 0.0f;
	}
    virtual Colour get_diffuse_BRDF() {
        return Colour();
    };

    virtual bool is_reflective() { return false; };
    virtual bool is_transparent() { return false; };
    virtual float get_eta() {return 1.0f;};
    virtual float get_kt() {return 0.0f;};
    virtual float get_kr() {return 0.0f;};
    virtual Colour get_tint() {return Colour(1.0f, 1.0f, 1.0f, 1.0f);}; 
};
