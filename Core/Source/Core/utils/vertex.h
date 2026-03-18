/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

#pragma once

#include <cmath>

class Vertex {
public:
	float x;
	float y;
	float z;
	float w;

	Vertex()
	{
		x = 0.0;
		y = 0.0;
		z = 0.0;
		w = 1.0f;
	}

	Vertex(float px, float py, float pz, float pw)
	{
		x = px;
		y = py;
		z = pz;
		w = pw;
	}

	Vertex(float px, float py, float pz)
	{
		x = px;
		y = py;
		z = pz;
		w = 1.0f;
	}

	void sub(Vertex &other)
	{
	  x -= other.x;
	  y -= other.y;
	  z -= other.z;
	}

	float distance(Vertex &other) const{
		float a = other.x - x;
		float b = other.y - y;
		float c = other.z - z;
		return sqrt(a*a + b*b + c*c);
	}

};
