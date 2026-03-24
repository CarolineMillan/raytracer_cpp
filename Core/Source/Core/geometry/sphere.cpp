/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

#include "sphere.h"
#include <math.h>


Sphere::Sphere(Vertex c, float r)
{
	center = c;
	radius = r;
}

void Sphere::intersection(Ray ray, Hit &hit)
{
	// offset ray by sphere position
	// equivalent to transforming ray into local sphere space
	Vector ro = ray.position - center;

	hit.flag = false;

	float a = (float)ray.direction.dot(ray.direction);
	float b = (float)(2.0 * ray.direction.dot(ro));
	float c = (float)ro.dot(ro) - radius*radius;

	float disc = b*b - 4 * a*c;

	// a negative value indicates no intersection
	if (disc < 0.0f) return; 

	float ds = sqrtf(disc);

	float t0 = (-b - ds) / 2.0f;
	float t1 = (-b + ds) / 2.0f;

	if (t1 < 0.0f) return;

	// if an intersection has been found, record details in hit object
	hit.what = this;

	//smallest root is positive.
	if (t0 > 0.0) {
		hit.t = t0;
        hit.position = ray.at(t0);
        hit.normal = hit.position - center;
		hit.normal.normalise();
        if(hit.normal.dot(ray.direction) > 0.0) {
            hit.normal.negate();
        }
		hit.flag = true;
		return;
	}
	
	hit.t = t1;
    hit.position = ray.at(t1);
    hit.normal = hit.position - center;
	hit.normal.normalise();
    if(hit.normal.dot(ray.direction) > 0.0) {
        hit.normal.negate();
    }
	hit.flag = true;
	return;
}
