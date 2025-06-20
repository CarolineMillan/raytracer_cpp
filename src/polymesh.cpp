/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "polymesh.h"

using namespace std;

PolyMesh::PolyMesh(char *file)
{
  Transform *transform = new Transform();

  this->do_construct(file, transform);
}

PolyMesh::PolyMesh(char *file, Transform *transform)
{
  this->do_construct(file, transform);
}

void PolyMesh::do_construct(char *file, Transform *transform)
{
  int count;
  ifstream meshfile(file);

  cerr << "Opening meshfile: " << file << endl;
  
  if (!meshfile.is_open())
  {
    cerr << "Problem reading meshfile (not found)." << endl;
    meshfile.close();
    exit(-1);
  }

  string line;

  try {
    getline(meshfile, line);
  } catch(ifstream::failure e)
  {
    cerr << "Problem reading meshfile (getline failed)." << endl;
  }

  if (line.compare("kcply") != 0)
  {
    cerr << "Problem reading meshfile (not kcply). [" << line << "]" << endl;
    meshfile.close();
    exit(-1);
  }

  try {
    getline(meshfile, line);
  } catch(ifstream::failure e)
  {
    cerr << "Problem reading meshfile (getline failed)." << endl;
    exit(-1);
  }

  istringstream vertex_iss(line);
  string vertex_element;
  string vertex_label;

  vertex_iss >> vertex_element >> vertex_label >> vertex_count;

  if ((vertex_element.compare("element") != 0)||(vertex_label.compare("vertex") != 0))
  {
    cerr << "Problem reading meshfile (element vertex)."<< endl;
    meshfile.close();
    exit(-1);
  }

  cerr << "Expect " << vertex_count << " vertices." << endl;
  
  try {
    getline(meshfile, line);
  } catch(ifstream::failure e)
  {
    cerr << "Problem reading meshfile (getline failed)." << endl;
    exit(-1);
  }

  istringstream triangle_iss(line);
  string triangle_element;
  string triangle_label;

  triangle_iss >> triangle_element >> triangle_label >> triangle_count;

  if ((triangle_element.compare("element") != 0)||(triangle_label.compare("face") != 0))
  {
    cerr << "Problem reading meshfile (element triangle)."<< endl;
    meshfile.close();
    exit(-1);
  }

  cerr << "Expect " << triangle_count << " triangles." << endl;

  vertex = new Vertex[vertex_count];
  
  triangle = new TriangleIndex[triangle_count];
  face_normal = new Vector[triangle_count];
  vertex_normal = new Vector[vertex_count];

  int i;

  for (i = 0; i < vertex_count; i += 1)
  {
    try {
      getline(meshfile, line);
    } catch(ifstream::failure e)
    {
      cerr << "Problem reading meshfile (getline failed)." << endl;
      exit(-1);
    }

    vertex_iss.clear();
    vertex_iss.str(line);

    vertex_iss >> vertex[i].x >> vertex[i].y >> vertex[i].z;

    vertex[i].w = 1.0f;

    transform->apply(vertex[i]);

    //std::cout << "vertex[i]: " << vertex[i].x << ", " << vertex[i].y << ", " << vertex[i].z << std::endl;
  }

  for (i = 0; i < triangle_count; i += 1)
  {
    try {
      getline(meshfile, line);
    } catch(ifstream::failure e)
    {
      cerr << "Problem reading meshfile (getline failed)." << endl;
      exit(-1);
    }

    triangle_iss.clear();
    triangle_iss.str(line);
    
    triangle_iss >> count >> triangle[i][0] >> triangle[i][1] >> triangle[i][2];
    /*    triangle[i][0] -= 1;
    triangle[i][1] -= 1;
    triangle[i][2] -= 1;*/

    if (count != 3)
    {
      cerr << "Problem reading meshfile (non-triangle present)." << endl;
      exit(-1);
    }

    compute_face_normal(i, face_normal[i]);
  }

  compute_vertex_normals();
  
  meshfile.close();
  cerr << "Meshfile read." << endl;
  next = 0;
}


// Moller-Trumbore
bool PolyMesh::rayTriangleIntersect(const Ray& ray, const Vector &v0, const Vector &v1, const Vector &v2, float &t)
{

  //std::cout << "ray: " << ray << std::endl;
  //printf("a");
  Vector p;
  Vector d;
  Vector e1,e2,h,s,q;
  float a,f,u,v;

  p.x = ray.position.x;
  p.y = ray.position.y;
  p.z = ray.position.z;
  d.x = ray.direction.x;
  d.y = ray.direction.y;
  d.z = ray.direction.z;

  //std::cout << "p: " << p.x << ", " << p.y << ", " << p.z << std::endl;
  //std::cout << "d: " << d.x << ", " << d.y << ", " << d.z << std::endl;

  e1.x = v1.x - v0.x;
  e1.y = v1.y - v0.y;
  e1.z = v1.z - v0.z;
  e2.x = v2.x - v0.x;
  e2.y = v2.y - v0.y;
  e2.z = v2.z - v0.z;

  d.cross(e2,h);
  a = e1.dot(h);

  //e1.cross(e2, h);
  //std::cout << "h: " << h.x << ", " << h.y << ", " << h.z << std::endl;
  //h.normalise();
  //a = d.dot(h);

  //std::cout << "v0: " << v0.x << ", " << v0.y << ", " << v0.z << std::endl;
  //std::cout << "v1: " << v1.x << ", " << v1.y << ", " << v1.z << std::endl;
  //std::cout << "v2: " << v2.x << ", " << v2.y << ", " << v2.z << std::endl;
  //std::cout << "d: " << d.x << ", " << d.y << ", " << d.z << std::endl;
/*
  std::cout << "v1: " << v1.x << ", " << v1.y << ", " << v1.z << std::endl;
  std::cout << "v2: " << v2.x << ", " << v2.y << ", " << v2.z << std::endl;
  std::cout << "d: " << d.x << ", " << d.y << ", " << d.z << std::endl;
  std::cout << "h normalised: " << h.x << ", " << h.y << ", " << h.z << std::endl;
  std::cout << "e1: " << e1.x << ", " << e1.y << ", " << e1.z << std::endl;
  std::cout << "e2: " << e2.x << ", " << e2.y << ", " << e2.z << std::endl;
  std::cout << "a: " << a << std::endl;
*/

  // check if the ray is parallel to the triangle
  if (a > -0.00001f && a < 0.00001f)
  {
    //printf("b");
    //std::cout << "a: " << a << std::endl;
    return false ;
  }

  f = 1/a;
  s.x = p.x - v0.x;
  s.y = p.y - v0.y;
  s.z = p.z - v0.z;

  //std::cout << "f: " << f << std::endl;
  //std::cout << "s: " << s.x << ", " << s.y << ", " << s.z << std::endl;

  u = f * s.dot(h);

  

  s.cross(e1,q);
  //std::cout << "q: " << q.x << ", " << q.y << ", " << q.z << std::endl;
  v = f * d.dot(q);
  //printf("    u=%.6f, v=%.6f\n", u, v);

  if (u < 0.0f || u > 1.0f)
  {
    //printf("c");
    return false;
  }

  if ((v < 0.0f) || ((u + v) > 1.0f))
  {
    //printf("d");
    return false;
  }

  // compute t
 
  t = f * e2.dot(q);

  //std::cout << "t: " << t << std::endl;

  if (t > 0.00001f)
  {
    //printf("e");
    return true; // it's in front ray start
  }

  //printf("f");

  // it's behind you
  return false;
}

void PolyMesh::compute_vertex_normals(void)
{
  int i,j;

  // The vertex_normal array is already zeroed.

  for (i = 0; i < triangle_count; i += 1)
  {
    for (j = 0; j < 3; j += 1)
    {
      vertex_normal[triangle[i][j]].add(face_normal[i]);
    }
  }

  for (i = 0; i < vertex_count; i += 1)
  {
    vertex_normal[i].normalise();
  }
}

void PolyMesh::compute_face_normal(int which_triangle, Vector &normal)
{
  Vector v0v1, v0v2;
  v0v1.x = vertex[triangle[which_triangle][1]].x - vertex[triangle[which_triangle][0]].x;
  v0v1.y = vertex[triangle[which_triangle][1]].y - vertex[triangle[which_triangle][0]].y;
  v0v1.z = vertex[triangle[which_triangle][1]].z - vertex[triangle[which_triangle][0]].z;

  v0v1.normalise();

  v0v2.x = vertex[triangle[which_triangle][2]].x - vertex[triangle[which_triangle][0]].x;
  v0v2.y = vertex[triangle[which_triangle][2]].y - vertex[triangle[which_triangle][0]].y;
  v0v2.z = vertex[triangle[which_triangle][2]].z - vertex[triangle[which_triangle][0]].z;

  v0v2.normalise();

  v0v1.cross(v0v2, normal);
  normal.normalise();
}

void PolyMesh::triangle_intersection(Ray ray, Hit &hit, int which_triangle)
{
  hit.flag = false;

  float ndotdir = face_normal[which_triangle].dot(ray.direction);

  //std::cout << "first ray: " << ray.direction.x << ", " << ray.direction.y << ", " << ray.direction.z << std::endl;
  //std::cout << "face_normal[i]: " << face_normal[which_triangle].x << ", " << face_normal[which_triangle].y << ", " << face_normal[which_triangle].z << std::endl;
  //std::cout << "ndotdir: " << ndotdir << std::endl;

  // ray is parallel to triangle so does not intersect
  if (fabs(ndotdir) < 0.000000001f) return;

  Vector v0,v1,v2;
  v0.x = vertex[triangle[which_triangle][0]].x;
  v1.x = vertex[triangle[which_triangle][1]].x;
  v2.x = vertex[triangle[which_triangle][2]].x;

  v0.y = vertex[triangle[which_triangle][0]].y;
  v1.y = vertex[triangle[which_triangle][1]].y;
  v2.y = vertex[triangle[which_triangle][2]].y;

  v0.z = vertex[triangle[which_triangle][0]].z;
  v1.z = vertex[triangle[which_triangle][1]].z;
  v2.z = vertex[triangle[which_triangle][2]].z;


  Vector o;

  o.x = ray.position.x;
  o.y = ray.position.y;
  o.z = ray.position.z;
  float t,u,v;

  hit.flag =  rayTriangleIntersect(ray, v0, v1, v2, t) ;

  //std::cout << "hit.flag: " << hit.flag << std::endl;

  if (hit.flag == false) return;
  

  // intersection is behind start of ray
  if (t <= 0.0f) return;

  //printf("made it");

  Vertex p;

  p.x = ray.position.x + t * ray.direction.x;
  p.y = ray.position.y + t * ray.direction.y;
  p.z = ray.position.z + t * ray.direction.z;

  hit.t = t;
  hit.what = this;
  hit.position = p;

  hit.normal = face_normal[which_triangle];

  hit.normal.normalise();

  //std::cout << "triangle intersection hit.normal: " << hit.normal.x << ", " << hit.normal.y << ", " << hit.normal.z << std::endl;

  return;
}

void PolyMesh::intersection(Ray ray, Hit &hit)
{
  //printf("in pm intersection");
  Hit current_hit;
  int i;

  hit.flag = false;

  // Check each triangle, find closest if any intersecion

  for(i = 0; i < triangle_count; i += 1)
  {
    //printf("A");
    triangle_intersection(ray, current_hit, i);

    if (current_hit.flag != false) {
      //printf("we have a triangle hit");
      //printf("B");
      if (hit.flag == false) {
        //printf("C");
	      hit = current_hit;

      } else if (current_hit.t < hit.t) {
        //printf("D");
        hit = current_hit;
      }
    }
  }

  if (hit.flag == true) {
    //printf("E");
    //printf("we have a pm hit!");
    if(hit.normal.dot(ray.direction) > 0.0) {
      hit.normal.negate();
    }
  }
  //std::cout << "pm intersection hit.normal: " << hit.normal.x << ", " << hit.normal.y << ", " << hit.normal.z << std::endl;
}
