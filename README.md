# Ray Tracing with Photon Mapping in C++

This is a ray tracer written in C++ that I wrote as courswork for my Advanced Computer Graphics module at the University of Bath.

It outputs a scene to a PPM image, which can be converted to a viewable PNG via your terminal using ImageMagick or a website like [convertio](https://convertio.co/).

## Features
- Whitted Ray Tracing
- supports spheres and Polymesh objects
- supports point lights, shadows, transparent and reflective materials
- phong shading (illumination, not interpoltion (yet...))
- photon mapping (global and caustic) using kd-trees and k-NN searches
- multi object scenes

## Development notes

I originally wrote this in 2019 as coursework and got (76%).

I have tidied it up (in 2025) before making it public here.

Changes I made:
- moved raytracer code from main.cpp to scene.cpp
- broke up the project into lots of smaller helper functions, rather than the handful of massive ones I had before

## Installation and Running the Path Tracer

## Example Images

## Future Plans:
- smart pointers
- proper importance sampling (right now I've just hard coded some photons to fire at specular objects, also for k-NN on globalTree)
- phong interpolation for meshes


## Acknowledgements
- my lecturer Ken Cameron provided the base code that I worked off for this project, and I've left his comments in 

## Licence