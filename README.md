# Ray Tracing with Photon Mapping in C++

This is a ray tracer written in C++ that I wrote as coursework for my Advanced Computer Graphics module at the University of Bath. I have since refactored the code into something more readable.

It outputs a scene to a PPM image, which can be converted to a viewable PNG via your terminal using ImageMagick or a website like [convertio](https://convertio.co/).

## Features

- Whitted Ray Tracing
- supports spheres and Polymesh objects
- supports point lights, shadows, transparent and reflective materials
- phong shading (illumination, not interpoltion (yet...))
- photon mapping (global and caustic) using kd-trees and k-NN searches
- multi object scenes

## Development notes

I originally wrote this in 2019 as coursework.

I have tidied it up (in 2025) before making it public here.

Changes I made:
- restructured the project
    - moved raytracer code from main.cpp to scene.cpp
    - broke up the project into lots of smaller helper functions, rather than the handful of massive ones I had before
    - moved code into subdirectories
    - used a project template to strucure the project
- added smart pointers to avoid memory leaks

## Installation and Running the Path Tracer

(currently just instructions for myself)
build using ```make config=debug```
run from the root of the project using ```./Binaries/macosx-x86_64/Debug/App/app```

## Example Images

## Future Plans:

- [X] smart pointers
- [ ] proper importance sampling (right now I've just hard coded some photons to fire at specular objects, also for k-NN on globalTree)
- [ ] phong interpolation for meshes
- [X] added in constructors and deconstructors
- [ ] helper functions for vector and vertex
- [ ] a cli
- [ ] store photon maps between renders (the whole point of photon mapping)
- [ ] add some features from Peter Shirley's ray tracing in one weekend (BVH, camera model, antialiasing, perhaps noise and procedural generation)

## Acknowledgements

- my lecturer Ken Cameron provided the base code that I worked off for this project, and I've left his comments in 
- I used [The Cherno's project template](https://github.com/TheCherno/ProjectTemplate) to structure the project
- kdtree.h and ```~Scene()``` were AI-assisted, the rest is my own work

## Licence
