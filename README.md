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
- removed gather_diffuse_reflection (this sped up the photon mapping and made it about 4x faster) (there was something else you did that drastically increased the speed -- look back and check what it was)
- added in gamma correction (uesd sqrt as an approximation for the actual formula)
- added in antialiasing
- fixed Phong normal interpolation (most of it was there, it was just adding the calculations into the final normal calculation) to make polymeshes smooth. I made this an optional thing via a boolean on PolyMesh. [scratchapixel](https://scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/moller-trumbore-ray-triangle-intersection.html) was useful for finding the fix here. 

## Installation and Running the Path Tracer

(currently just instructions for myself)
build using ```make config=debug```
run from the root of the project using ```./Binaries/macosx-x86_64/Debug/App/app```

## Example Images

## Future Plans:

- [X] smart pointers
- [ ] proper importance sampling (right now I've just hard coded some photons to fire at specular objects, also for k-NN on globalTree)
- [X] phong interpolation for meshes
- [X] added in constructors and deconstructors
- [ ] helper functions for vector and vertex
- [ ] a cli
- [ ] store photon maps between renders (the whole point of photon mapping)
- [ ] add some features from Peter Shirley's ray tracing in one weekend (BVH, camera model, antialiasing, perhaps noise and procedural generation)
    - [X] antialiasing
    - [X] gamma correction
    - [X] more material subclasses
    - [ ] BVH
    - [ ] area lights
     - [ ] perlin noise (maybe look at procedural generation? this is more long term)
- [X] delete ```BRDF_s```, which should be used as the specular BRDF for the caustic photon map. I used the diffuse BRDF here instead, as an approximation/simplification. Perhaps add in specular BRDF in future.
- [X] add in other material subclasses
- [ ] sort ```use namespace std``` and ```pragma once``` statements

## Acknowledgements

- my lecturer Ken Cameron provided the base code that I worked off for this project, and I've left his comments in 
- I used [The Cherno's project template](https://github.com/TheCherno/ProjectTemplate) to structure the project
- kdtree.h and ```~Scene()``` were AI-assisted, the rest is my own work

## Licence


## notes

- whitted ray tracer part uses fresnel equations, but the photon mapping uses russian roulette to decide on reflection vs transmission. fresnel equations are physically correct and russian roulette uses randomness and is considered simpler. Maybe change this to make both parts of the ray tracer consistant with each other.
- Phong Illumination model is the ambient + diffuse + specular equation that you use to calculate the colour. Phong Shading (or normal interpolation) is the inerpolation you do using barycentric coordinates on triangles to smooth the polymesh. Two different concepts named after the same guy.
- glass currently reflected light back at the light source (currently a point light, a mathematical point that isn't really a visible light source), so there's a caustic pattern that is physically correct but looks unnatural because it is usually reflected back onto an area light source and therefore obscured. The fix is to add emmisive materials/area lights, which is on the to do list.
- global photon map currently does not bounce, so there is no colour bleeding. Adding in diffuse reflection is on the to do list. the ```break;``` in ```photon_trace()``` stops the photon from being reflected any further. you need to fix ```gather_diffuse_reflection()``` and ```g_russian_roulette``` to get colour bleeding working.



materials stuff:
- added in get_diffuse_BRDF() method on Material
- added in is_reflective() and is_transparent() methods on Material (this means that each material has it's own built in property, and they can't be set manually in the scene setup ie once i've defined a material it can't be tampered with outside of the class)
- delete BRDF_d on Material
- delete reflective and transparent on Material
- created Glass and Metal child classes of Material
- made hard coded importance sampling a circle rather than a square -- still need to do it properly
- added a tint to glass and scaled for it whenever refraction takes place
- made a better constructor for Phong
