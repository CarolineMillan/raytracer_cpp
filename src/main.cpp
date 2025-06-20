// g++ -o mainexecutable main.cpp framebuffer.cpp polymesh.cpp directional_light.cpp sphere.cpp phong.cpp scene.cpp -lm
// ./mainexecutable

#include <iostream>

#include "hit.h"
#include "ray.h"
#include "framebuffer.h"
#include "transform.h"
#include "polymesh.h"
#include "directional_light.h"
#include "phong.h"
#include "sphere.h"
#include "scene.h"

int main() {

    // create a framebuffer
    int width = 1280;
	int height = 1280;
    FrameBuffer *fb = new FrameBuffer(width,height);

	//create a scene 
	Scene scene = Scene();
	//scene.test();
	scene.teapot_box();

    // create a ray starting at (0,0,0) to use for the camera
	Ray ray;
	ray.position.x = 0.0001f;
	ray.position.y = 0.0f;
	ray.position.z = 0.0f;

	std::cout << "(width, height) = (" << width << ", " << height << ")" << std::endl;

    // loop through every pixel in the screen
    for (int y = 0; y < height; y += 1) {
		for (int x = 0; x < width; x += 1) {
			// fire a ray into the scene
            float fx = (float)x/(float)width;
			float fy = (float)y/(float)height;
			Vector direction;
			ray.direction.x = (fx-0.5f);
			ray.direction.y = (0.5f-fy);
			ray.direction.z = 0.5f;
			ray.direction.normalise();
			Colour colour;
			float depth = 0;

            // do a raytrace
			scene.raytrace(ray, colour, depth, 4);

            // plot it in the framebuffer
			fb->plotPixel(x, y, colour.r, colour.g, colour.b);
			fb->plotDepth(x,y, depth);
		}
		std::cerr << "*" << flush;
	}
	// write framebuffer to a ppm file
	fb->writeRGBFile((char *)"test.ppm");
	printf("finished");
	return 0;
}