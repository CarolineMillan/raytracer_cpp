// quick run using:
// g++ -o mainexecutable main.cpp framebuffer.cpp polymesh.cpp directional_light.cpp sphere.cpp phong.cpp scene.cpp point_light.cpp photon.cpp -lm
// ./mainexecutable
// magick test.ppm test.png

#include <iostream>
#include <ctime>

#include "hit.h"
#include "ray.h"
#include "framebuffer.h"
#include "transform.h"
#include "polymesh.h"
#include "directional_light.h"
#include "point_light.h"
#include "phong.h"
#include "sphere.h"
#include "scene.h"

int main() {

	// Get the timestamp for the current date and time
	time_t timestamp;
	time(&timestamp);

	// Display the date and time represented by the timestamp
	cout << ctime(&timestamp);

    // create a framebuffer
    int width = 128; 
	int height = 128; // 1280
    FrameBuffer *fb = new FrameBuffer(width,height);

	srand (time(NULL)); //initialises random seed

	//create a scene 
	Scene scene = Scene();
	//scene.test();
	scene.teapot_box();
	scene.create_photon_maps();

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

			Hit h = Hit();

            // do a raytrace
			scene.raytrace(ray, colour, depth, 4, h);

            // plot it in the framebuffer
			fb->plotPixel(x, y, colour.r, colour.g, colour.b);
			fb->plotDepth(x,y, depth);
		}
		if (y % (width/10) == 0) {
			time(&timestamp);
			cout << ctime(&timestamp);
		}
		cerr << "Progress: " << y << "/" << height << "\n" << flush;
	}
	// write framebuffer to a ppm file
	fb->writeRGBFile((char *)"test.ppm");

	time(&timestamp);

	// Display the date and time represented by the timestamp
	cout << "Finished at: " << ctime(&timestamp);
	return 0;
}
