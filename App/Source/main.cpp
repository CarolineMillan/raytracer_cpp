// quick run using:
// g++ -o mainexecutable main.cpp framebuffer.cpp polymesh.cpp directional_light.cpp sphere.cpp phong.cpp scene.cpp point_light.cpp photon.cpp -lm
// ./mainexecutable
// magick test.ppm test.png

#include <iostream>
#include <ctime>
#include <filesystem>

//#include "hit.h"
//#include "ray.h"
//#include "framebuffer.h"
//#include "transform.h"
//#include "polymesh.h"
//#include "directional_light.h"
//#include "point_light.h"
//#include "phong.h"
//#include "sphere.h"
//#include "scene.h"


#include "Core/core/ray.h"
#include "Core/core/framebuffer.h"
#include "Core/core/scene.h"
#include "Core/geometry/hit.h"
#include "Core/utils/transform.h"
#include "Core/geometry/polymesh.h"
#include "Core/geometry/lights/directional_light.h"
#include "Core/geometry/lights/point_light.h"
#include "Core/materials/phong.h"
#include "Core/geometry/sphere.h"

int main() {

	// Get the timestamp for the current date and time
	time_t timestamp;
	time(&timestamp);

	// Display the date and time represented by the timestamp
	cout << ctime(&timestamp);

    // create a framebuffer
    int width = 100; 
	int height = 100; // 1280
    FrameBuffer fb = FrameBuffer(width,height);

	srand (time(NULL)); //initialises random seed

	//create a scene 
	Scene scene = Scene();
	//scene.test();
	scene.teapot_box();
    //scene.cornell_box();
	scene.create_photon_maps();

    // create a ray starting at (0,0,0) to use for the camera
	Ray ray;
	ray.position.x = 0.0001f;
	ray.position.y = 0.0f;
	ray.position.z = 0.0f;

	std::cout << "(width, height) = (" << width << ", " << height << ")" << std::endl;
    std::filesystem::create_directories("images");

    // loop through every pixel in the screen
    for (int y = 0; y < height; y += 1) {
		for (int x = 0; x < width; x += 1) {
            // for antialiasing, take several samples then average the result
            int samples = 4;
            Colour total;
            float depth = 0;
            
            for (int i = 0; i < samples; i++) {
                // fire a ray into the scene, add a random offset
                float fx = ((float)x + (rand() % 1000)/1000.0f)/(float)width;
                float fy = ((float)y + (rand() % 1000)/1000.0f)/(float)height;
                Vector direction;
                ray.direction.x = (fx-0.5f);
                ray.direction.y = (0.5f-fy);
                ray.direction.z = 0.5f;
                ray.direction.normalise();
                Colour colour;
                depth = 0;

                Hit h = Hit();

                // do a raytrace
                scene.raytrace(ray, colour, depth, 4, h);
                total.add(colour);
            }

            float invert = 1.0f / samples;
            total.scale(invert);

            // plot it in the framebuffer
			fb.plotPixel(x, y, total.r, total.g, total.b);
			fb.plotDepth(x,y, depth);
		}
		if (y % (width/10) == 0) {
			time(&timestamp);
			cout << ctime(&timestamp);
		}
		cerr << "Progress: " << y << "/" << height << "\n" << flush;
	}
    // used AI to get a timestamp in filename
    char filename[64];
    std::strftime(filename, sizeof(filename), "images/render_%Y%m%d_%H%M%S.ppm", localtime(&timestamp));
    // end of AI
	// write framebuffer to a ppm file
	fb.writeRGBFile(filename);

	time(&timestamp);

	// Display the date and time represented by the timestamp
	cout << "Finished at: " << ctime(&timestamp);
	return 0;
}
