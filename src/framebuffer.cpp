/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

#include <stdlib.h>
#include <float.h>
#include <iostream>
#include <fstream>
using namespace std;

#include "framebuffer.h"

FrameBuffer::FrameBuffer(int w, int h)
{
  if ((w<0)||(w>2048)||(h<0)||(h>2048))
  {
    exit(-1);
  }

  this->framebuffer = new Pixel[w * h];

  if (!this->framebuffer)
  {
    exit(-2);
  }

  this->width = w;
  this->height = h;

  for (int i = 0; i < w*h; i += 1)
  {
    this->framebuffer[i].red   = 0.0f;
    this->framebuffer[i].green = 0.0f;
    this->framebuffer[i].blue  = 0.0f;
    this->framebuffer[i].depth = 0.0f;
  }
  
  return;
}

int FrameBuffer::plotPixel(int x, int y, float red, float green, float blue)
{
  if ((x<0)||(x>=this->width)||(y<0)||(y>=this->height))
  {
    printf("fb returned -1");
    return -1;
  }

  //  if ((red > 1.0f) || (red < 0.0f)) cerr<<"out of range\n";
  
  this->framebuffer[y * this->width + x].red = red;
  this->framebuffer[y * this->width + x].green = green;
  this->framebuffer[y * this->width + x].blue = blue;

  return 0;
}

int FrameBuffer::plotDepth(int x, int y, float depth)
{
  if ((x<0)||(x>=this->width)||(y<0)||(y>=this->height))
  {
    return -1;
  }

  this->framebuffer[y * this->width + x].depth = depth;

  return 0;
}

int FrameBuffer::getDepth(int x, int y, float &depth)
{
  if ((x<0)||(x>=this->width)||(y<0)||(y>=this->height))
  {
    return -1;
  }
  
  depth = this->framebuffer[y * this->width + x].depth;

  return 0;
}

int FrameBuffer::getPixel(int x, int y, float &red, float &green, float &blue)
{
  if ((x<0)||(x>=this->width)||(y<0)||(y>=this->height))
  {
    return -1;
  }

  red = this->framebuffer[y * this->width + x].red;
  green = this->framebuffer[y * this->width + x].green;
  blue = this->framebuffer[y * this->width + x].blue;

  return 0;
}

int FrameBuffer::writeRGBFile(char *filename)
{
  float min = 0.0f;
  float max = 0.0f;

  ofstream outfile;

  outfile.open(filename, ofstream::binary);

  if (!outfile.is_open())
  {
    return -1;
  }

  for (int i = 0; i < this->width*this->height; i += 1)
  {
    if (this->framebuffer[i].red > max) max = this->framebuffer[i].red;
    if (this->framebuffer[i].green > max) max = this->framebuffer[i].green;
    if (this->framebuffer[i].blue > max) max = this->framebuffer[i].blue;
    if (this->framebuffer[i].red < min) min = this->framebuffer[i].red;
    if (this->framebuffer[i].green < min) min = this->framebuffer[i].green;
    if (this->framebuffer[i].blue < min) min = this->framebuffer[i].blue;
  }

  float diff = max - min;
  if (diff == 0.0f) diff = 1.0f;

  outfile << "P6\n";
  outfile << this->width << " " << this->height << "\n255\n";

  std::cout << "min/max/diff: " << min << "/" << max << "/" << diff << std::endl;

  
  for (int j = 0; j<  this->width*this->height; j += 1)
  {
    unsigned char pixel[3];
    pixel[0] = static_cast<unsigned char>(((framebuffer[j].red   - min) / diff) * 255.0f);
    pixel[1] = static_cast<unsigned char>(((framebuffer[j].green - min) / diff) * 255.0f);
    pixel[2] = static_cast<unsigned char>(((framebuffer[j].blue  - min) / diff) * 255.0f);

    outfile.write(reinterpret_cast<char*>(pixel), 3);

  /*  
    if (framebuffer[j].red !=0 || framebuffer[j].green != 0 || framebuffer[j].blue != 0) {
      std::cout << "j: " << j << std::endl;
      std::cout << "r: " << framebuffer[j].red << std::endl;
      std::cout << "g: " << framebuffer[j].green << std::endl;
      std::cout << "b: " << framebuffer[j].blue << std::endl;
    }
    
    outfile << (unsigned char)(((this->framebuffer[j].red-min)/diff)*255.0);
    outfile << (unsigned char)(((this->framebuffer[j].green-min)/diff)*255.0);
    outfile << (unsigned char)(((this->framebuffer[j].blue-min)/diff)*255.0);
    */  
  }
  

  
  
  outfile.close();
  return 0;
}

int FrameBuffer::writeDepthFile(char *filename)
{
  float max = 0;
  float min = FLT_MAX;

  ofstream outfile;

  outfile.open(filename, ofstream::binary);

  if (!outfile.is_open())
  {
    return -1;
  }

  for (int i = 0; i < this->width*this->height; i += 1)
  {
    if (this->framebuffer[i].depth > max) max = this->framebuffer[i].depth;
    if (this->framebuffer[i].depth < min) min = this->framebuffer[i].depth;
  }

  float diff = max - min;
  if (diff == 0.0f) diff = 1.0f;

  //  cerr << "Min/max/diff" << min << "/" << max << "/" << diff << endl;

  outfile << "P6\n";
  outfile << this->width << " " << this->height << "\n255\n";

  for (int j = 0; j<  this->width*this->height; j += 1)
  {
    unsigned char pd = (unsigned char)(((this->framebuffer[j].depth-min)/diff)*255.0);
    outfile << pd << pd << pd;
  }
  
  outfile.close();
  return 0;
}
