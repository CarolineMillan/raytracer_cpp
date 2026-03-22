// Metal is a child class of Material

#pragma once

#include "material.h"

class Metal : public Material {
public:
    bool is_reflective() override { return true; };
    Metal(float kr = 0.95f, Colour colour = Colour(0.1f, 0.05f, 0.0f, 1.0f)) : kr(kr), colour(colour) {};

    float get_kr() override {return kr;};

    void compute_base_colour(Colour &result) override {
        result.r = colour.r;
        result.g = colour.g;
        result.b = colour.b;
    };
private:
    float kr;
    Colour colour;
};
