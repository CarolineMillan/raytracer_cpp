// Glass is a child class of Material

#pragma once

#include "material.h"

class Glass : public Material {
public:

    bool is_reflective() override { return true; };
    bool is_transparent() override { return true; };

    Glass(float eta = 1.5f, float kt = 0.9f, Colour tint = Colour(1.0f, 1.0f, 1.0f, 1.0f)) : eta(eta), kt(kt), tint(tint) {};

    float get_eta() override {return eta;};

    float get_kt() override {return kt;};

    Colour get_tint() override {
        return tint;
    };

private:
    float eta;
    float kt; // probability photon transmits in russian roulette
    Colour tint;
};
