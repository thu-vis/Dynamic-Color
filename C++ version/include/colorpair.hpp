/*
 *  colorpair.hpp
 *  Header file for color pair harmony judgements
 */

#pragma once
#include "color.hpp"
#include <cmath>

inline double pairColorHarmony(const Color &color1, const Color &color2)
{
    // implementaion of color harmony model in paper: Universal models of colour emotion and colour harmony 2018
    // return -0.7 * tanh(-0.7 + 0.04 * fabs(color1.hue - color2.hue)) - 0.3 * tanh(-1.1 + 0.05 * fabs(color1.chroma - color2.chroma)) +
    //        0.4 * tanh(-0.8 + 0.05 * fabs(color1.lab_color[0] - color2.lab_color[0])) +
    //        0.3 + 0.6 * tanh(-4.2 + 0.028 * (color1.lab_color[0] + color2.lab_color[0]));

    // norm to 0 - 1
    return 0.5 - 0.175 * tanh(-0.7 + 0.04 * fabs(color1.hue - color2.hue)) - 0.075 * tanh(-1.1 + 0.05 * fabs(color1.chroma - color2.chroma)) +
           0.1 * tanh(-0.8 + 0.05 * fabs(color1.lab_color[0] - color2.lab_color[0])) +
           +0.15 * tanh(-4.2 + 0.028 * (color1.lab_color[0] + color2.lab_color[0]));
}