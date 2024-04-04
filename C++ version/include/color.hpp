/*
 *  color.hpp
 *  Header file for color: a summary class for different color spaces
 */

#pragma once
#include <color-util/RGB_to_XYZ.hpp>
#include <color-util/XYZ_to_Lab.hpp>
#include <color-util/XYZ_to_RGB.hpp>
#include <color-util/Lab_to_XYZ.hpp>
#include <color-util/RGB_to_HSL.hpp>
#include <color-util/CIEDE2000.hpp>
#include <color-util/type.hpp>
#include <iostream>
#include <vector>
#include <cmath>
#include "utils.hpp"

typedef colorutil::RGB rgb;
typedef colorutil::XYZ xyz;
typedef colorutil::HSL hsl;
typedef colorutil::Lab lab;

const double epsilon = 1.0;

inline void normrgb(rgb &rgb_color)
{
    rgb_color[0] = std::max(0.0, std::min(1.0, rgb_color[0]));
    rgb_color[1] = std::max(0.0, std::min(1.0, rgb_color[1]));
    rgb_color[2] = std::max(0.0, std::min(1.0, rgb_color[2]));
}

class Color
{
public:
    rgb rgb_color;
    lab lab_color;
    hsl hsl_color;
    double hue, chroma;

private:
    void norm()
    {
        normrgb(this->rgb_color);
        this->lab_color = colorutil::convert_XYZ_to_Lab(colorutil::convert_RGB_to_XYZ(rgb_color));
        double a = this->lab_color[1], b = this->lab_color[2];
        double hue = atan2(b, a) * 180 / M_PI;
        if (hue < 0)
            hue += 360;
        this->hue = hue;
        this->chroma = sqrt(a * a + b * b);
        this->hsl_color = colorutil::convert_RGB_to_HSL(rgb_color);
    }

public:
    Color() : rgb_color(0, 0, 0), lab_color(0, 0, 0), hue(0), chroma(0) {}

    Color(const Color &color) : rgb_color(color.rgb_color), lab_color(color.lab_color), hsl_color(color.hsl_color), hue(color.hue), chroma(color.chroma) {}

    Color(const rgb &rgb_color) : rgb_color(rgb_color), hsl_color(colorutil::convert_RGB_to_HSL(rgb_color)), lab_color(colorutil::convert_XYZ_to_Lab(colorutil::convert_RGB_to_XYZ(rgb_color)))
    {
        double a = this->lab_color[1], b = this->lab_color[2];
        double hue = atan2(b, a) * 180 / M_PI;
        if (hue < 0)
            hue += 360;
        this->hue = hue;
        this->chroma = sqrt(a * a + b * b);
    }

    Color(const lab &lab_color) : lab_color(lab_color)
    {
        this->rgb_color = colorutil::convert_XYZ_to_RGB(colorutil::convert_Lab_to_XYZ(lab_color));
        norm();
    }

    Color(double hue, double chroma, double luminance) : hue(hue), chroma(chroma)
    {
        double hue_rad = hue * M_PI / 180;
        this->lab_color = lab(luminance, chroma * cos(hue_rad), chroma * sin(hue_rad));
        this->rgb_color = colorutil::convert_XYZ_to_RGB(colorutil::convert_Lab_to_XYZ(lab_color));
        norm();
    }

    double dist(const Color &color) const
    {
        return colorutil::calculate_CIEDE2000(this->lab_color, color.lab_color);
    }

    void modifyInRange(const Color &base_color, double range, double hue_min = 0.0, double hue_max = 0.0)
    {
        if (range < 1e-5)
        {
            this->lab_color = lab(base_color.lab_color[0], base_color.lab_color[1], base_color.lab_color[2]);
            this->rgb_color = colorutil::convert_XYZ_to_RGB(colorutil::convert_Lab_to_XYZ(lab_color));
            this->norm();
            return;
        }
        double dis = this->dist(base_color);
        double ratio = dis > 0 ? range / dis : 1;
        double hue = this->hue, chroma = this->chroma, lumi = this->lab_color[0];
        int ct = 0;
        while (dis > range)
        {
            ct += 1;
            hue = hue - (hue - base_color.hue) * ratio;
            chroma = chroma - (chroma - base_color.chroma) * ratio;
            lumi = lumi - (lumi - base_color.lab_color[0]) * ratio;
            double hue_rad = hue * M_PI / 180;
            this->lab_color = lab(lumi, chroma * cos(hue_rad), chroma * sin(hue_rad));
            dis = this->dist(base_color);
            if (ct > 20)
                break;
        }

        if (abs(hue_max - hue_min) > 1e-5)
        {
            bool in_range = false;
            double delta1, delta2;
            if (hue_max > hue_min)
            {
                if (hue >= hue_min && hue <= hue_max)
                    in_range = true;
                else if (hue < hue_min)
                {
                    delta1 = hue_min - hue;
                    delta2 = 360 - (hue_max - hue);
                }
                else
                {
                    delta2 = hue - hue_max;
                    delta1 = 360 - (hue - hue_min);
                }
            }
            else
            {
                if (hue >= hue_min || hue <= hue_max)
                    in_range = true;
                else
                {
                    delta1 = hue_min - hue;
                    delta2 = hue - hue_max;
                }
            }

            if (!in_range)
            {
                if (delta1 < delta2)
                    hue = hue_min;
                else
                    hue = hue_max;
            }

            double hue_rad = hue * M_PI / 180;
            this->lab_color = lab(lumi, chroma * cos(hue_rad), chroma * sin(hue_rad));
        }

        this->rgb_color = colorutil::convert_XYZ_to_RGB(colorutil::convert_Lab_to_XYZ(lab_color));
        norm();
    }

    void modifyInConstraint(const int &base_id, const std::vector<Color> &base_colors)
    {
        auto base_color = base_colors[base_id];
        double dist = this->dist(base_color);

        double min_dist = dist + epsilon;
        int min_id = base_id;
        for (int i = 0; i < base_colors.size(); i++)
        {
            if (i == base_id)
                continue;
            double dist2 = this->dist(base_colors[i]);
            if (dist2 < min_dist)
            {
                min_dist = dist2;
                min_id = i;
            }
        }
        if (min_id != base_id)
        {
            double hue = this->hue, chroma = this->chroma, lumi = this->lab_color[0];
            double ratio = (dist - min_dist + epsilon) / dist;
            if (ratio > 1.0)
            {
                ratio = 1.0;
            }

            hue = hue - (hue - base_color.hue) * ratio;
            chroma = chroma - (chroma - base_color.chroma) * ratio;
            lumi = lumi - (lumi - base_color.lab_color[0]) * ratio;
            double hue_rad = hue * M_PI / 180;
            this->lab_color = lab(lumi, chroma * cos(hue_rad), chroma * sin(hue_rad));
            this->rgb_color = colorutil::convert_XYZ_to_RGB(colorutil::convert_Lab_to_XYZ(lab_color));
            norm();
        }
    }
};

inline double dist(const Color &color1, const Color &color2)
{
    return colorutil::calculate_CIEDE2000(color1.lab_color, color2.lab_color);
}

typedef std::vector<Color> Palette;

inline double hsv_saturation(const Color &color)
{
    double l = color.hsl_color[2];
    double s = color.hsl_color[1];
    double v = l + s * min(l, 1 - l);
    return ((v > 0) ? (2 - 2 * l / v) : 0);
}

inline Palette modifyMinHueDelta(const Palette &palette, double hue_delta)
{
    // get all hue
    vector<double> hues;
    for (auto &color : palette)
    {
        hues.push_back(color.hue);
    }
    auto arg_sort = argsort(hues);
    sort(hues.begin(), hues.end());
    vector<int> resort(arg_sort.size(), 0);
    for (int i = 0; i < arg_sort.size(); i++)
    {
        resort[arg_sort[i]] = i;
    }

    vector<double> new_hues;
    new_hues.push_back(hues[0]);
    double delta0 = hues[0] + 360 - hues[hues.size() - 1];
    for (int i = 1; i < hues.size(); i++)
    {
        double delta = max(hue_delta, hues[i] - hues[i - 1]);
        double new_hue = new_hues[i - 1] + delta;
        new_hues.push_back(new_hue);
    }
    double delta = (hue_delta > delta0) ? 0.0 : hue_delta;
    double new_hue = hues[hues.size() - 1] + delta;
    double r = 1.0;
    if (new_hue - hues[0] > 360.0)
        r = 360.0 / (new_hue - hues[0]);
    if (r < 1.0)
    {
        for (auto &h : new_hues)
        {
            h *= r;
        }
    }

    // rotate
    double sum_delta = 0;
    for (int i = 0; i < new_hues.size(); i++)
    {
        sum_delta += (new_hues[i] - hues[i]);
    }
    double average_delta = sum_delta / new_hues.size();
    for (auto &h : new_hues)
    {
        h -= average_delta;
        h = fmod(h, 360.0);
    }

    // modify
    Palette new_palette;
    for (int i = 0; i < palette.size(); i++)
    {
        auto &color = palette[i];
        new_palette.push_back(Color(new_hues[resort[i]], color.chroma, color.lab_color[0]));
    }
    return new_palette;
}