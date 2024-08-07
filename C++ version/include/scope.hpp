/*
 *  scope.hpp
 *  Header file for color scope and function for disturb color in scope
 */

#pragma once
#include "color.hpp"
#include "utils.hpp"

class Scope
{
public:
    static constexpr double hue_scope_default[2] = {0, 360};
    static constexpr double chroma_scope_default[2] = {40, 85};
    static constexpr double lumi_scope_default[2] = {40, 85}; 
    // If you intend to use it in a PDF document or on a darker scree, we recommend a higher lightness range, such as [55, 85]
    // If you intend to use it for data visualization on a brighter screen, we recommend a more looser range, such as [40, 85] which is the same with our paper.

    double hue_scope[2];
    double chroma_scope[2];
    double lumi_scope[2];
    bool all_hue = false;
    double margin = 10.0; // add margin to avoid overflow in disturb
    // we often found color overflow range when in boundary (when hcl to rgb)
    // so we add margin which only work when disturb color

    Scope() noexcept
    {
        std::copy(std::cbegin(hue_scope_default), std::cend(hue_scope_default), std::begin(hue_scope));
        std::copy(std::cbegin(chroma_scope_default), std::cend(chroma_scope_default), std::begin(chroma_scope));
        std::copy(std::cbegin(lumi_scope_default), std::cend(lumi_scope_default), std::begin(lumi_scope));
        if (hue_scope[0] == 0 && hue_scope[1] == 360)
        {
            all_hue = true;
        }
    }

    Scope(double *hue_scope, double *chroma_scope, double *lumi_scope) noexcept
    {
        std::copy(hue_scope, hue_scope + 2, std::begin(this->hue_scope));
        std::copy(chroma_scope, chroma_scope + 2, std::begin(this->chroma_scope));
        std::copy(lumi_scope, lumi_scope + 2, std::begin(this->lumi_scope));
        if (hue_scope[0] == 0 && hue_scope[1] == 360)
        {
            all_hue = true;
        }
    }

    double randomHue() const
    {
        double h;
        if (all_hue)
        {
            h = randomDouble(0, 360);
        }
        h = randomDouble(hue_scope[0], hue_scope[1]);
        if (h > 84 && h < 115)
        {
            if (h < 99.5)
                h = 84;
            else
                h = 115;
        }
        return h;
    }

    double randomChroma() const
    {
        return randomDouble(chroma_scope[0], chroma_scope[1]);
    }

    double randomLumi() const
    {
        return randomDouble(lumi_scope[0], lumi_scope[1]);
    }
};

inline Color disturbColorInScope(const Color &color, const Scope &scope, const double *range)
{
    // range: array of 3 for h c and l
    double h = color.hue, c = color.chroma, l = color.lab_color[0];
    c = randomDouble(std::max(c - range[1], scope.chroma_scope[0] + margin), std::min(c + range[1] + margin, scope.chroma_scope[1]));
    l = randomDouble(std::max(l - range[2], scope.lumi_scope[0] + margin), std::min(l + range[2] + margin, scope.lumi_scope[1]));
    if (!scope.all_hue)
    {
        h = randomDouble(std::max(h - range[0], scope.hue_scope[0]), std::min(h + range[0], scope.hue_scope[1]));
    }
    else
    {
        h = h + randomDoubleSym(range[0]);
        h = std::fmod(h, 360.0);
        if (h < 0)
            h += 360.0;
    }

    // remove square [85, 114]
    if (h > 84 - margin && h < 115 + margin)
    {
        if (h < 99.5)
            h = 84 - margin;
        else
            h = 115 + margin;
    }

    // ToDo: personal perference
    // such as lighter colors with green color
    // here is an example with little change, which will appear more probs
    if (h > 110 && h < 150) {
        c += 5;
        l += 5;
    }

    return Color(h, c, l);
}

inline Color modifyColorInScope(const Color &color, const Scope &scope)
{
    // range: array of 3 for h c and l
    double h = color.hue, c = color.chroma, l = color.lab_color[0];
    c = std::max(std::min(c, scope.chroma_scope[1]), scope.chroma_scope[0]);
    l = std::max(std::min(l, scope.lumi_scope[1]), scope.lumi_scope[0]);
    if (!scope.all_hue)
    {
        h = std::max(std::min(h, scope.hue_scope[1]), scope.hue_scope[0]);
    }
    else
    {
        h = std::fmod(h, 360.0);
        if (h < 0)
            h += 360.0;
    }

    // remove square [85, 114]
    if (h > 84 && h < 115)
    {
        if (h < 99.5)
            h = 84;
        else
            h = 115;
    }

    return Color(h, c, l);
}
