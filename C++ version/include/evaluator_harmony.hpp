/*
 *  evaluator_harmony.hpp
 *  Header file for color evaluation of harmony
 */

#pragma once
#include "color.hpp"
#include "utils.hpp"
#include "cluster.hpp"
#include <unordered_map>
#include <unordered_set>
#include <array>
#include <vector>
#include <cmath>
#include <optional>
#include <algorithm>
using namespace std;

class HueTemplate
{
public:
    struct Bound
    {
        double val, del;
    };

private:
    using range_t = pair<double, double>;
    vector<range_t> ranges, gaps;
    vector<Bound> bounds;

public:
    HueTemplate() = default;

    HueTemplate(const vector<range_t> &raw_ranges) : ranges(raw_ranges)
    {
        gaps.reserve(ranges.size());
        bounds.reserve(ranges.size() * 3);

        for (size_t i = 0; i + 1 < ranges.size(); i++)
        {
            double l = ranges[i].second, r = ranges[i + 1].first;
            double mid = (l + r) / 2;
            gaps.push_back({l, r});
            bounds.push_back({l, 1});
            bounds.push_back({mid, -2});
            bounds.push_back({r, 1});
        }

        double last_l = ranges.back().second;
        double last_r = ranges.front().first + 360;
        double last_mid = (last_l + last_r) / 2;
        gaps.push_back({last_l, last_r});

        if (last_l >= 360)
        {
            last_l -= 360;
        }
        if (last_mid >= 360)
        {
            last_mid -= 360;
        }
        last_r -= 360;

        bounds.push_back({last_l, 1});
        bounds.push_back({last_mid, -2});
        bounds.push_back({last_r, 1});
    }

    double hue_min_dist(double hue) const noexcept
    {
        // find the first range whose left bound >= hue
        auto r = lower_bound(ranges.begin(), ranges.end(), range_t{0, hue},
                             [](const range_t &a, const range_t &b)
                             { return a.second < b.second; });

        // happens when hue <= the right bound of the last range <= 360
        if (r == ranges.end())
        {
            return min(hue - ranges.back().second, ranges.front().first + 360 - hue);
        }

        if (r->first <= hue)
        {
            return 0;
        }

        auto pre = r != ranges.begin() ? r - 1 : ranges.end() - 1;

        double pre_second = pre->second;
        // happens when r == ranges.begin()
        if (pre_second > hue)
        {
            pre_second -= 360;
        }

        // hue falls into the first range
        if (pre_second > hue)
        {
            return 0;
        }

        return min(r->first - hue, hue - pre_second);
    }

    pair<double, double> checkMinDist(const vector<double> &hues, const optional<vector<double>> &saturs) const noexcept
    {
        constexpr double deg2rad = M_PI / 180.0;

        vector<Bound> all_bounds;
        all_bounds.reserve(bounds.size() * hues.size());

        // the derivative of total dist at current rotation
        double deriv = 0, cur_val = 0;
        for (size_t i = 0; i < hues.size(); i++)
        {
            double p = hues[i];
            double w = saturs.has_value() ? saturs.value()[i] : 1;

            for (const auto &b : bounds)
            {
                double new_val = b.val - p;

                if (new_val == 0)
                {
                    continue;
                }
                if (new_val < 0)
                {
                    new_val += 360;
                }

                all_bounds.push_back({new_val, b.del * w});
            }

            if (p < ranges.front().first)
            {
                p += 360;
            }

            // By using upper_bound, we ensure that g->second > p.
            auto g = upper_bound(gaps.begin(), gaps.end(), range_t{0, p},
                                 [](const range_t &a, const range_t &b)
                                 { return a.second < b.second; });

            assert(g != gaps.end());
            if (g->first <= p)
            {
                double mid = (g->first + g->second) / 2;
                deriv += p < mid ? w : -w;
            }

            cur_val += hue_min_dist(hues[i]) * w;
        }

        sort(all_bounds.begin(), all_bounds.end(), [](const Bound &b1, const Bound &b2)
             { return b1.val < b2.val; });

        double cur_rotate = 0;

        double min_val = cur_val;
        double min_rotate = cur_rotate;

        // enumerate all rotation angles where deriv changes
        for (const auto &b : all_bounds)
        {
            cur_val += deriv * (b.val - cur_rotate);
            deriv += b.del;
            cur_rotate = b.val;

            if (cur_val < min_val)
            {
                min_val = cur_val;
                min_rotate = cur_rotate;
            }
        }

        if (min_rotate >= 360)
        {
            min_rotate -= 360;
        }

        return {min_val * deg2rad, min_rotate};
    }
};

class MatsudaTemplates
{
public:
    unordered_map<char, HueTemplate> alter_templates;
    vector<char> used_names = {'i', 'V', 'L', 'I', 'T', 'Y', 'X'};
    unordered_set<char> all_names = {'i', 'V', 'L', 'I', 'T', 'Y', 'X'};
    vector<HueTemplate> templates;
    MatsudaTemplates() noexcept
    {
        initTemplates();
    }
    MatsudaTemplates(const vector<char> &use_names)
    {
        initTemplates();
        used_names.clear();
        for (const auto &name : use_names)
        {
            if (all_names.find(name) == all_names.end())
            {
                throw invalid_argument("Invalid template name");
            }
            used_names.push_back(name);
            templates.push_back(alter_templates.at(name));
        }
    }

    void changeUsedTemplates(const vector<char> &use_names)
    {
        used_names.clear();
        templates.clear();
        for (const auto &name : use_names)
        {
            if (all_names.find(name) == all_names.end())
            {
                throw invalid_argument("Invalid template name");
            }
            used_names.push_back(name);
            templates.push_back(alter_templates.at(name));
        }
    }

    double checkTemplatesMin(const vector<double> &hues, const optional<vector<double>> &saturs) const noexcept
    {
        double min_dist = 360.0;
        for (int i = 0; i < used_names.size(); i++)
        {
            const auto res = templates[i].checkMinDist(hues, saturs);
            auto value = get<0>(res);
            if (get<0>(res) < min_dist)
            {
                min_dist = value;
            }
        }
        // cout << min_dist << endl;
        return -min_dist;
    }

private:
    void initTemplates() noexcept
    {
        // ranges need to be non-overlapping and sorted ascendingly
        alter_templates['i'] = HueTemplate({{360 - 9.0, 360 + 9.0}});
        alter_templates['V'] = HueTemplate({{360 - 46.8, 360 + 46.8}});
        alter_templates['L'] = HueTemplate({{50.4, 129.6}, {360 - 9.0, 360 + 9.0}});
        alter_templates['I'] = HueTemplate({{171.0, 189.0}, {360 - 9.0, 360 + 9.0}});
        alter_templates['T'] = HueTemplate({{0.0, 180.0}});
        alter_templates['Y'] = HueTemplate({{171.0, 189.0}, {360 - 46.8, 360 + 46.8}});
        alter_templates['X'] = HueTemplate({{133.2, 226.8}, {360 - 46.8, 360 + 46.8}});
    }
};

class GeometricLC
{
public:
    double kc = 1.0;
    double kl = 2.0;
    GeometricLC() noexcept {}
    double linearLoss(const Palette &palette) const noexcept
    {
        double _c = 0.0, _l = 0.0, _w = 0.0;
        vector<double> ws;
        for (const auto &color : palette)
        {
            double sc = 1.0 + 0.045 * color.chroma;
            double lpow = pow(color.lab_color[0] - 50, 2);
            double sl = 1.0 + 0.015 * lpow / sqrt(20 + lpow);
            double w = pow(kc * sc + kl * sl, -2);
            _c += w * color.chroma;
            _l += w * color.lab_color[0];
            _w += w;
            ws.push_back(w);
        }
        _c /= _w;
        _l /= _w;
        double tp1 = 0.0, tp2 = 0.0;
        for (int i = 0; i < palette.size(); i++)
        {
            tp1 += ws[i] * (_l - palette[i].lab_color[0]) * (_c - palette[i].chroma);
            tp2 += ws[i] * (pow(palette[i].lab_color[0] - _l, 2) - pow(palette[i].chroma - _c, 2));
        }
        double _phi = 0.5 * atan2(-2 * tp1, tp2);
        double cos_phi = cos(_phi), sin_phi = sin(_phi);
        double _r = _c * cos_phi + _l * sin_phi;

        double _MD = 0.0;
        for (int i = 0; i < palette.size(); i++)
        {
            double MD = fabs(palette[i].chroma * cos_phi + palette[i].lab_color[0] * sin_phi - _r);
            _MD += max(MD - 15.0, 0.0); // MD +
        }
        return _MD;
    }
};

class HarmonyEvaluator
{
public:
    MatsudaTemplates templates;
    GeometricLC geometric_lc;
    HarmonyEvaluator(vector<char> used_names = {'L', 'T', 'X'}) noexcept // {'i', 'V', 'L', 'I', 'T', 'Y', 'X'}
    {
        templates.changeUsedTemplates(used_names);
    };

    double evaluateTemplates(const Palette &palette, bool use_saturs = true) const noexcept
    {
        vector<double> hues, saturs;

        // cout << "hues and saturations" << endl;
        for (const auto &color : palette)
        {
            hues.push_back(fmod(color.hsl_color[0] * 360.0, 360.0));
            // saturations.push_back(hsv_saturation(color));
            // cout << color.hsl_color[0] * 360 << " " << hsv_saturation(color) << endl;
            // hues.push_back(color.hue);
        }

        if (use_saturs)
        {
            for (const auto &color : palette)
            {
                saturs.push_back(hsv_saturation(color));
            }

            return templates.checkTemplatesMin(hues, saturs);
        }

        return templates.checkTemplatesMin(hues, {});
    };

    double evaluateLc(const Palette &palette) const noexcept
    {
        return -geometric_lc.linearLoss(palette);
    };
    ~HarmonyEvaluator() noexcept {}
};
