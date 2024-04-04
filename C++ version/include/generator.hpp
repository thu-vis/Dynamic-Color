/*
 *  generator.hpp
 *  Header file for palette and color generator
 */

#pragma once
#include "scope.hpp"
#include "evaluator.hpp"
#include <array>
#include <fstream>
using namespace std;

class Generator
{
public:
    // 0. basic arguments
    bool gen_color = true;
    double estimate_delta = 10.0; // only 1 fixed parent color --> a smaller space for consistency
    // larger space for pcolors but use constraints in related distance
    static constexpr double init_color_min_delta = 20.0;
    static constexpr double init_hue_min_delta = 20.0;

    // 1. color arguments
    Scope scope;
    static constexpr double color_min_range = 10.0;
    static constexpr double hue_min_range = 5.0;
    static constexpr double color_max_range = 100.0;
    static constexpr double hue_max_range = 180.0;

    double disturb_step[3] = {60.0, 60.0, 60.0};
    double default_disturb_step[3] = {60.0, 60.0, 60.0};
    static constexpr double global_init_dist = 10.0;
    double global_color_dist = 10.0;
    bool modify_pcolors = false;
    // 2. SA arguments
    double shift_p = 1.0;
    double global_dec = 0.99;
    int iter_num = 48;
    double tradeoff = 1.0; // trade off betweeen in-class and inter-class
    // 3. for evaluation
    PaletteEvaluator evaluator = PaletteEvaluator();
    bool use_data = false;
    Data::DrivenMode driven_mode = Data::DrivenMode::discrimination;

    // for test
    vector<double> best_value_save;

    Generator() {}

    // before version
    Palette findBasicColors(const Palette &palette, vector<int> &palette_sizes, vector<shared_ptr<Data>> &driven_data);

    Palette findPaletteGlobal(const Palette &init_colors, vector<int> &palette_sizes, vector<shared_ptr<Data>> &driven_data, bool init = false);

    void outputPalette(const Palette &palette, string filepath)
    {
        ofstream out(filepath);
        for (auto &color : palette)
        {
            out << color.rgb_color[0] * 255 << " " << color.rgb_color[1] * 255 << " " << color.rgb_color[2] * 255 << std::endl;
        }
    }

private:
    Palette disturbPalette(const Palette &palette, const vector<ClusterInfo> &clusters, const unordered_map<int, int> &idx2id, const vector<Color> &basic_colors);
    Palette disturbValue(const Palette &palette, const vector<ClusterInfo> &clusters, const unordered_map<int, int> &idx2id, const vector<Color> &basic_colors);
    Palette disturbPosition(const Palette &palette, const vector<ClusterInfo> &clusters, const unordered_map<int, int> &idx2id);
    Palette disturnPositionWheel(const Palette &palette, const vector<ClusterInfo> &clusters);
    void set_disturb_range(int stage, bool modify_parent_color = false)
    {
        switch (stage)
        {
        case 0:
            // init start palette with blue noise
            disturb_step[0] = 60.0;
            disturb_step[1] = 60.0;
            disturb_step[2] = 60.0;
            break;
        case 1:
            // modify discriminability
            disturb_step[0] = 30.0;
            disturb_step[1] = 30.0;
            disturb_step[2] = 30.0;
            break;
        case 2:
            // modify harmony: finetune
            disturb_step[0] = 10.0;
            disturb_step[1] = 10.0;
            disturb_step[2] = 10.0;
            break;
        case 3:
        default:
            // modify data: change little
            disturb_step[0] = 2.0;
            disturb_step[1] = 2.0;
            disturb_step[2] = 2.0;
            break;
        }
    }
};
