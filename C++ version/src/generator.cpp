/*
 *  generator.cpp
 *  Implementation of color generator.
 */

#include "generator.hpp"
#include "mtl_adjustor.hpp"
#include <cassert>
#include <unordered_map>
#include <omp.h>
#include <memory>
using namespace std;

constexpr double Generator::color_min_range;
constexpr double Generator::color_max_range;
constexpr double Generator::hue_min_range;
constexpr double Generator::hue_max_range;
constexpr double Generator::global_init_dist;
constexpr double Generator::init_hue_min_delta;
constexpr double Generator::init_color_min_delta;

constexpr double Scope::chroma_scope_default[2];
constexpr double Scope::hue_scope_default[2];
constexpr double Scope::lumi_scope_default[2];

const double center_range[2] = {45, 80}; // check the environment
const double relax = 1.0;                // Set the strength of pareto specific preferences
const double harg = 1.0;                 // Harmony score weight (hue templates with l-c)

Palette Generator::findBasicColors(const Palette &basic_palette, vector<int> &palette_sizes, vector<shared_ptr<Data>> &driven_data)
{
    if (basic_palette.size() == 1)
    {
        gen_color = false;
        return basic_palette;
    }
    double l0 = scope.lumi_scope[0];
    double l1 = scope.lumi_scope[1];
    double c0 = scope.chroma_scope[0];
    double c1 = scope.chroma_scope[1];
    scope.lumi_scope[0] = center_range[0];
    scope.lumi_scope[1] = center_range[1];
    scope.chroma_scope[0] = center_range[0];
    scope.chroma_scope[1] = center_range[1];

    Palette res;
    if (basic_palette.size() == 0 && palette_sizes.size() > 0)
    {
        if (palette_sizes.size() == 1)
        {
            // if there is only one color, we can use a default color
            res.push_back(Color(rgb(0.996, 0.894, 0.702)));
        }
        else
        {
            Palette temp_palette;
            vector<int> temp_sizes;
            vector<shared_ptr<Data>> temp_data;
            if (driven_data.size() > 0)
            {
                auto cluster_data = driven_data[driven_data.size() - 1];
                temp_data.push_back(cluster_data);
                temp_data.push_back(cluster_data);
            }
            temp_palette.push_back(Color(scope.randomHue(), 0.5 * (l0 + l1), 0.5 * (c0 + c1)));
            temp_sizes.push_back(palette_sizes.size());
            res = findPaletteGlobal(temp_palette, temp_sizes, temp_data, true);
        }
        gen_color = true;
    }
    else
    {
        vector<int> temp_sizes(basic_palette.size(), 1);
        res = modifyMinHueDelta(basic_palette, init_hue_min_delta);
        res = findPaletteGlobal(res, temp_sizes, driven_data, true);
        gen_color = false;
    }

    // reset scope
    scope.lumi_scope[0] = l0;
    scope.lumi_scope[1] = l1;
    scope.chroma_scope[0] = c0;
    scope.chroma_scope[1] = c1;
    return res;
}

Palette Generator::findPaletteGlobal(const Palette &init_colors, vector<int> &palette_sizes, vector<shared_ptr<Data>> &driven_data, bool modify_pcolors)
{
    // ------------------------------------------- Init basic color info -------------------------
    Palette basic_colors;
    for (auto &color : init_colors)
    {
        basic_colors.push_back(modifyColorInScope(color, scope));
    }

    // set start constraint
    global_color_dist = global_init_dist;
    evaluator.global_dist = global_color_dist;
    this->modify_pcolors = modify_pcolors;
    // Previous literature often believed that hue was more important in harmony,
    // so we slightly adjusted some proportions. Of course, setting it to 1 is not a problem at all.
    double bharg = harg * 0.8;

    // init cluster info and estimate delta
    vector<ClusterInfo> clusters;
    vector<double> sqrt_nums;
    unordered_map<int, int> idx2id;
    vector<double> all_hues;
    int start = 0, end = 0;
    for (auto i = 0; i < basic_colors.size(); i++)
    {
        start = end;
        end += palette_sizes[i];
        clusters.push_back(ClusterInfo(i, start, end, palette_sizes[i], basic_colors[i]));
        all_hues.push_back(basic_colors[i].hue);
        sqrt_nums.push_back(sqrt(palette_sizes[i]));
        for (auto j = start; j < end; j++)
        {
            idx2id[j] = i;
        }
    }

    const double dist_ratio = tradeoff;
    if (all_hues.size() > 1)
    {
        auto arg_hues = argsort(all_hues);
        auto length = all_hues.size();
        int fi = length - 1;
        // cout << all_hues[arg_hues[fi]] << " " << all_hues[arg_hues[0]];
        for (auto i = 0; i < length; i++)
        {
            // by neighbor delta hue
            double min_range;
            if (i == 0)
            {
                double gap_hue = max(sqrt_nums[arg_hues[i]], sqrt_nums[arg_hues[fi]]) * dist_ratio;
                min_range = (360.0 - all_hues[arg_hues[fi]] + all_hues[arg_hues[i]]) * (sqrt_nums[arg_hues[i]] / (gap_hue + sqrt_nums[arg_hues[i]] + sqrt_nums[arg_hues[fi]]));
            }
            else
            {
                double gap_hue = max(sqrt_nums[arg_hues[i]], sqrt_nums[arg_hues[i - 1]]) * dist_ratio;
                min_range = (all_hues[arg_hues[i]] - all_hues[arg_hues[i - 1]]) * (sqrt_nums[arg_hues[i]] / (gap_hue + sqrt_nums[arg_hues[i]] + sqrt_nums[arg_hues[i - 1]]));
            }

            double max_range;
            if (i == fi)
            {
                double gap_hue = max(sqrt_nums[arg_hues[i]], sqrt_nums[arg_hues[0]]) * dist_ratio;
                max_range = (360.0 - all_hues[arg_hues[fi]] + all_hues[arg_hues[i]]) * (sqrt_nums[arg_hues[i]] / (gap_hue + sqrt_nums[arg_hues[i]] + sqrt_nums[arg_hues[0]]));
            }
            else
            {
                double gap_hue = max(sqrt_nums[arg_hues[i]], sqrt_nums[arg_hues[i + 1]]) * dist_ratio;
                max_range = (all_hues[arg_hues[i + 1]] - all_hues[arg_hues[i]]) * (sqrt_nums[arg_hues[i]] / (gap_hue + sqrt_nums[arg_hues[i]] + sqrt_nums[arg_hues[i + 1]]));
            }
            double eq_range = min(min_range, max_range);
            if (modify_pcolors)
                eq_range = max(eq_range, init_hue_min_delta);
            else
                eq_range = max(eq_range, hue_min_range);
            if (!modify_pcolors && clusters[arg_hues[i]].num == 1)
            {
                eq_range = 0.0;
            }

            clusters[arg_hues[i]].hue_min = fmod(all_hues[arg_hues[i]] - eq_range, 360.0);
            if (clusters[arg_hues[i]].hue_min < 0.0)
                clusters[arg_hues[i]].hue_min += 360.0;

            clusters[arg_hues[i]].hue_max = fmod(all_hues[arg_hues[i]] + eq_range, 360.0);
            if (clusters[arg_hues[i]].hue_max < 0.0)
                clusters[arg_hues[i]].hue_max += 360.0;

            clusters[arg_hues[i]].hue_range = eq_range;
        }
    }
    else
    {
        clusters[0].hue_min = 0.0;
        clusters[0].hue_max = 0.0;
        clusters[0].hue_range = 180;
        if ((!modify_pcolors && !gen_color))
            clusters[0].hue_range = min(estimate_delta * sqrt_nums[0], hue_max_range);
    }

    vector<vector<double>> dist_matrix(basic_colors.size(), vector<double>(basic_colors.size(), 0.0));
    double min_basic_dist = 30.0;
    for (auto i = 0; i < basic_colors.size(); i++)
    {
        for (auto j = i + 1; j < basic_colors.size(); j++)
        {
            auto ds = dist(basic_colors[i], basic_colors[j]);
            dist_matrix[i][j] = ds;
            dist_matrix[j][i] = ds;
            if (ds < min_basic_dist)
                min_basic_dist = ds;
        }
    }

    for (auto i = 0; i < basic_colors.size(); i++)
    {

        double min_range = 10000.0;
        for (auto j = 0; j < basic_colors.size(); j++)
        {
            if (i == j)
                continue;
            auto gap = max(sqrt_nums[i], sqrt_nums[j]) * dist_ratio;
            auto r = dist_matrix[i][j] * sqrt_nums[i] / (gap + sqrt_nums[i] + sqrt_nums[j]);
            if (r < min_range)
                min_range = r;
        }
        if (min_range > color_max_range)
        {
            if (!modify_pcolors && !gen_color)
            {
                // only 1 color to expand, estimate a smaller range for consistency
                min_range = min(estimate_delta * sqrt_nums[0], color_max_range);
            }
            else
            {
                min_range = color_max_range;
            }
        }

        if (modify_pcolors)
            min_range = max(min_range, init_color_min_delta);
        else
            min_range = max(min_range, color_min_range);

        if (!modify_pcolors && clusters[i].num == 1)
        {
            min_range = 0.0;
        }
        clusters[i].range = min_range;
        clusters[i].initrange = min_range;
    }

    shared_ptr<Data> total_data = nullptr;
    use_data = false;
    if (driven_data.size() > 0)
    {
        use_data = true;
        assert(driven_data.size() == (clusters.size() + 1));
        for (int i = 0; i < clusters.size(); i++)
        {
            clusters[i].data = driven_data[i];
        }
        driven_mode = driven_data[0]->mode;
        total_data = driven_data[driven_data.size() - 1];
    }
    evaluator.setClusters(clusters, total_data);

    // ------------------------------------------- Set Arguments ----------------------------------------
    double alpha = 0.0;

    const bool use_early_stop = true;
    double early_stop = 1e-6;
    int stop_iter_step0 = 50;
    int stop_iter = 10;
    int cur_stop_iter = 0;

    const double start_temper_step0 = 1e-3;
    const double end_temper_step0 = 1e-9;
    const double start_temper = 1e-3;
    const double end_temper = 1e-6;

    // ------------------------------------------- Formally starting optimization -------------------------
    // 1. first step: balance discri
    shift_p = 1.0;
    // evaluate start palette and start sa
    Palette palette;
    for (auto i = 0; i < clusters.size(); i++)
    {
        double range[3] = {clusters[i].hue_range, 60.0, 60.0};
        if (clusters[i].num == 1)
        {
            palette.push_back(Color(clusters[i].basic_color));
            continue;
        }
        for (auto j = 0; j < clusters[i].num; j++)
        {
            Color c = disturbColorInScope(clusters[i].basic_color, this->scope, range);
            if (modify_pcolors)
                c.modifyInConstraint(i, basic_colors);
            c.modifyInRange(clusters[i].basic_color, clusters[i].range, clusters[i].hue_min, clusters[i].hue_max);
            palette.push_back(c);
        }
    }

    double cur_min_dist = 0.0;
    double cur_score = evaluator.evaluate_discriminate(palette, cur_min_dist);
    double cur_template_score = evaluator.evaluate_templates(palette);

    Palette start_palette = palette;
    double start_score = cur_score;
    set_disturb_range(0, modify_pcolors);

    double dart_throw_r = 30.0;
    const double lr = 0.25;
    if (basic_colors.size() > 1)
        dart_throw_r = max(min(min_basic_dist, dart_throw_r), global_color_dist);

    double max_min_dist = 0.0;
    int rotate_n = 0;
    while (dart_throw_r >= 0.0)
    {
        bool find_satisfied = true;
#pragma omp parallel for
        for (int j = 0; j < iter_num; j++)
        {
            Palette tp_palette = start_palette;
            // Palette tp_palette = palette;
            int count = 0, sign = 0;
            while ((sign = evaluator.isDiscriminative(tp_palette, dart_throw_r)) > 0)
            {
                count += 1;
                if (count >= 100)
                {
                    find_satisfied = false;
                    break;
                }
                Color re_color = disturbColorInScope(tp_palette[sign], scope, disturb_step);
                int cluster_id = idx2id.at(sign);
                if (modify_pcolors)
                    re_color.modifyInConstraint(cluster_id, basic_colors);
                re_color.modifyInRange(clusters[cluster_id].basic_color, clusters[cluster_id].range, clusters[cluster_id].hue_min, clusters[cluster_id].hue_max);
                tp_palette[sign] = re_color;
            }
            double min_dist = 0.0;
            double new_score = evaluator.evaluate_discriminate(tp_palette, min_dist);

#pragma omp critical
            {
                if (new_score > start_score)
                {
                    start_score = new_score;
                    start_palette = tp_palette;
                }
                if (min_dist > max_min_dist)
                {
                    max_min_dist = min_dist;
                }
            }
        }

        if (find_satisfied || rotate_n > 50)
        {
            break;
        }

        dart_throw_r -= (dart_throw_r - max_min_dist) * lr;
        rotate_n += 1;
    }
    if (dart_throw_r < global_color_dist)
    {
        global_color_dist = dart_throw_r;
        stop_iter_step0 = 100;
    }

    evaluator.initDataArgs(start_palette);

    double best_score = start_score;
    Palette best_palette = start_palette;
    palette = start_palette;
    cur_score = start_score;
    set_disturb_range(1, modify_pcolors);

    double cur_temper = start_temper_step0;
    double dec = global_dec;
    int cur_iter = 0;

    double up_bound[2] = {best_score, 0.0};
    double low_bound[2] = {up_bound[0], cur_template_score}; // first save dis and template score
    double lc_bound[2] = {-1e-5, 0.0};                       // save lc score
    double data_bound[2] = {0.0, 0.0};                       // save data score

    while (cur_temper > end_temper_step0)
    {
        Palette temper_palette = palette;
        double temp_score = cur_score;
        double temp_best_score = best_score;
#pragma omp parallel for
        for (int j = 0; j < iter_num; j++)
        {
            Palette new_palette = disturbPalette(palette, clusters, idx2id, basic_colors);
            double min_dist = 0.0;
            double new_score = evaluator.evaluate_discriminate(new_palette, min_dist);
            double template_score = evaluator.evaluate_templates(new_palette);
            double lc_score = evaluator.evaluate_lc(new_palette);
            double data_score = 0.0;
            if (use_data)
                data_score = evaluator.evaluate_data(new_palette);
            double delta_score = cur_score - new_score;
#pragma omp critical
            {
                // estimate up bound and low bound
                if (new_score > best_score)
                {
                    best_score = new_score;
                    best_palette = new_palette;
                }

                if (delta_score <= 0 || randomOne() < exp((-delta_score) / cur_temper))
                {
                    temper_palette = new_palette;
                    temp_score = new_score;
                }

                if (min_dist > global_color_dist && new_score < low_bound[0])
                {
                    low_bound[0] = new_score;
                }

                if (template_score < low_bound[1])
                {
                    low_bound[1] = template_score;
                }

                if (lc_score < lc_bound[0])
                {
                    lc_bound[0] = lc_score;
                }

                if (data_score > data_bound[1])
                {
                    data_bound[1] = data_score;
                }

                if (data_score < data_bound[0])
                {
                    data_bound[0] = data_score;
                }
            }
        }

        palette = temper_palette;
        cur_score = temp_score;

        if (use_early_stop)
        {

            if (abs(temp_best_score - best_score) < early_stop)
            {
                cur_stop_iter += 1;
            }
            else
            {
                cur_stop_iter = 0;
            }
            if (cur_stop_iter > stop_iter_step0)
            {
                break;
            }
        }

        cur_iter += 1;
        cur_temper *= dec;
    }

    auto cur_items = evaluator.evaluate_items(best_palette);
    double cur_score_without_penalty = 2 * cur_items[0] + 0.1 * cur_items[1];
    up_bound[0] = cur_score_without_penalty;
    up_bound[1] = 0.0;
    low_bound[0] = low_bound[0] - 1e-5;
    low_bound[1] = low_bound[1] - 1e-5;
    double tp_bound[2] = {low_bound[1], up_bound[1]};
    double harmony_args[2] = {lc_bound[1] - lc_bound[0], bharg * (tp_bound[1] - tp_bound[0])};
    low_bound[1] = 1.0 * (tp_bound[0] * harmony_args[0] + lc_bound[0] * harmony_args[1]);

    // 2. second step: adjust palette hue harmony
    set_disturb_range(2, modify_pcolors);

    cur_temper = start_temper;
    palette = best_palette;
    double best_dis_score = evaluator.evaluate_discriminate(best_palette, cur_min_dist);
    double start_dis_score = best_dis_score;
    double best_harmony_score = evaluator.evaluate_templates(best_palette) * harmony_args[0] + evaluator.evaluate_lc(best_palette) * harmony_args[1];
    double cur_score1 = best_dis_score;
    double cur_score2 = best_harmony_score;

    cur_iter = 0;
    int batch_num = 5;
    double mt = 0.0;
    double min_value = 0.1;
    double args[2] = {1.0, 0.0};
    cur_stop_iter = 0;

    while (cur_temper > end_temper)
    {
        // MTL change Args
        if (cur_iter % batch_num == 0)
        {
            double norm_score1 = (best_dis_score - low_bound[0]) / (up_bound[0] - low_bound[0]);
            double norm_score2 = (best_harmony_score - low_bound[1]) / (up_bound[1] - low_bound[1]);
            double arg0 = 1.0 - norm_score1;
            double arg1 = 1.0 - norm_score2;
            if (arg0 < min_value)
            {
                arg0 = min_value;
            }
            if (arg1 < min_value)
            {
                arg1 = min_value;
            }
            double sum = arg0 + arg1;
            arg0 /= sum;
            arg1 /= sum;

            args[0] = mt * args[0] + (1 - mt) * arg0;
            args[1] = mt * args[1] + (1 - mt) * arg1;

            palette = best_palette;
            best_score = args[0] * norm_score1 + args[1] * norm_score2;
            cur_score = best_score;
        }

        Palette temper_palette = palette;
        double temper_score = cur_score;
        double temper_best_score = best_score;

        double temper_dis_score = cur_score1;
        double temper_harmony_score = cur_score2;
        double norm_start_dis_score = (start_dis_score - low_bound[0]) / (up_bound[0] - low_bound[0]);
        double norm_temp_dis_score = (temper_dis_score - low_bound[0]) / (up_bound[0] - low_bound[0]);
        double norm_temp_harmony_score = (temper_harmony_score - low_bound[1]) / (up_bound[1] - low_bound[1]);

#pragma omp parallel for
        for (int j = 0; j < iter_num; j++)
        {
            Palette new_palette = disturbPalette(palette, clusters, idx2id, basic_colors);
            double min_dist = 0.0;
            double dis_score = evaluator.evaluate_discriminate(new_palette, min_dist);
            double harmony_score = evaluator.evaluate_templates(new_palette) * harmony_args[0] + evaluator.evaluate_lc(new_palette) * harmony_args[1];
            double norm_dis_score = (dis_score - low_bound[0]) / (up_bound[0] - low_bound[0]);
            double norm_harmony_score = (harmony_score - low_bound[1]) / (up_bound[1] - low_bound[1]);

            double new_score = args[0] * norm_dis_score + args[1] * norm_harmony_score;
            double delta_score = cur_score - new_score;

            bool importance_factor = (norm_temp_dis_score > norm_temp_harmony_score * relax)
                                         ? (norm_dis_score > norm_harmony_score * relax)
                                         : (norm_dis_score >= norm_start_dis_score);
#pragma omp critical
            {
                if (importance_factor)
                {
                    if (new_score > best_score)
                    {
                        best_score = new_score;
                        best_palette = new_palette;
                        best_dis_score = dis_score;
                        best_harmony_score = harmony_score;
                    }

                    if (delta_score <= 0 || randomOne() < exp((-delta_score) / cur_temper))
                    {
                        temper_palette = new_palette;
                        temper_score = new_score;
                        temper_dis_score = dis_score;
                        temper_harmony_score = harmony_score;
                    }
                }

                if (dis_score > up_bound[0])
                {
                    up_bound[0] = dis_score;
                }
            }
        }

        palette = temper_palette;
        cur_score = temper_score;
        cur_score1 = temper_dis_score;
        cur_score2 = temper_harmony_score;

        cur_iter += 1;
        cur_temper *= dec;

        if (use_early_stop)
        {

            if (abs(temper_best_score - best_score) < early_stop)
            {
                cur_stop_iter += 1;
            }
            else
            {
                cur_stop_iter = 0;
            }
            if (cur_stop_iter > stop_iter)
            {
                break;
            }
        }
    }

    best_value_save.clear();
    if (!modify_pcolors)
    {
        best_value_save.push_back(tp_bound[0]);
        best_value_save.push_back(tp_bound[1]);
        best_value_save.push_back(lc_bound[0]);
        best_value_save.push_back(lc_bound[1]);
    }
    if (!use_data)
        return best_palette;

    // 3. third step: data driven score
    set_disturb_range(3, modify_pcolors);
    alpha = args[1] / args[0];

    cur_temper = start_temper;
    palette = best_palette;

    start_dis_score = evaluator.evaluate_discriminate(palette, cur_min_dist);
    double start_harmony_score = evaluator.evaluate_templates(palette) * harmony_args[0] + evaluator.evaluate_lc(palette) * harmony_args[1];
    double best_score_1 = start_dis_score + alpha * start_harmony_score;
    double best_score_2 = evaluator.evaluate_data(best_palette);

    cur_stop_iter = 0;
    cur_iter = 0;
    shift_p = 0.1;

    double dis_bound[2] = {low_bound[0], up_bound[0]};
    double harmony_bound[2] = {low_bound[1], up_bound[1]};

    args[0] = 1.0;
    args[1] = 0.0;
    up_bound[0] = best_score_1;
    low_bound[0] += low_bound[1] * alpha;
    up_bound[1] = data_bound[1];
    low_bound[1] = data_bound[0];

    while (cur_temper > end_temper)
    {
        // TODO: MTL change Args
        if (cur_iter % batch_num == 0)
        {
            double norm_best_score1 = (best_score_1 - low_bound[0]) / (up_bound[0] - low_bound[0]);
            double norm_best_score2 = ((best_score_2 - low_bound[1]) / (up_bound[1] - low_bound[1]));
            double arg0 = 1.0 - norm_best_score1;
            double arg1 = 1.0 - norm_best_score2;
            if (arg0 < min_value)
            {
                arg0 = min_value;
            }
            if (arg1 < min_value)
            {
                arg1 = min_value;
            }
            double sum = arg0 + arg1;
            arg0 /= sum;
            arg1 /= sum;

            args[0] = mt * args[0] + (1 - mt) * arg0;
            args[1] = mt * args[1] + (1 - mt) * arg1;

            palette = best_palette;

            best_score = args[0] * norm_best_score1 + args[1] * norm_best_score2;
            cur_score = best_score;
        }
        Palette temper_palette = palette;
        double temp_score = cur_score;
        double temp_best_score = best_score;
        double temp_dis_score = evaluator.evaluate_discriminate(palette, cur_min_dist);
        double temp_harmony_score = evaluator.evaluate_templates(palette) * harmony_args[0] + evaluator.evaluate_lc(palette) * harmony_args[1];
        double temp_data_score = evaluator.evaluate_data(palette);

        double norm_temp_dis_score = (temp_dis_score - dis_bound[0]) / (dis_bound[1] - dis_bound[0]);
        double norm_temp_harmony_score = (temp_harmony_score - harmony_bound[0]) / (harmony_bound[1] - harmony_bound[0]);
        double norm_start_dis_score = (start_dis_score - dis_bound[0]) / (dis_bound[1] - dis_bound[0]);
        double norm_start_harmony_score = (start_harmony_score - harmony_bound[0]) / (harmony_bound[1] - harmony_bound[0]);
        double norm_temp_data_score = (temp_data_score - low_bound[1]) / (up_bound[1] - low_bound[1]);

#pragma omp parallel for
        for (int j = 0; j < iter_num; j++)
        {
            Palette new_palette = disturbPalette(palette, clusters, idx2id, basic_colors);
            double min_dist = 0.0;
            double dis_score = evaluator.evaluate_discriminate(new_palette, min_dist);
            double harmony_score = evaluator.evaluate_templates(new_palette) * harmony_args[0] + evaluator.evaluate_lc(new_palette) * harmony_args[1];
            double score1 = dis_score + alpha * harmony_score;
            double score2 = evaluator.evaluate_data(new_palette);

            double norm_score1 = (score1 - low_bound[0]) / (up_bound[0] - low_bound[0]);
            double norm_score2 = (score2 - low_bound[1]) / (up_bound[1] - low_bound[1]);
            double norm_dis = (dis_score - dis_bound[0]) / (dis_bound[1] - dis_bound[0]);
            double norm_harmony = (harmony_score - harmony_bound[0]) / (harmony_bound[1] - harmony_bound[0]);

            double new_score = args[0] * norm_score1 + args[1] * norm_score2;
            double delta_score = cur_score - new_score;

            bool importance_factor =
                (norm_temp_dis_score > norm_temp_harmony_score * relax)
                    ? (norm_dis > norm_harmony * relax) || ((norm_temp_harmony_score > norm_temp_data_score * relax) ? (norm_harmony > norm_score2 * relax) : (norm_harmony >= norm_start_harmony_score))
                    : (norm_dis >= norm_start_dis_score);

#pragma omp critical
            {
                if (importance_factor)
                {
                    if (new_score > best_score)
                    {
                        best_score = new_score;
                        best_palette = new_palette;
                        best_score_1 = score1;
                        best_score_2 = score2;
                    }

                    if (delta_score <= 0 || randomOne() < exp((-delta_score) / cur_temper))
                    {
                        temper_palette = new_palette;
                        temp_score = new_score;
                    }
                }
                if (score2 > up_bound[1])
                {
                    up_bound[1] = score2;
                }
            }
        }

        palette = temper_palette;
        cur_score = temp_score;

        cur_iter += 1;
        cur_temper *= dec;
    }
    return best_palette;
}

Palette Generator::disturbPalette(const Palette &palette, const vector<ClusterInfo> &clusters, const unordered_map<int, int> &idx2id, const vector<Color> &basic_colors)
{
    if (randomOne() < shift_p)
    {
        return disturbValue(palette, clusters, idx2id, basic_colors);
    }
    else
    {
        if (driven_mode == Data::DrivenMode::discrimination)
        {
            return disturbPosition(palette, clusters, idx2id);
        }
        else
        {
            static int iter_random = 5;
            Palette p = disturnPositionWheel(palette, clusters);
            double v = clusterDataDrivenValue(p, clusters);
            for (int i = 0; i < iter_random; i++)
            {
                Palette p1 = disturbPosition(p, clusters, idx2id);
                double v1 = clusterDataDrivenValue(p1, clusters);
                if (v1 > v)
                {
                    p = p1;
                    v = v1;
                }
            }
            if (v < clusterDataDrivenValue(palette, clusters))
            {
                return disturbPosition(palette, clusters, idx2id);
            }
            else
            {
                return p;
            }
        }
    }
}

Palette Generator::disturbValue(const Palette &palette, const vector<ClusterInfo> &clusters, const unordered_map<int, int> &idx2id, const vector<Color> &basic_colors)
{
    Palette new_palette = palette;
    int idx = randomInt(0, new_palette.size() - 1);
    Color dis_color = disturbColorInScope(new_palette[idx], scope, disturb_step);
    int cluster_id = idx2id.at(idx);
    if (modify_pcolors)
        dis_color.modifyInConstraint(cluster_id, basic_colors);
    dis_color.modifyInRange(clusters[cluster_id].basic_color, clusters[cluster_id].range, clusters[cluster_id].hue_min, clusters[cluster_id].hue_max);
    new_palette[idx] = dis_color;
    return new_palette;
}

Palette Generator::disturbPosition(const Palette &palette, const vector<ClusterInfo> &clusters, const unordered_map<int, int> &idx2id)
{
    Palette new_palette = palette;
    int idx1 = randomInt(0, new_palette.size() - 1);
    int swap_part = idx2id.at(idx1);
    if (clusters[swap_part].num == 1)
        return new_palette;
    int idx2 = randomInt(clusters[swap_part].start, clusters[swap_part].end - 1);
    int ct = 0;
    while (idx1 == idx2)
    {
        idx2 = randomInt(clusters[swap_part].start, clusters[swap_part].end - 1);
        ct += 1;
        if (ct > 5)
            return new_palette;
    }
    swap(new_palette[idx1], new_palette[idx2]);
    return new_palette;
}

Palette Generator::disturnPositionWheel(const Palette &palette, const vector<ClusterInfo> &clusters)
{
    Palette new_palette = palette;
    for (auto &cluster : clusters)
    {
        if (cluster.num == 1)
            continue;
        auto wheel = cluster.data->getLabelSequence();
        vector<double> huelist;
        for (int i = cluster.start; i < cluster.end; i++)
        {
            huelist.push_back(palette[i].hsl_color[0]);
        }
        auto arghue = argsort(huelist);
        for (int i = 0; i < cluster.num; i++)
        {
            new_palette[cluster.start + wheel[i]] = palette[cluster.start + arghue[i]];
        }
    }
    return new_palette;
}