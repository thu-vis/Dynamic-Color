#include "palette_evaluator.hpp"

#include "evaluator_discri.hpp"
#include "evaluator_harmony.hpp"
#include "colorpair.hpp"

#define USE_OMP true

#if USE_OMP
#include "omp.h"

#define _TO_STR(x) #x
#define TO_STR(x) _TO_STR(x)

#define OMP_OP(x) omp x

#define PRAGMA_OMP(x) _Pragma(TO_STR(OMP_OP(x)))

#else
#define PRAGMA_OMP(x)
#endif

static pair<double, double> eval_perceptual_dist(const Palette &palette)
{
    double sum_dist = 0;
    double min_dist = numeric_limits<double>::max();

    PRAGMA_OMP(parallel for reduction(+ : sum_dist) reduction(min : min_dist))
    for (size_t i = 0; i < palette.size() - 1; i++)
    {
        double sum_dist_ = 0;
        double min_dist_ = numeric_limits<double>::max();

        PRAGMA_OMP(parallel for reduction(+ : sum_dist_) reduction(min : min_dist_))
        for (size_t j = i + 1; j < palette.size(); j++)
        {
            double dis = dist(palette[i], palette[j]);
            sum_dist_ += dis;
            min_dist_ = min(min_dist_, dis);
        }

        sum_dist += sum_dist_;
        min_dist = min(min_dist, min_dist_);
    }

    sum_dist /= (palette.size() * (palette.size() - 1) / 2);
    return {sum_dist, min_dist};
}

static pair<double, double> eval_name_dist(const Palette &palette)
{
    vector<int> color_indexes;
    vector<unsigned int> indices = {0, 0};
    vector<double> kddists = {0.0, 0.0};
    auto kdtree = get_kdtree();

    for (const auto &color : palette)
    {
        vector<double> lab = {color.lab_color[0], color.lab_color[1], color.lab_color[2]};
        kdtree->knnSearch(&lab[0], 1, &indices[0], &kddists[0]);
        color_indexes.push_back(indices[0]);
    }

    const auto &cosine_diff = get_cosine_diff();
    double sum_dist = 0;
    double min_dist = numeric_limits<double>::max();

    PRAGMA_OMP(parallel for reduction(+ : sum_dist) reduction(min : min_dist))
    for (size_t i = 0; i < palette.size() - 1; i++)
    {
        double sum_dist_ = 0;
        double min_dist_ = numeric_limits<double>::max();

        PRAGMA_OMP(parallel for reduction(+ : sum_dist_) reduction(min : min_dist_))
        for (size_t j = i + 1; j < palette.size(); j++)
        {
            double dis = cosine_diff[color_indexes[i]][color_indexes[j]];
            sum_dist_ += dis;
            min_dist_ = min(min_dist_, dis);
        }

        sum_dist += sum_dist_;
        min_dist = min(min_dist, min_dist_);
    }

    sum_dist /= (palette.size() * (palette.size() - 1) / 2);
    return {sum_dist, min_dist};
}

static pair<double, double> eval_palette_harmony(const Palette &palette)
{
    static const HarmonyEvaluator evaluator;
    return {evaluator.evaluateTemplates(palette), evaluator.evaluateLc(palette)};
}

static double eval_pair_harmony(const Palette &palette)
{
    double sum_harmo = 0;

    PRAGMA_OMP(parallel for reduction(+ : sum_harmo))
    for (size_t i = 0; i < palette.size() - 1; i++)
    {
        double sum_harmo_ = 0;

        PRAGMA_OMP(parallel for reduction(+ : sum_harmo_))
        for (size_t j = i + 1; j < palette.size(); j++)
        {
            double dis = pairColorHarmony(palette[i], palette[j]);
            sum_harmo_ += dis;
        }

        sum_harmo += sum_harmo_;
    }

    sum_harmo /= (palette.size() * (palette.size() - 1) / 2);
    return sum_harmo;
}

map<string, double> evaluate_palette(const vector<vector<double>> &in_palette)
{
    if (in_palette.size() <= 1)
    {
        throw invalid_argument("Too few colors in palette");
    }

    Palette palette;
    for (const auto &color : in_palette)
    {
        if (color.size() != 3)
        {
            throw invalid_argument("Not a valid palette");
        }

        palette.push_back(Color(rgb(color[0], color[1], color[2])));
    }

    auto [avg_per_dis, min_per_dis] = eval_perceptual_dist(palette);
    auto [avg_name_dis, min_name_dis] = eval_name_dist(palette);
    auto [tp, lc] = eval_palette_harmony(palette);
    double avg_pair_harmo = eval_pair_harmony(palette);

    return {
        {"avg_per_dis", avg_per_dis},
        {"min_per_dis", min_per_dis},
        {"avg_name_dis", avg_name_dis},
        {"min_name_dis", min_name_dis},
        {"templates", tp},
        {"lc", lc},
        {"avg_pair_harmo", avg_pair_harmo}};
}

map<string, double> evaluate_palette_grid(const vector<vector<double>> &in_palette,
                                          const int width, const int height, const vector<int> &data, const vector<int> &labels, const vector<int> &label_list, vector<vector<double>> &similarity)
{
    if (in_palette.size() <= 1)
    {
        throw invalid_argument("Too few colors in palette");
    }

    Palette palette;
    for (const auto &color : in_palette)
    {
        if (color.size() != 3)
        {
            throw invalid_argument("Not a valid palette");
        }

        palette.push_back(Color(rgb(color[0], color[1], color[2])));
    }

    GridData totaldata(Data::DrivenMode::discrimination);
    totaldata.initGridData(width, height, data, labels, label_list, similarity, 0, 1);

    auto [avg_per_dis, min_per_dis] = eval_perceptual_dist(palette);
    auto [avg_name_dis, min_name_dis] = eval_name_dist(palette);
    auto [tp, lc] = eval_palette_harmony(palette);
    double avg_pair_harmo = eval_pair_harmony(palette);
    double spatial_pair_harmony = totaldata.getDataDrivenValue(palette, 0, palette.size(), 0.0, 1.0);

    return {
        {"avg_per_dis", avg_per_dis},
        {"min_per_dis", min_per_dis},
        {"avg_name_dis", avg_name_dis},
        {"min_name_dis", min_name_dis},
        {"templates", tp},
        {"lc", lc},
        {"avg_pair_harmo", avg_pair_harmo},
        {"spatial_pair_harmo", spatial_pair_harmony}};
}