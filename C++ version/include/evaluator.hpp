/*
 *  evaluator.hpp
 *  Header file for color evaluation
 */

#pragma once
#include "color.hpp"
#include "utils.hpp"
#include "data.hpp"
#include "evaluator_discri.hpp"
#include "evaluator_harmony.hpp"
#include <omp.h>

#include "colorpair.hpp"

class PaletteEvaluator
{
public:
    double A = 1.0;
    const DiscriEvaluator discri_evaluator;
    const HarmonyEvaluator harmony_evaluator;

    array<double, 2> data_args = {1.0, 1.0}; // data-driven-first, data-driven-harmony: default we use normalization
    shared_ptr<Data> total_data = nullptr;
    std::vector<ClusterInfo> *clusters = nullptr;
    bool use_data = false;
    double global_dist = 10.0;

    PaletteEvaluator() noexcept {}
    ~PaletteEvaluator() noexcept {}

    void initDataArgs(const Palette &init_palette)
    {
        if (use_data)
        {
            auto data_score = total_data->getDataDrivenValue(init_palette, 0, init_palette.size(), 1.0, 0.0);
            auto pair_harmony_score = total_data->getDataDrivenValue(init_palette, 0, init_palette.size(), 0.0, 1.0);
            // use start two score to normlize these two items
            data_args[0] = min(fabs(floor2(pair_harmony_score)), 1.0);
            data_args[1] = min(fabs(floor2(data_score)), 1.0);
            if (total_data->mode == Data::DrivenMode::similarity)
            {
                // performance is better when similarity with larger args
                // but of course, it's also ok to set it to 1.0
                data_args[1] *= 3.0;
                // data_args[1] *= 1.0;
            }
        }
    }

    void setClusters(std::vector<ClusterInfo> &cluster_info, shared_ptr<Data> data = nullptr)
    {
        clusters = &cluster_info;
        total_data = data;
        if (data != nullptr)
            use_data = true;
    }

    int items_num = 5;
    vector<double> evaluate_items(const Palette &palette) const
    {
        vector<double> items;
        items.push_back(discri_evaluator.getNameDifference(palette));
        items.push_back(discri_evaluator.getValueDifference(palette));
        items.push_back(harmony_evaluator.evaluateTemplates(palette));
        items.push_back(harmony_evaluator.evaluateLc(palette));
        if (use_data)
            items.push_back(total_data->getDataDrivenValue(palette, 0, palette.size(), data_args[0], data_args[1]));
        return items;
    }

    double evaluate_discriminate(const Palette &palette, double &dist) const
    {
        dist = discri_evaluator.getValueDifference(palette);
        return discri_evaluator.getNameDifference(palette) * 2 + dist * 0.1 + A * min(dist - global_dist, 0.0);
    }

    double evaluate_templates(const Palette &palette) const
    {
        return harmony_evaluator.evaluateTemplates(palette);
    }

    double evaluate_lc(const Palette &palette) const
    {
        return harmony_evaluator.evaluateLc(palette);
    }

    double evaluate_data(const Palette &palette) const
    {
        return total_data->getDataDrivenValue(palette, 0, palette.size(), data_args[0], data_args[1]);
    }

    int isDiscriminative(const Palette &palette, const double &global_dist) const noexcept
    {
        return discri_evaluator.isDiscriminative(palette, global_dist);
    }
};