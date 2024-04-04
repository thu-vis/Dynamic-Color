/*
 *  evaluator_discri.hpp
 *  Header file for color evaluation of discrimination
 */

#pragma once
#include "color.hpp"
#include "utils.hpp"
#include "kdtree.hpp"
#include <unordered_map>
#include <array>
using namespace std;

void loadDifferenceInfo();
shared_ptr<KDTree> get_kdtree();
const vector<vector<float>> &get_cosine_diff();

class DiscriEvaluator
{
public:
    DiscriEvaluator();
    ~DiscriEvaluator();

    double evaluate(const Palette &palette, double alpha, double beta, double &rt_min_dist) const noexcept;
    double getNameDifference(const Palette &palette) const;
    double getValueDifference(const Palette &palette) const noexcept;
    int isDiscriminative(const Palette &palette, const double &global_dist) const noexcept;

private:
    double cosinSimilarity(int index1, int index2) const noexcept;
};
