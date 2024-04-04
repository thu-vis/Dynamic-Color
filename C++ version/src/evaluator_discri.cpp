/*
 *  evaluator_discri.cpp
 *  Implementation of color evaluation of discrimination.
 */

#include "evaluator_discri.hpp"
#include "matrix.hpp"
#include <nlohmann/json.hpp>
#include <fstream>
#include <cmath>
#include <omp.h>
using json = nlohmann::json;
using namespace std;

static bool init = false;
static PointCloud cloud;
static shared_ptr<KDTree> kdtree;
static vector<vector<float>> cosin_difference;
static vector<string> c3names;
static unordered_map<int, int> T;
static int C, W;

float cosinSimilarity(int index1, int index2)
{
    double sa = 0, sb = 0, sc = 0;
    double ta, tb;
    int id1, id2;
    for (auto w = 0; w < W; w++)
    {
        id1 = index1 * W + w;
        id2 = index2 * W + w;
        if (T.find(id1) == T.end())
            ta = 0;
        else
            ta = T.at(id1);
        if (T.find(id2) == T.end())
            tb = 0;
        else
            tb = T.at(id2);
        sa += ta * ta;
        sb += tb * tb;
        sc += ta * tb;
    }
    return sc / sqrt(sa * sb);
}

void loadDifferenceInfo()
{
    init = true;
    ifstream c3file("static/c3_data.json");
    json c3data = json::parse(c3file);
    auto color = c3data["color"];
    for (auto i = 0; i < color.size(); i += 3)
    {
        vector<double> c3color = {color[i], color[i + 1], color[i + 2]};
        cloud.points.push_back(c3color);
    }
    kdtree = make_shared<KDTree>(3, cloud, KDTreeSingleIndexAdaptorParams(10));
    kdtree->buildIndex();
    cout << cloud.points.size() << endl;

    auto names = c3data["terms"];
    for (auto i = 0; i < names.size(); i++)
    {
        c3names.push_back(names[i]);
    }
    auto Tjson = c3data["T"];
    for (auto i = 0; i < Tjson.size(); i += 2)
    {
        T[Tjson[i]] = Tjson[i + 1];
    }
    C = cloud.points.size();
    W = c3names.size();

    cosin_difference.resize(C);
    for (auto i = 0; i < C; i++)
    {
        cosin_difference[i].resize(C);
    }

    auto start2 = clock();
    cout << "start init" << C << endl;

    loadBinaryMatrix(cosin_difference, "static/cosin_difference.bin");
    auto end2 = clock();
    cout << "init time: " << (double)(end2 - start2) / CLOCKS_PER_SEC << "s" << cosin_difference[6125][5397] << " " << cosinSimilarity(6125, 5397) << endl;
}

shared_ptr<KDTree> get_kdtree()
{
    if (!init)
    {
        loadDifferenceInfo();
    }

    return kdtree;
}

const vector<vector<float>> &get_cosine_diff()
{
    if (!init)
    {
        loadDifferenceInfo();
    }

    return cosin_difference;
}

DiscriEvaluator::DiscriEvaluator()
{
    // init kdtree
    if (!init)
    {
        loadDifferenceInfo();
    }
}

DiscriEvaluator::~DiscriEvaluator()
{
}

double DiscriEvaluator::evaluate(const Palette &palette, double alpha, double beta, double &min_dist) const noexcept
{
    min_dist = getValueDifference(palette);
    return alpha * getNameDifference(palette) + beta * min_dist;
}

double DiscriEvaluator::getNameDifference(const Palette &palette) const
{
    vector<int> color_indexes;
    vector<unsigned int> indices = {0, 0};
    vector<double> kddists = {0.0, 0.0};
    for (auto &color : palette)
    {
        vector<double> lab = {color.lab_color[0], color.lab_color[1], color.lab_color[2]};
        kdtree->knnSearch(&lab[0], 1, &indices[0], &kddists[0]);
        color_indexes.push_back(indices[0]);
    }

    double total_sum = 0;
    for (auto i = 0; i < color_indexes.size(); i++)
    {
        double sum = 0;
        for (auto j = i + 1; j < color_indexes.size(); j++)
        {
            sum += cosin_difference[color_indexes[i]][color_indexes[j]];
        }
        total_sum += sum;
    }
    total_sum /= (palette.size() * (palette.size() - 1) * 0.5);
    return total_sum;
}

double DiscriEvaluator::getValueDifference(const Palette &palette) const noexcept
{
    double min_dist = 100000.0;
    for (auto i = 0; i < palette.size(); i++)
    {
        double min_dist_private = min_dist;
        for (auto j = i + 1; j < palette.size(); j++)
        {
            double tp = dist(palette[i], palette[j]);
            if (tp < min_dist_private)
                min_dist_private = tp;
        }

        if (min_dist_private < min_dist)
            min_dist = min_dist_private;
    }
    return min_dist;
}

double DiscriEvaluator::cosinSimilarity(int index1, int index2) const noexcept
{
    return cosinSimilarity(index1, index2);
}

int DiscriEvaluator::isDiscriminative(const Palette &palette, const double &global_dist) const noexcept
{
    int idx = -1;
    for (auto i = 0; i < palette.size(); i++)
    {
        for (auto j = i + 1; j < palette.size(); j++)
        {
            double color_dis = dist(palette[i], palette[j]);
            if (color_dis < global_dist)
            {
                return j;
            }
        }
    }
    return idx;
}