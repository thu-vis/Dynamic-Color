/*
 *  scatter_data.cpp
 *  Implementation of data driven for scatterplot data.
 *
 *  Created by Jiashu Chen. 2024-02-25
 */

#include "data.hpp"
#include "utils.hpp"
#include "kdtree.hpp"
#include "colorpair.hpp"
#include <unordered_set>
#include <array>
using namespace std;

void ScatterData::initScatterData(const vector<vector<double>> &points, const vector<int> &labels, const vector<int> &filter_labels,
                                  const vector<vector<double>> &label_similarity, int matrix_start, int use_mode)
{
    // basic info and norm labels
    unordered_set<int> label_set(filter_labels.begin(), filter_labels.end());
    unordered_map<int, int> labelnorm;
    vector<double> area(filter_labels.size(), 0);
    double num_sum = 0.0;
    label_centers.clear();
    for (int i = 0; i < filter_labels.size(); i++)
    {
        labelnorm[filter_labels[i]] = i;
        label_centers.push_back({0.0, 0.0, 0.0});
    }

    PointCloud cloud;
    vector<int> normlabels;
    for (auto i = 0; i < points.size(); i++)
    {
        int label = labels[i];
        if (label_set.find(label) == label_set.end())
            continue;
        area[labelnorm[label]] += 1.0;
        cloud.points.push_back(points[i]);
        normlabels.push_back(labelnorm[label]);
    }
    for (auto &a : area)
    {
        num_sum += a;
    }
    label_num = filter_labels.size();

    // calculate grid centers
    if (use_mode != 1)
    {
        center[0] = 0.0;
        center[1] = 0.0;
        center[2] = 0.0;

        for (auto i = 0; i < cloud.points.size(); i++)
        {
            auto &point = cloud.points[i];
            auto &normlabel = normlabels[i];
            label_centers[normlabel][0] += point[0];
            label_centers[normlabel][1] += point[1];
            label_centers[normlabel][2] += 1.0;
            center[0] += point[0];
            center[1] += point[1];
            center[2] += 1.0;
        }

        for (int i = 0; i < label_num; i++)
        {
            label_centers[i][0] /= label_centers[i][2];
            label_centers[i][1] /= label_centers[i][2];
            label_centers[i][2] /= num_sum;
        }
        center[0] /= center[2];
        center[1] /= center[2];
        for (int i = 0; i < label_num; i++)
        {
            label_alphas.push_back(atan2(label_centers[i][1] - center[1], label_centers[i][0] - center[0]));
        }
        arg_sort_alpha = argsort(label_alphas);
    }

    // calculate args
    if (use_mode != 0)
    {
        KDTree2 kdtree(2, cloud, KDTreeSingleIndexAdaptorParams(10));
        kdtree.buildIndex();

        args.resize(filter_labels.size());
        for (int i = 0; i < filter_labels.size(); i++)
        {
            args[i].resize(filter_labels.size(), 0);
        }
        hargs.resize(filter_labels.size());
        for (int i = 0; i < filter_labels.size(); i++)
        {
            hargs[i].resize(filter_labels.size(), 0);
        }

        // solve data driven args
        for (auto i = 0; i < cloud.points.size(); i++)
        {
            vector<double> point(cloud.points[i]);
            int n = 0;
            int label_id = normlabels[i];
            vector<int> res_labels;
            vector<double> res_dists;
            auto judgeLabel = [&res_labels, &res_dists, &label_id, &n](int normedlabel, double dist)
            {
                if (normedlabel != label_id)
                {
                    res_labels.push_back(normedlabel);
                    res_dists.push_back(dist);
                    n++;
                }
            };

            // get knn index and distance
            vector<unsigned int> knn_indices(k + 1, 0);
            vector<double> knn_dists(k + 1, 0.0);
            kdtree.knnSearch(&point[0], k + 1, &knn_indices[0], &knn_dists[0]);
            for (int i = 1; i < k + 1; i++)
            {
                judgeLabel(normlabels[knn_indices[i]], knn_dists[i]);
            }

            if (n > 0)
            {
                double delta = 1.0 / n;
                for (int i = 0; i < n; i++)
                {
                    auto add = delta * res_dists[i];
                    hargs[label_id][res_labels[i]] += add;
                    hargs[res_labels[i]][label_id] += add;
                }
            }
        }

        // get avg basic harmony args
        for (int i = 0; i < label_num; i++)
        {
            for (int j = i + 1; j < label_num; j++)
            {
                hargs[i][j] /= num_sum;
            }
        }

        // get data driven args
        if (mode == DrivenMode::discrimination)
        {
            for (int i = 0; i < filter_labels.size(); i++)
            {
                for (int j = i + 1; j < filter_labels.size(); j++)
                {
                    args[i][j] = hargs[i][j];
                }
            }
        }
        if (mode == DrivenMode::similarity)
        {
            for (int i = 0; i < filter_labels.size(); i++)
            {
                for (int j = i + 1; j < filter_labels.size(); j++)
                {
                    args[i][j] = hargs[i][j] * g(label_similarity[i + matrix_start][j + matrix_start]);
                }
            }
        }
    }
}

double ScatterData::calDDValue(const Palette &palette, int start, int end)
{
    double dis = 0;
    for (int i = 0; i < args.size(); i++)
    {
        for (int j = i + 1; j < args.size(); j++)
        {
            dis += args[i][j] * dist(palette[i + start], palette[j + start]);
        }
    }
    return dis;
}

double ScatterData::calHValue(const Palette &palette, int start, int end)
{
    double harmony = 0;
    for (int i = 0; i < args.size(); i++)
    {
        for (int j = i + 1; j < args.size(); j++)
        {
            harmony += hargs[i][j] * pairColorHarmony(palette[i + start], palette[j + start]);
        }
    }
    return harmony;
}

vector<int> ScatterData::getRandomWheelSequence()
{
    int rotate = randomInt(0, label_num - 1);
    vector<int> res;
    for (int i = 0; i < label_num; i++)
    {
        res.push_back(arg_sort_alpha[(i + rotate) % label_num]);
    }
    if (randomOne() < 0.5)
        reverse(res.begin(), res.end());
    return res;
}