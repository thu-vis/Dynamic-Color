/*
 *  grid_data.cpp
 *  Implementation of data driven for grid layout data.
 */

#include "data.hpp"
#include "utils.hpp"
#include "colorpair.hpp"
#include <unordered_set>
#include <array>

void GridData::initGridData(int width, int height, const vector<int> &data, const vector<int> &labels, const vector<int> &filter_labels,
                            const vector<vector<double>> &label_similarity, int matrix_start, int use_mode)
{
    // basic info and norm labels
    this->width = width;
    this->height = height;
    unordered_set<int> label_set(filter_labels.begin(), filter_labels.end());
    unordered_map<array<int, 2>, int, PosHash, PosEqual> pos2lpos;
    unordered_map<int, int> labelnorm;
    vector<double> area(filter_labels.size(), 0);
    double num_sum = 0.0;
    label_centers.clear();
    for (int i = 0; i < filter_labels.size(); i++)
    {
        labelnorm[filter_labels[i]] = i;
        label_centers.push_back({0.0, 0.0, 0.0});
    }
    int cur_pos = 0;
    const int max_id = labels.size() - 1;
    for (auto i = 0; i < width; i++)
    {
        for (auto j = 0; j < height; j++)
        {
            int id = i * height + j;
            if (data[id] > max_id)
                continue;
            int num = data[id];
            int label = labels[num];
            if (label_set.find(label) == label_set.end())
                continue;
            area[labelnorm[label]] += 1.0;
            grids.push_back(GridCell(i, j, id, num, cur_pos, label, labelnorm[label]));
            pos2lpos[{i, j}] = cur_pos;
            cur_pos += 1;
        }
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
        for (auto &grid : grids)
        {
            label_centers[grid.normlabel][0] += grid.x;
            label_centers[grid.normlabel][1] += grid.y;
            label_centers[grid.normlabel][2] += 1.0;
            center[0] += grid.x;
            center[1] += grid.y;
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
        for (auto &grid : grids)
        {
            int x = grid.x, y = grid.y;
            if (pos2lpos.find({x - 1, y}) != pos2lpos.end())
                grid.left = pos2lpos.at({x - 1, y});
            if (pos2lpos.find({x + 1, y}) != pos2lpos.end())
                grid.right = pos2lpos.at({x + 1, y});
            if (pos2lpos.find({x, y - 1}) != pos2lpos.end())
                grid.top = pos2lpos.at({x, y - 1});
            if (pos2lpos.find({x, y + 1}) != pos2lpos.end())
                grid.bottom = pos2lpos.at({x, y + 1});
            if (pos2lpos.find({x - 1, y - 1}) != pos2lpos.end())
                grid.lefttop = pos2lpos.at({x - 1, y - 1});
            if (pos2lpos.find({x - 1, y + 1}) != pos2lpos.end())
                grid.leftbottom = pos2lpos.at({x - 1, y + 1});
            if (pos2lpos.find({x + 1, y - 1}) != pos2lpos.end())
                grid.righttop = pos2lpos.at({x + 1, y - 1});
            if (pos2lpos.find({x + 1, y + 1}) != pos2lpos.end())
                grid.rightbottom = pos2lpos.at({x + 1, y + 1});
        }

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
        // int num = 0;
        for (auto &grid : grids)
        {
            int n = 0;
            int label_id = grid.normlabel;
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
            if (grid.top != -1)
                judgeLabel(grids[grid.top].normlabel, 1);
            if (grid.bottom != -1)
                judgeLabel(grids[grid.bottom].normlabel, 1);
            if (grid.left != -1)
                judgeLabel(grids[grid.left].normlabel, 1);
            if (grid.right != -1)
                judgeLabel(grids[grid.right].normlabel, 1);
            if (grid.lefttop != -1)
                judgeLabel(grids[grid.lefttop].normlabel, sqrt1_2);
            if (grid.righttop != -1)
                judgeLabel(grids[grid.righttop].normlabel, sqrt1_2);
            if (grid.leftbottom != -1)
                judgeLabel(grids[grid.leftbottom].normlabel, sqrt1_2);
            if (grid.rightbottom != -1)
                judgeLabel(grids[grid.rightbottom].normlabel, sqrt1_2);

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

double GridData::calDDValue(const Palette &palette, int start, int end)
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

double GridData::calHValue(const Palette &palette, int start, int end)
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

vector<int> GridData::getRandomWheelSequence()
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