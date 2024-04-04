/*
 *  data.hpp
 *  Header file for data driven evaluation
 */

#pragma once
#include "color.hpp"
#include <iostream>
#include <unordered_map>
#include <cassert>
using namespace std;

inline double g(double similartiy)
{
    // return -0.5 * (1 + similartiy);
    return -similartiy;
}

class Data
{
public:
    enum DataType
    {
        scatter = 0,
        line = 1,
        grid = 2
    };
    enum DrivenMode
    {
        discrimination = 1,
        similarity = 2
    };
    DataType type = DataType::grid;
    DrivenMode mode = DrivenMode::discrimination;
    Data(const DataType type, const DrivenMode mode) : type(type), mode(mode) {}

    virtual double getDataDrivenValue(const Palette &palette, int start, int end, double data_driven = 1.0, double harmony = 1.0)
    {
        cout << "Error: getDataDrivenValue() is not implemented" << endl;
        return 0;
    }
    virtual vector<int> getLabelSequence()
    {
        cout << "Error: getLabelSequence() is not implemented" << endl;
        return vector<int>();
    }
    int label_num = 0;
};

class ScatterData : public Data
{
public:
    ScatterData(DrivenMode mode, int k = 7) : Data(DataType::scatter, mode)
    {
        this->k = k;
    }
    double getDataDrivenValue(const Palette &palette, int start, int end, double data_driven = 1.0, double harmony = 1.0) override
    {
        return data_driven * calDDValue(palette, start, end) + harmony * calHValue(palette, start, end);
    }
    void initScatterData(const vector<vector<double>> &points, const vector<int> &labels, const vector<int> &filter_labels,
                         const vector<vector<double>> &label_similarity, int matrix_start = 0, int use_mode = 0);
    vector<int> getLabelSequence() override
    {
        return getRandomWheelSequence();
    }

private:
    int k;
    vector<int> getRandomWheelSequence();
    vector<vector<double>> args;
    vector<vector<double>> hargs;
    double calDDValue(const Palette &palette, int start, int end);
    double calHValue(const Palette &palette, int start, int end);
    vector<array<double, 3>> label_centers;
    array<double, 3> center = {0.0, 0.0, 0.0};
    vector<double> label_alphas;
    vector<int> arg_sort_alpha;
};

class LineData : public Data
{
public:
    LineData(DrivenMode mode, int k = 7) : Data(DataType::scatter, mode)
    {
        this->k = k;
    }
    double getDataDrivenValue(const Palette &palette, int start, int end, double data_driven = 1.0, double harmony = 1.0) override
    {
        return data_driven * calDDValue(palette, start, end) + harmony * calHValue(palette, start, end);
    }
    void initLineData(const vector<vector<double>> &points, const vector<int> &labels, const vector<int> &filter_labels,
                      const vector<vector<double>> &label_similarity, int matrix_start = 0, int use_mode = 0);
    vector<int> getLabelSequence() override
    {
        return getRandomWheelSequence();
    }

private:
    int k;
    vector<int> getRandomWheelSequence();
    vector<vector<double>> args;
    vector<vector<double>> hargs;
    double calDDValue(const Palette &palette, int start, int end);
    double calHValue(const Palette &palette, int start, int end);
    vector<vector<double>> label_centers;
    vector<double> label_heights;
    vector<int> arg_sort_heights;
};

class GridCell
{
public:
    int x, y, id;
    int num, pos;
    int label, normlabel;
    int left = -1, right = -1, top = -1, bottom = -1;
    int lefttop = -1, righttop = -1, leftbottom = -1, rightbottom = -1;
    GridCell(int x, int y, int id, int num, int pos, int label, int normlabel) : x(x), y(y), id(id), num(num), pos(pos), label(label), normlabel(normlabel) {}
};

class GridData : public Data
{
public:
    GridData(DrivenMode mode) : Data(DataType::grid, mode) {}
    double getDataDrivenValue(const Palette &palette, int start, int end, double data_driven = 1.0, double harmony = 1.0) override
    {
        return data_driven * calDDValue(palette, start, end) + harmony * calHValue(palette, start, end);
    }
    void initGridData(int width, int height, const vector<int> &data, const vector<int> &labels, const vector<int> &filter_labels,
                      const vector<vector<double>> &label_similarity, int matrix_start = 0, int use_mode = 0);
    vector<int> getLabelSequence() override
    {
        return getRandomWheelSequence();
    }

private:
    vector<int> getRandomWheelSequence();
    int width, height;
    vector<GridCell> grids;
    vector<vector<double>> args;
    vector<vector<double>> hargs;
    double calDDValue(const Palette &palette, int start, int end);
    double calHValue(const Palette &palette, int start, int end);
    vector<array<double, 3>> label_centers;
    array<double, 3> center = {0.0, 0.0, 0.0};
    vector<double> label_alphas;
    vector<int> arg_sort_alpha;
};