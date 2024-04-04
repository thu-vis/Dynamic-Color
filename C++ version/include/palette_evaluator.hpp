#pragma once

#include "color.hpp"
#include <map>

map<string, double> evaluate_palette(const vector<vector<double>> &in_palette);

map<string, double> evaluate_palette_grid(const vector<vector<double>> &in_palette,
                                          const int width, const int height, const vector<int> &data, const vector<int> &labels, const vector<int> &label_list, vector<vector<double>> &similarity);