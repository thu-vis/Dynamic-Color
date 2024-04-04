/*
 *  mtl_adjuster.hpp
 *  Header file for multitask learning adjust args
 */

#pragma once
#include <vector>
#include <cmath>
using namespace std;

class ArgsAdjuster
{
public:
    int arg_nums = 0;
    int T = 0;
    int point = 0;
    vector<vector<double>> loss;
    vector<vector<double>> args;
    vector<vector<double>> checkpoints;

    ArgsAdjuster(int arg_nums) : arg_nums(arg_nums)
    {
        loss.push_back(vector<double>(arg_nums, 0));
        args.push_back(vector<double>(arg_nums, 1));
        checkpoints.push_back(vector<double>(arg_nums, 0));
    }

    double getValueBf(vector<double> &lossvalues) const noexcept
    {
        double value = 0;
        for (int i = 0; i < arg_nums; i++)
        {
            value += lossvalues[i] * args[T - 1][i];
        }
        return value;
    }

    double getValueCur(vector<double> &lossvalues) const noexcept
    {
        double value = 0;
        for (int i = 0; i < arg_nums; i++)
        {
            value += lossvalues[i] * args[T][i];
        }
        return value;
    }

    void adjust(vector<double> &t1values) noexcept
    {
        auto t0values = loss[T];
        T += 1;
        vector<double> w(arg_nums, 0);
        double tmpw = 0;
        double sumw = 0;
        for (int i = 0; i < arg_nums; i++)
        {
            tmpw = t0values[i] == 0 ? T : t1values[i] / t0values[i];
            w[i] = exp(tmpw / T);
            sumw += w[i];
        }
        for (int i = 0; i < arg_nums; i++)
        {
            w[i] = w[i] / sumw * arg_nums;
        }
        loss.push_back(t1values);
        args.push_back(w);
    }

    int findMinIndex() noexcept
    {
        checkpoints.push_back(loss[T]);
        point += 1;
        int min_index = 0;
        double min_value = checkpoints[point][0] - checkpoints[point - 1][0];
        for (int i = 1; i < arg_nums; i++)
        {
            double tp_value = checkpoints[point][i] - checkpoints[point - 1][i];
            if (tp_value < min_value)
            {
                min_index = i;
                min_value = tp_value;
            }
        }
        return min_index;
    }
};