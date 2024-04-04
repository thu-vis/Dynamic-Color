#pragma once
#include "data.hpp"
#include <memory>
using namespace std;

class ClusterInfo
{
    // cluster info for global generator
public:
    const int id, start, end, num;
    const Color basic_color;
    double ratio, range, initrange, hue_range;
    double hue_min = 0.0, hue_max = 0.0;
    shared_ptr<Data> data = nullptr;
    ClusterInfo(int id, int start, int end, int num, Color basic_color) : id(id), start(start), end(end), num(num), basic_color(basic_color)
    {
        ratio = initrange = range = 0.0;
    }
};

inline double clusterDataDrivenValue(const Palette &palette, const vector<ClusterInfo> &clusters, bool is_print = false, double data_driven = 0, double harmony = 1)
{
    int start = 0, end = 0;
    double driven_value = 0;
    if (is_print)
        cout << "cluster num: " << clusters.size() << " " << palette.size() << endl;
    for (auto &cluster : clusters)
    {
        end += cluster.num;
        if (is_print)
        {
            cout << start << " " << end << endl;
            cout << cluster.start << " " << cluster.end << " " << cluster.num << endl;
        }
        driven_value += cluster.data->getDataDrivenValue(palette, start, end, data_driven, harmony) * cluster.ratio;
        if (is_print)
            cout << "cluster " << cluster.num << " " << cluster.ratio << " " << cluster.data->getDataDrivenValue(palette, start, end) << endl;
        start = end;
    }
    return driven_value;
}