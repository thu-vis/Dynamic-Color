#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include "generator.hpp"
#include "palette_evaluator.hpp"
#include <ctime>
#include <memory>
#include <omp.h>
#include <vector>
using namespace std;
namespace py = pybind11;

void preloadFunction()
{
    loadDifferenceInfo();
    PyRun_SimpleString("print('Preloading content...')");
}

vector<vector<double>> findPaletteGlobal(vector<vector<double>> &rgb_colors, vector<int> &palette_sizes,
                                         vector<double> &hue_range, const double tradeoff = 1.0, const int iter_num = 16, const double global_dec = 0.99, const int thread = 12)
{
    // set openmp
    omp_set_num_threads(thread);

    // convert inputs to Palette
    resetSeed();
    Generator generator;
    generator.scope.hue_scope[0] = hue_range[0];
    generator.scope.hue_scope[1] = hue_range[1];
    generator.evaluator.A = 1.0;
    generator.tradeoff = tradeoff;
    if (hue_range[0] == 0 && hue_range[1] == 360)
    {
        generator.scope.all_hue = true;
    }
    else
    {
        generator.scope.all_hue = false;
    }
    generator.iter_num = iter_num;
    generator.global_dec = global_dec;

    Palette palette;
    for (auto &color : rgb_colors)
    {
        palette.push_back(Color(rgb(color[0], color[1], color[2])));
    }
    vector<shared_ptr<Data>> driven_data;
    Palette basic_palette = generator.findBasicColors(palette, palette_sizes, driven_data);
    Palette new_palette = generator.findPaletteGlobal(basic_palette, palette_sizes, driven_data);

    vector<vector<double>> result;
    for (auto &color : basic_palette)
    {
        result.push_back({color.rgb_color[0], color.rgb_color[1], color.rgb_color[2]});
    }
    for (auto &color : new_palette)
    {
        result.push_back({color.rgb_color[0], color.rgb_color[1], color.rgb_color[2]});
    }
    return result;
}

vector<vector<double>> findPaletteGlobalScatter(vector<vector<double>> &rgb_colors, vector<int> &palette_sizes,
                                                vector<vector<double>> &data, vector<int> &labels, vector<int> &label_list, vector<vector<double>> &similarity,
                                                vector<double> &hue_range, int mode, const double tradeoff = 1.0, const int iter_num = 48, const double global_dec = 0.99, const int thread = 12)
{
    // --------------------------------------init setting---------------------------------------------------------------------------------------------
    omp_set_num_threads(thread);
    const double A = 1.0;
    resetSeed();
    Generator generator;
    generator.scope.hue_scope[0] = hue_range[0];
    generator.scope.hue_scope[1] = hue_range[1];
    generator.evaluator.A = A;
    generator.tradeoff = tradeoff;
    if (hue_range[0] == 0 && hue_range[1] == 360)
    {
        generator.scope.all_hue = true;
    }
    else
    {
        generator.scope.all_hue = false;
    }
    generator.iter_num = iter_num;
    generator.global_dec = global_dec;

    Palette palette;
    for (auto &color : rgb_colors)
    {
        palette.push_back(Color(rgb(color[0], color[1], color[2])));
    }
    int start = 0, end = 0;
    auto data_mode = Data::DrivenMode::discrimination;
    if (mode == 1)
    {
        data_mode = Data::DrivenMode::similarity;
    }

    // --------------------------------------convert cluster data info and generate basic color----------------------------------------------------------
    // convert cluster average data
    unordered_map<int, int> class_num;
    for (auto &label : labels)
    {
        class_num[label]++;
    }

    start = 0;
    vector<int> clabel_list;
    vector<int> cluster_num(palette_sizes.size(), 0);
    unordered_map<int, int> class_2_cluster;
    for (int i = 0; i < palette_sizes.size(); i++)
    {
        for (int j = 0; j < palette_sizes[i]; j++)
        {
            cluster_num[i] += class_num[label_list[start + j]];
            class_2_cluster[label_list[start + j]] = i;
        }
        clabel_list.push_back(i);
        start += palette_sizes[i];
    }

    vector<vector<double>> cluster_similarity(palette_sizes.size(), vector<double>(palette_sizes.size(), 0.0));
    for (int i = 0; i < label_list.size(); i++)
    {
        for (int j = 0; j < label_list.size(); j++)
        {
            int c1 = class_2_cluster[label_list[i]];
            int c2 = class_2_cluster[label_list[j]];
            int num1 = class_num[label_list[i]];
            int num2 = class_num[label_list[j]];
            cluster_similarity[c1][c2] += similarity[i][j] * num1 * num2;
        }
    }
    for (int i = 0; i < palette_sizes.size(); i++)
    {
        for (int j = 0; j < palette_sizes.size(); j++)
        {
            cluster_similarity[i][j] /= cluster_num[i] * cluster_num[j];
        }
    }

    vector<int> clabels;
    for (int i = 0; i < labels.size(); i++)
    {
        clabels.push_back(class_2_cluster[labels[i]]);
    }
    vector<shared_ptr<Data>> cluster_driven_data;
    for (int i = 0; i < palette_sizes.size(); i++)
    {
        ScatterData data1(data_mode);
        vector<int> filter_clist = {i};
        data1.initScatterData(data, clabels, filter_clist, cluster_similarity, i);
        cluster_driven_data.push_back(make_shared<ScatterData>(data1));
    }
    ScatterData clusterdata(data_mode);
    clusterdata.initScatterData(data, clabels, clabel_list, cluster_similarity, 0, 2);
    cluster_driven_data.push_back(make_shared<ScatterData>(clusterdata));

    // init basic colors
    Palette basic_palette = generator.findBasicColors(palette, palette_sizes, cluster_driven_data);

    // --------------------------------------convert driven data info and generate final color----------------------------------------------------------
    // convert driven data
    vector<shared_ptr<Data>>
        driven_data;
    start = 0;
    end = 0;
    for (int i = 0; i < palette_sizes.size(); i++)
    {
        end += palette_sizes[i];
        vector<int> filter_labels(label_list.begin() + start, label_list.begin() + end);

        ScatterData data1(data_mode);
        data1.initScatterData(data, labels, filter_labels, similarity, start);
        driven_data.push_back(make_shared<ScatterData>(data1));
        start = end;
    }
    ScatterData totaldata(data_mode);
    totaldata.initScatterData(data, labels, label_list, similarity, 0, 1);
    driven_data.push_back(make_shared<ScatterData>(totaldata));

    Palette new_palette = generator.findPaletteGlobal(basic_palette, palette_sizes, driven_data);

    // --------------------------------------get output、clear data and debug logs in best_value_save----------------------------------------------------------
    // convert new_palette to vector<vector<double>>
    vector<vector<double>> result;
    for (auto &color : basic_palette)
    {
        result.push_back({color.rgb_color[0], color.rgb_color[1], color.rgb_color[2]});
    }
    for (auto &color : new_palette)
    {
        result.push_back({color.rgb_color[0], color.rgb_color[1], color.rgb_color[2]});
    }

    // clear driven_data
    for (auto &data : cluster_driven_data)
    {
        data.reset();
    }
    for (auto &data : driven_data)
    {
        data.reset();
    }
    cluster_driven_data.clear();
    driven_data.clear();
    return result;
}

vector<vector<double>> findPaletteGlobalLine(vector<vector<double>> &rgb_colors, vector<int> &palette_sizes,
                                             vector<vector<double>> &data, vector<int> &labels, vector<int> &label_list, vector<vector<double>> &similarity,
                                             vector<double> &hue_range, int mode, const double tradeoff = 1.0, const int iter_num = 48, const double global_dec = 0.99, const int thread = 12)
{
    // --------------------------------------init setting---------------------------------------------------------------------------------------------
    omp_set_num_threads(thread);
    const double A = 1.0; // penalty factor in Ldiscrimination
    resetSeed();
    Generator generator;
    generator.scope.hue_scope[0] = hue_range[0];
    generator.scope.hue_scope[1] = hue_range[1];
    generator.evaluator.A = A;
    generator.tradeoff = tradeoff;
    if (hue_range[0] == 0 && hue_range[1] == 360)
    {
        generator.scope.all_hue = true;
    }
    else
    {
        generator.scope.all_hue = false;
    }
    generator.iter_num = iter_num;
    generator.global_dec = global_dec;

    Palette palette;
    for (auto &color : rgb_colors)
    {
        palette.push_back(Color(rgb(color[0], color[1], color[2])));
    }
    int start = 0, end = 0;
    auto data_mode = Data::DrivenMode::discrimination;
    if (mode == 1)
    {
        data_mode = Data::DrivenMode::similarity;
    }

    // --------------------------------------convert cluster data info and generate basic color----------------------------------------------------------
    // convert cluster average data
    unordered_map<int, int> class_num;
    for (auto &label : labels)
    {
        class_num[label]++;
    }

    start = 0;
    vector<int> clabel_list;
    vector<int> cluster_num(palette_sizes.size(), 0);
    unordered_map<int, int> class_2_cluster;
    for (int i = 0; i < palette_sizes.size(); i++)
    {
        for (int j = 0; j < palette_sizes[i]; j++)
        {
            cluster_num[i] += class_num[label_list[start + j]];
            class_2_cluster[label_list[start + j]] = i;
        }
        clabel_list.push_back(i);
        start += palette_sizes[i];
    }

    vector<vector<double>> cluster_similarity(palette_sizes.size(), vector<double>(palette_sizes.size(), 0.0));
    for (int i = 0; i < label_list.size(); i++)
    {
        for (int j = 0; j < label_list.size(); j++)
        {
            int c1 = class_2_cluster[label_list[i]];
            int c2 = class_2_cluster[label_list[j]];
            int num1 = class_num[label_list[i]];
            int num2 = class_num[label_list[j]];
            cluster_similarity[c1][c2] += similarity[i][j] * num1 * num2;
        }
    }
    for (int i = 0; i < palette_sizes.size(); i++)
    {
        for (int j = 0; j < palette_sizes.size(); j++)
        {
            cluster_similarity[i][j] /= cluster_num[i] * cluster_num[j];
        }
    }

    vector<int> clabels;
    for (int i = 0; i < labels.size(); i++)
    {
        clabels.push_back(class_2_cluster[labels[i]]);
    }
    vector<shared_ptr<Data>> cluster_driven_data;
    for (int i = 0; i < palette_sizes.size(); i++)
    {
        LineData data1(data_mode);
        vector<int> filter_clist = {i};
        data1.initLineData(data, clabels, filter_clist, cluster_similarity, i);
        cluster_driven_data.push_back(make_shared<LineData>(data1));
    }
    LineData clusterdata(data_mode);
    clusterdata.initLineData(data, clabels, clabel_list, cluster_similarity, 0, 2);
    cluster_driven_data.push_back(make_shared<LineData>(clusterdata));

    // init basic colors
    Palette basic_palette = generator.findBasicColors(palette, palette_sizes, cluster_driven_data);

    // --------------------------------------convert driven data info and generate final color----------------------------------------------------------
    // convert driven data
    vector<shared_ptr<Data>> driven_data;
    start = 0;
    end = 0;
    for (int i = 0; i < palette_sizes.size(); i++)
    {
        end += palette_sizes[i];
        vector<int> filter_labels(label_list.begin() + start, label_list.begin() + end);

        LineData data1(data_mode);
        data1.initLineData(data, labels, filter_labels, similarity, start);
        driven_data.push_back(make_shared<LineData>(data1));
        start = end;
    }
    LineData totaldata(data_mode);
    totaldata.initLineData(data, labels, label_list, similarity, 0, 1);
    driven_data.push_back(make_shared<LineData>(totaldata));

    Palette new_palette = generator.findPaletteGlobal(basic_palette, palette_sizes, driven_data);

    // --------------------------------------get output、clear data and debug logs in best_value_save----------------------------------------------------------
    // convert new_palette to vector<vector<double>>
    vector<vector<double>> result;
    for (auto &color : basic_palette)
    {
        result.push_back({color.rgb_color[0], color.rgb_color[1], color.rgb_color[2]});
    }
    for (auto &color : new_palette)
    {
        result.push_back({color.rgb_color[0], color.rgb_color[1], color.rgb_color[2]});
    }

    // clear driven_data
    for (auto &data : cluster_driven_data)
    {
        data.reset();
    }
    for (auto &data : driven_data)
    {
        data.reset();
    }
    cluster_driven_data.clear();
    driven_data.clear();
    return result;
}

vector<vector<double>> findPaletteGlobalGrid(vector<vector<double>> &rgb_colors, vector<int> &palette_sizes,
                                             int width, int height, vector<int> &data, vector<int> &labels, vector<int> &label_list, vector<vector<double>> &similarity,
                                             vector<double> &hue_range, int mode, const double tradeoff = 1.0, const int iter_num = 48, const double global_dec = 0.99, const int thread = 12)
{
    // --------------------------------------init setting---------------------------------------------------------------------------------------------
    omp_set_num_threads(thread);
    const double A = 1.0; // penalty factor in Ldiscrimination
    resetSeed();
    Generator generator;
    generator.scope.hue_scope[0] = hue_range[0];
    generator.scope.hue_scope[1] = hue_range[1];
    generator.evaluator.A = A;
    if (hue_range[0] == 0 && hue_range[1] == 360)
    {
        generator.scope.all_hue = true;
    }
    else
    {
        generator.scope.all_hue = false;
    }
    generator.iter_num = iter_num;
    generator.global_dec = global_dec;
    generator.tradeoff = tradeoff;

    Palette palette;
    for (auto &color : rgb_colors)
    {
        palette.push_back(Color(rgb(color[0], color[1], color[2])));
    }
    int start = 0, end = 0;
    auto data_mode = Data::DrivenMode::discrimination;
    if (mode == 1)
    {
        data_mode = Data::DrivenMode::similarity;
    }

    // --------------------------------------convert cluster data info and generate basic color----------------------------------------------------------
    // convert cluster average data
    unordered_map<int, int> class_num;
    for (auto &label : labels)
    {
        class_num[label]++;
    }

    start = 0;
    vector<int> clabel_list;
    vector<int> cluster_num(palette_sizes.size(), 0);
    unordered_map<int, int> class_2_cluster;
    for (int i = 0; i < palette_sizes.size(); i++)
    {
        for (int j = 0; j < palette_sizes[i]; j++)
        {
            cluster_num[i] += class_num[label_list[start + j]];
            class_2_cluster[label_list[start + j]] = i;
        }
        clabel_list.push_back(i);
        start += palette_sizes[i];
    }

    vector<vector<double>> cluster_similarity(palette_sizes.size(), vector<double>(palette_sizes.size(), 0.0));
    for (int i = 0; i < label_list.size(); i++)
    {
        for (int j = 0; j < label_list.size(); j++)
        {
            int c1 = class_2_cluster[label_list[i]];
            int c2 = class_2_cluster[label_list[j]];
            int num1 = class_num[label_list[i]];
            int num2 = class_num[label_list[j]];
            cluster_similarity[c1][c2] += similarity[i][j] * num1 * num2;
        }
    }
    for (int i = 0; i < palette_sizes.size(); i++)
    {
        for (int j = 0; j < palette_sizes.size(); j++)
        {
            cluster_similarity[i][j] /= cluster_num[i] * cluster_num[j];
        }
    }

    vector<int> clabels;
    for (int i = 0; i < labels.size(); i++)
    {
        clabels.push_back(class_2_cluster[labels[i]]);
    }
    vector<shared_ptr<Data>> cluster_driven_data;
    for (int i = 0; i < palette_sizes.size(); i++)
    {
        GridData data1(data_mode);
        vector<int> filter_clist = {i};
        data1.initGridData(width, height, data, clabels, filter_clist, cluster_similarity, i);
        cluster_driven_data.push_back(make_shared<GridData>(data1));
    }
    GridData clusterdata(data_mode);
    clusterdata.initGridData(width, height, data, clabels, clabel_list, cluster_similarity, 0, 2);
    cluster_driven_data.push_back(make_shared<GridData>(clusterdata));

    // init basic colors
    Palette basic_palette = generator.findBasicColors(palette, palette_sizes, cluster_driven_data);

    // --------------------------------------convert driven data info and generate final color----------------------------------------------------------
    // convert driven data
    vector<shared_ptr<Data>> driven_data;
    start = 0;
    end = 0;
    for (int i = 0; i < palette_sizes.size(); i++)
    {
        end += palette_sizes[i];
        vector<int> filter_labels(label_list.begin() + start, label_list.begin() + end);

        GridData data1(data_mode);
        data1.initGridData(width, height, data, labels, filter_labels, similarity, start);
        driven_data.push_back(make_shared<GridData>(data1));
        start = end;
    }
    GridData totaldata(data_mode);
    totaldata.initGridData(width, height, data, labels, label_list, similarity, 0, 1);
    driven_data.push_back(make_shared<GridData>(totaldata));

    Palette new_palette = generator.findPaletteGlobal(basic_palette, palette_sizes, driven_data);

    // --------------------------------------get output、clear data and debug logs in best_value_save----------------------------------------------------------
    // convert new_palette to vector<vector<double>>
    vector<vector<double>> result;
    for (auto &color : basic_palette)
    {
        result.push_back({color.rgb_color[0], color.rgb_color[1], color.rgb_color[2]});
    }
    for (auto &color : new_palette)
    {
        result.push_back({color.rgb_color[0], color.rgb_color[1], color.rgb_color[2]});
    }

    // clear driven_data
    for (auto &data : cluster_driven_data)
    {
        data.reset();
    }
    for (auto &data : driven_data)
    {
        data.reset();
    }
    cluster_driven_data.clear();
    driven_data.clear();
    return result;
}

PYBIND11_MODULE(dynamic_color, m)
{
    preloadFunction();
    m.doc() = "dynamic_color";

    py::list hue_range = py::list();
    hue_range.append(py::int_(0));
    hue_range.append(py::int_(360));

    m.def("findPalette", &findPaletteGlobal, "Find global palette",
          py::arg_v("rgb_colors", py::list(), "List of RGB colors in the range [0, 1]，[[0,0,1]] or [] when initializing"),
          py::arg_v("palette_sizes", py::list(), "List representing the number of children nodes under each parent node, such as [5]"),
          py::arg_v("hue_range", hue_range, "Hue range, default [0, 360]"),
          py::arg_v("tradeoff", 1.0, "Tradeoff between intra-cluster and inter-cluster distances, default 1.0"),
          py::arg_v("iter_num", 48, "Number of iterations for the SA internal loop, default 48"),
          py::arg_v("global_dec", 0.99, "Global decay rate for SA, default 0.99"),
          py::arg_v("thread", 12, "Number of threads for OpenMP parallelization, default 12"));

    m.def("findPaletteForScatter",
          &findPaletteGlobalScatter,
          "A global color is generated for scatter data.",
          py::arg_v("rgb_colors", py::list(), "List of RGB colors in the range [0, 1], [[]] when initializing"),
          py::arg_v("palette_sizes", py::list(), "List representing the number of children nodes under each parent node, such as [5]"),
          py::arg_v("data", py::list(), "Line list, such as [[0, 0], [1, 1], [2, 2]]"),
          py::arg_v("labels", py::list(), "Array recording the category of each sample."),
          py::arg_v("label_list", py::list(), "Array recording all labels in the same order as before"),
          py::arg_v("similarity", py::list(), "A matrix of size n*n, recording the similarity of each category"),
          py::arg_v("hue_range", hue_range, "Hue range, default [0, 360]"),
          py::arg_v("mode", 1, "Mode parameter, 0 for discrimination, 1 for similarity, default 1"),
          py::arg_v("tradeoff", 1.0, "Tradeoff between intra-cluster and inter-cluster distances, default 1.0"),
          py::arg_v("iter_num", 48, "Number of iterations for the SA internal loop, default 48"),
          py::arg_v("global_dec", 0.99, "Global decay rate for SA, default 0.99"),
          py::arg_v("thread", 12, "Number of threads for OpenMP parallelization, default 12"));

    m.def("findPaletteForline",
          &findPaletteGlobalLine,
          "A global color is generated for line data.",
          py::arg_v("rgb_colors", py::list(), "List of RGB colors in the range [0, 1], [[]] when initializing"),
          py::arg_v("palette_sizes", py::list(), "List representing the number of children nodes under each parent node, such as [5]"),
          py::arg_v("data", py::list(), "Line list, such as [[0, 0], [1, 1], [2, 2]]"),
          py::arg_v("labels", py::list(), "Array recording the category of each sample."),
          py::arg_v("label_list", py::list(), "Array recording all labels in the same order as before"),
          py::arg_v("similarity", py::list(), "A matrix of size n*n, recording the similarity of each category"),
          py::arg_v("hue_range", hue_range, "Hue range, default [0, 360]"),
          py::arg_v("mode", 1, "Mode parameter, 0 for discrimination, 1 for similarity, default 1"),
          py::arg_v("tradeoff", 1.0, "Tradeoff between intra-cluster and inter-cluster distances, default 1.0"),
          py::arg_v("iter_num", 48, "Number of iterations for the SA internal loop, default 48"),
          py::arg_v("global_dec", 0.99, "Global decay rate for SA, default 0.99"),
          py::arg_v("thread", 12, "Number of threads for OpenMP parallelization, default 12"));

    m.def("findPaletteForGrid",
          &findPaletteGlobalGrid,
          "Find global palette for grid",
          py::arg_v("rgb_colors", py::list(), "List of RGB colors in the range [0, 1], [[]] when initializing"),
          py::arg_v("palette_sizes", py::list(), "List representing the number of children nodes under each parent node, such as [5]"),
          py::arg_v("width", py::int_(), "Width of the grid"),
          py::arg_v("height", py::int_(), "Height of the grid"),
          py::arg_v("data", py::list(), "Array recording the index of each element placed in the grid"),
          py::arg_v("labels", py::list(), "Array recording the category of each index, if it exceeds then it represents an empty grid"),
          py::arg_v("label_list", py::list(), "Array recording all labels in the same order as before"),
          py::arg_v("similarity", py::list(), "A matrix of size n*n, recording the similarity of each category"),
          py::arg_v("hue_range", hue_range, "Hue range, default [0, 360]"),
          py::arg_v("mode", 1, "Mode parameter, 0 for discrimination, 1 for similarity, default 1"),
          py::arg_v("tradeoff", 1.0, "Tradeoff between intra-cluster and inter-cluster distances, default 1.0"),
          py::arg_v("iter_num", 48, "Number of iterations for the SA internal loop, default 48"),
          py::arg_v("global_dec", 0.99, "Global decay rate for SA, default 0.99"),
          py::arg_v("thread", 12, "Number of threads for OpenMP parallelization, default 12"));
}
