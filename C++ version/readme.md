# Dynamic-Color-C++/Python

This is a pybind11 project, which uses C++ to generate libraries usable in Python. Of course, you can also refer to the invocation method in dynamic_color.cpp to rewrite a main.cpp and compile it. You can modify the cmake file to directly use it in C++.

### Code dependency

This project utilizes cmake to generate pybind11 programs, which can be directly called in Python for C code execution; please clone the code from the corresponding folder of pybind11: https://github.com/pybind/pybind11

Environment requirements: g++ cmake
Open-source library dependencies used in the code: eigen, color-util, nanoflann, nlohmann

Most of these dependencies are already included in the include folder (in source code form). If you wish to install these dependencies in other forms, please refer to:

- eigen: https://eigen.tuxfamily.org/index.php?title=Main_Page#Documentation for vector calculations
- color-util: https://github.com/yuki-koyama/color-util for color systems and color difference calculation // Fixed a small error and raised an issue to the original project
- nanoflann: https://github.com/jlblancoc/nanoflann for fast nearest neighbor search using kdtree
- nlohmann/json: https://github.com/nlohmann/json for JSON reading and storage using C++

### Code start

```
mkdir build
cd build
cmake ..
make
```

Then you will find the packaged library in the build directory, for example, under python version 3.8 in linux as dynamic-color.cpython-38-x86_64-linux-gnu.so

You can use it like a library in Python under version 3.8.

**Very important: Please place the "static" folder in the root directory of the calling environment to generate results quickly.**
You can download it from Google Drive:
https://drive.google.com/drive/folders/1XtNUwa650u4l0qPGOCvZ6DbmdTABvz8M?usp=drive_link.


### Code use
We provide four functions: findPalette, findPaletteForScatter, findPaletteForLine, and findPaletteForGrid, which respectively generate colorings for scenarios including non-distribution, scatter plots, PCP, and grid visualizations. Specific parameters can be viewed using the help function in Python or by referring to src/dynamic_color.cpp(as follow).
```cpp
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
```
